#' A function to create data for NEHA with discrete time EHA data
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param node A character string name of the variable that gives the node id.
#' @param time A character string name of the variable that gives the time, in integers.
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param cascade A character string name of the variable that gives the cascade id
#' @return A data frame in which, in addition to all of the variables in 'eha_data', there is one column for each directed dyad, named 'i_j', where 'i' and 'j' are node ids, in which the value indicates the number of time points before 'time' that 'i' experienced the event in the respective 'cascade'. If 'i' did not experience the cascade before 'time', the value is 0.
#' @export
data_neha_discrete <- function(eha_data,node,time,event,cascade,covariates){

  # find out if node is numeric
  n1c <- substr(as.character(eha_data[,node]),1,1)

  if(!all(is.element(tolower(n1c),letters))){
  	print("appending n_ to beginning of all node id's since at least one seems numeric")
  	eha_data[,node] <- paste("n_",eha_data[,node],sep="")
  }

  # extract node ids
  node <- eha_data[,node]

  unique.nodes <- unique(node)

  # make set of possible diffusion relationships
  node.pairs <- combinat::combn(unique.nodes,2)
  node.pairs <- sort(c(paste(node.pairs[1,],node.pairs[2,],sep="_"),paste(node.pairs[2,],node.pairs[1,],sep="_")))

  # create the matrix in which to hold the diffusion tie variables
  diff_var_mat <- matrix(0,length(node),length(node.pairs))

  # define column labels of diffusion tie variables
  colnames(diff_var_mat) <- node.pairs

  # add indicator variables to the diffusion tie variables matrix
  for(i in 1:nrow(eha_data)){
    # find the adopting node
    node_i <- node[i]
    # find the policy id
    cascade_i <- eha_data[i,cascade]
    # find the time
    time_i <- eha_data[i,time]
    # find previous events in the same cascade
    previous.events <- node[which( (eha_data[,time] < time_i) & ( eha_data[,cascade] == cascade_i) & (eha_data[,event]==1)  & (node != node_i))  ]
    # calculate time since previous events
    if(length(previous.events) > 0){
      times.since.events <- time_i-eha_data[which( (eha_data[,time] < time_i) & ( eha_data[,cascade] == cascade_i) & (eha_data[,event]==1) & (node != node_i) ),time  ]
      col.ids <- paste(previous.events,node_i,sep="_")
      diff_var_mat[i,col.ids] <- times.since.events
    }

  }
  network_betas <- numeric(ncol(diff_var_mat))

  edge_sum <- apply(diff_var_mat,2,sum)
  diff_var_mat <- diff_var_mat[,which(edge_sum > 0)]
  # throw out edges with exactly 0 correlation
  cor_e_x <- cor(cbind(eha_data[,event],diff_var_mat))[1,-1]
  diff_var_mat <- diff_var_mat[,which(cor_e_x != 0)]
  data.frame(eha_data,diff_var_mat,stringsAsFactors=F)
}


#' A function to simulate a single cascade from the discrete-time logistic specification NEHA data generating process.
#' @param x The [(number of nodes)x(number of observed time points)]x[number of covariates + 2] data frame of covariates to be used in simulating cascades. Note, each node must be at risk for the event until the maximum time. This function does not currently permit dyadic covariates that are themselves a function of past events (e.g., the number of geographically adjacent nodes that have experienced an event)
#' @param node A character string name of the variable in x that gives the node id
#' @param time A character string name of the variable in x that gives the time, in integers
#' @param beta A (ncol(x)-2)x1 matrix of regression coefficients with row names that match column names in x
#' @param gamma A (number of diffusion ties)x1 matrix with row names of the form "sending-node_receiving-node" such that sending and receiving node ids ## match elements of x[,node]. Elements should be non-negatave numeric values.
#' @param a A non-negative numeric value that models the exponential decay of sender influence. The effect of a previous event experienced by a diffusion tie sender on the log odds of an event at time t is gamma*exp(-exp(a)*(t-source_time)), where source_time is the time the sender experienced the event.
#' @return A data frame in discrete-time EHM format for a single cascade, in which there is at most one event occurrence for each node, indicated by the 'event' variable.
#' @export
simulate_neha_discrete <- function(x,node,time,beta,gamma,a=-8){
  times <- sort(unique(x[,time]))
  nodes <- sort(unique(x[,node]))
  x_linear_predictor <- as.matrix(x[,row.names(beta)])%*%beta
  event_times <- rep(NA,length(nodes))
  source <- do.call('rbind',strsplit(rownames(gamma),"_"))[,1]
  follower <- do.call('rbind',strsplit(rownames(gamma),"_"))[,2]
  event <- rep(NA,nrow(x))
  for(t in times){
    nodes_at_risk <- nodes[is.na(event_times)]
    for(n in nodes_at_risk){
      xlp_nt <- x_linear_predictor[which((x[,time]==t) & (x[,node]==n))]
      n_sources <- source[which(follower==n)]
      source_effects <- gamma[which(follower==n),1]
      source_times <- event_times[match(n_sources,nodes)]
      # a = .25
      # plot(0:20,exp(a*-(0:20)),ylim=c(0,1))
      time_effects <- source_effects*exp(-exp(a)*(t-source_times))
      time_effects <- time_effects[which(!is.na(time_effects))]
      linear_predictor <- xlp_nt
      if(length(time_effects)>0){
        linear_predictor <- linear_predictor + sum(time_effects)
      }
      pr_nt <- 1/(1+exp(-linear_predictor))
      event_nt <- 1*(pr_nt > runif(1))
      event[which((x[,time]==t) & (x[,node]==n))] <- 0
      if(event_nt==1){
        event_times[which(nodes==n)] <- t
        event[which((x[,time]==t) & (x[,node]==n))] <- 1
      }
    }

  }
  x <- data.frame(x,event,stringsAsFactors=F)
  data.frame(na.omit(x),stringsAsFactors=F)
}


estimate_a <- function(dv,edges){

  one_times <- c(unlist(edges[which(dv==1),]))
  one_times <- one_times[which(one_times>0)]
  min_time <- min(one_times)
  max_time <- max(one_times)

  #exp(-exp(a)*max_time)/exp(-exp(a)*min_time) = 0.5
  #-exp(a)*max_time+exp(a)*min_time = log(0.5)
  #exp(a) = log(0.5)/(min_time-max_time)
  a = log(log(0.9)/(min_time-max_time))

  for(i in 1:ncol(edges)){
    edges[,i] <- (edges[,i] > 0 )*exp(-exp(a)*edges[,i])
  }

  list(a,edges)

}


#' @import boot foreach doParallel parallel
update_a <- function(data_with_aest,covariates,edges_subset,event,old_a,edge_vars,ncore=2){

  multa <- c(1.1,1.05,0.95,0.9)

  cvglm_a <- function(a,data_with_aest,covariates,edges_subset,event){

    hol <- function(y,prob){
      sum(log(prob^y*(1-prob)^(1-y)))
    }

    effect_names <- c(covariates,edges_subset)
    data_with_aest[,edges_subset] <- data_with_aest[,edges_subset]^a

    form <- as.formula(paste(event,"~",paste(effect_names,collapse="+"),sep=""))

    glm.est <- glm(form,family=binomial,data=data_with_aest)

    set.seed(9202011)

    cv.glm(data_with_aest,glm.est,hol,10)$delta[2]

  }

  cl <- makeCluster(ncore)
  registerDoParallel(cl)

  cv_scores <- foreach(a =  multa,.packages=c("boot")) %dopar% {
    cvglm_a(a,data_with_aest,covariates,edges_subset,event)
  }

  stopCluster(cl)

  cv_scores <- c(unlist(cv_scores))

  besta <- multa[which.max(cv_scores)]

  data_with_aest[, edge_vars] <- data_with_aest[, edge_vars]^besta

  new_a <- old_a*besta

  list(new_a,data_with_aest)

}


#' A function to estimate NEHA parameters using greedy edge addition
#' @import foreach doParallel parallel
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param A character string name of the variable that gives the node id
#' @param time A character string name of the variable that gives the time, in integers
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param cascade A character string name of the variable that gives the cascade id
#' @param covariates character vector of covariate names to include in the neha, excluding the intercept.
#' #' @param ncore integer indicating the number of cores to use in parallel computation.
#' @export
neha_geta <- function(eha_data,node,time,event,cascade,covariates,ncore=2){

  data_for_neha <- data_neha_discrete(eha_data,node,time,event,cascade,covariates)

  edge_vars <- names(data_for_neha)[which(!is.element(names(data_for_neha),names(eha_data)))]

  a.estimate <- estimate_a(data_for_neha[,event],data_for_neha[,edge_vars])

  print("a estimation complete")

  data_for_neha[,edge_vars] <- a.estimate[[2]]

  # results
  list(a_est = a.estimate[[1]],data_for_neha)

}


#' A function to select edges and prepare data for NEHA estimation
#' @import doParallel foreach glmulti parallel
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param node A character string name of the variable that gives the node id
#' @param time A character string name of the variable that gives the time, in integers
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param cascade A character string name of the variable that gives the cascade id
#' @param covariates character vector of covariate names to include in the neha, excluding the intercept.
#' @param ncore an integer giving the number of cores to use in parallel computation.
#' @param negative logical, experimental indicating whether to go thorugh a phase of negative tie inference
#' @return A list with five elements.
#' \itemize{
#'   \item a_est - The estimated value of alpha, the edge effect decay parameter.
#'   \item edges - A character vector giving the names of the edges inferred 'sender_receiver'.
#'   \item data_for_neha - A dataframe that can be used to find the NEHA estimates using a function for logistic regression.
#'   \item combined_formula - NEHA formula to use if you want a single gamma estimate.
#'   \item separate_formula - NEHA formula to use if you want a separate gamma estimate for each edge.
#' }
#' @examples
#' library(neha)

#' \dontrun{
#' # Simulate data for NEHA
#' # basic data parameters
#' cascades <- 50
#' nodes <- 20
#' times <- 30
#' nties <- 25

#' # generate dataframe
#' time <- sort(rep(1:times,nodes))
#' node <- paste("n",as.character(rep(1:nodes,times)),sep="")
#' intercept <- rep(1,length(time))
#' covariate <- runif(length(time))-2
#' data_for_sim <- data.frame(time,node,intercept,covariate,stringsAsFactors=FALSE)
#'
#' # regression parameters
#' beta <- cbind(c(-2.5,.25))
#' rownames(beta) <- c("intercept","covariate")
#'
#' # generate network effects
#' possible_ties <- rbind(t(combn(1:nodes,2)),t(combn(1:nodes,2))[,c(2,1)])
#' possible_ties <- paste(paste("n",possible_ties[,1],sep=""),paste("n",possible_ties[,2],sep=""),sep="_")
#' ties <- sample(possible_ties,nties)
#' gamma <- cbind(rep(1.5,length(ties)))
#' rownames(gamma) <- ties
#'
#' # initiate simulated data object
#' simulated_data <- NULL
#'
#' # generate the data one cascade at a time
#' for(c in 1:cascades){
#'   simulated_cascade <- simulate_neha_discrete(x=data_for_sim,node="node",time="time",beta=beta,gamma=gamma,a=-6)
#'  simulated_cascade <- data.frame(simulated_cascade,cascade=c,stringsAsFactors=F)
#'   simulated_data <- rbind(simulated_data,simulated_cascade)
#' }
#'
#' # infer edges
#' neha_results <- neha(simulated_data,node="node",time="time",event="event",cascade="cascade",covariates="covariate",ncore=3)
#'
#' # estimate NEHA logistic regression
#' neha_estimate <- glm(neha_results$combined_formula,data=neha_results$data_for_neha,family=binomial)
#' summary(neha_estimate)
#' }
#'
#'
#' @export
neha <- function(eha_data,node,time,event,cascade,covariates,
                 ncore=2, negative=F){

  data_with_aest <- neha:::neha_geta(eha_data,node=node,time=time,event=event,cascade=cascade,covariates=covariates,ncore=ncore)

  a.estimate <- data_with_aest[[1]]

  a_est <- a.estimate

  data_with_aest <- data_with_aest[[2]]

  edge_vars <- names(data_with_aest)[which(!is.element(names(data_with_aest),names(eha_data)))]

  within_r_corr <- NULL
  for(i in 1:length(edge_vars)){
    ri <- strsplit(edge_vars[i],"_")[[1]][2]
    xi <- data_with_aest[which(data_with_aest[,node]==ri), edge_vars[i]]
    yi <- data_with_aest[which(data_with_aest[,node]==ri),event]
    #within_r_corr <- c(within_r_corr, wilcox.test(as.numeric(xi[which(yi==0)]),as.numeric(xi[which(yi==1)]),alternative="greater")$p.value)
    cor_ye <- suppressWarnings(cor(yi,data_with_aest[which(data_with_aest[,node]==ri), edge_vars]))
    cor_ye <- cor_ye[which(!is.na(cor_ye))]
    within_r_corr <- c(within_r_corr,(suppressWarnings(cor(yi,xi))-mean(cor_ye))/sd(cor_ye))

  }

  names(within_r_corr) <- edge_vars

  unodes <- unique(data_with_aest[,node])

  edges_subset <-NULL

  converged <- F

  iteration <- 1

  while(!converged){

    old_edges <- edges_subset

    effect_names <- c(covariates,edges_subset)

    if(!is.null(edges_subset)){
      new_a <- neha:::update_a(data_with_aest=data_with_aest,covariates=covariates,edges_subset=edges_subset,event=event,old_a=a_est,edge_vars=edge_vars,ncore=ncore)
      a_est <- new_a[[1]]
      data_with_aest <- new_a[[2]]
    }

    off <- rep(0,nrow(data_with_aest))

    if(length(covariates) > 0){

    full_estimate <- glm(data_with_aest[,event] ~ as.matrix(data_with_aest[,effect_names]),family=binomial)

    off <- as.matrix(data_with_aest[,covariates])%*%(coef(full_estimate)[2:(length(covariates)+1)])

    }

    data_with_aest$off <- off

    cl <- makeCluster(ncore)
    registerDoParallel(cl)

    edges_subset <- foreach(reciever = unodes,.packages=c("glmulti")) %dopar% {

      edger <- do.call('rbind',strsplit(names(within_r_corr),"_"))[,2]

      corr_r <- within_r_corr[which(edger==reciever)]

      #corr_r <- corr_r[which(corr_r>1.96)]
      if(length(corr_r) > 10) corr_r <- corr_r[order(-corr_r)[1:10]]
      screenedr <- names(corr_r)
      if(length(screenedr)>0){
        #Xyr <- data_with_aest[which(data_with_aest[,node]==reciever),c(screenedr,event)]
        Xyr <- data_with_aest[,c(screenedr,event)]
        names(Xyr)[ncol(Xyr)] <- "y"
        offr <- off[which(data_with_aest[,node]==reciever)]

        sdxy <- apply(Xyr,2,sd)

        Xyr <- Xyr[,which(sdxy > 0)]

        increased <- T

        est0 <- glm(Xyr[,"y"]~1,family=binomial)

        BICfit <- BIC(est0)

        vars <- 1
        res.best.logistic <- NULL
        while(increased){

          res.best.logistic.v <-
            glmulti(y="y",xr=names(Xyr)[-ncol(Xyr)], data = Xyr,
                    level = 1,               # No interaction considered
                    method = "h",            # Exhaustive approach
                    crit = "bic",            # AIC as criteria
                    confsetsize = 1,         # Keep 5 best models
                    plotty = F, report = F,  # No plot or interim reports
                    fitfunction = "glm",     # glm function
                    maxsize = vars,
                    minsize = ifelse(vars==1,0,vars),
                    family = binomial,offset=off)
          BICvars <- BIC(attributes(res.best.logistic.v)$objects[[1]])
          increased <- BICvars < BICfit

          if(increased){
            res.best.logistic <- res.best.logistic.v
            BICfit <- BICvars
            vars <- vars + 1
          }

        }

        best_vars <- names(coef(res.best.logistic))

        best_edges <- best_vars[which(is.element(best_vars,screenedr))]

        best_edges

      }

    }

    stopCluster(cl)

    edges_subset <- unique(c(unlist(edges_subset)))

    converged <- all(is.element(edges_subset,old_edges)) & all(is.element(old_edges,edges_subset))

  }

  edges_inferred <- edges_subset

  if(negative){

    edges_subset_pos <- edges_subset

    edges_subset <-NULL

    converged <- F

    iteration <- 1

    while(!converged){

      old_edges <- edges_subset

      effect_names <- c(covariates,edges_subset,edges_subset_pos)

      off <- rep(0,nrow(data_with_aest))

      if(length(covariates) > 0){

        full_estimate <- glm(data_with_aest[,event] ~ as.matrix(data_with_aest[,effect_names]),family=binomial)

        off <- as.matrix(data_with_aest[,covariates])%*%(coef(full_estimate)[2:(length(covariates)+1)])

      }

      data_with_aest$off <- off

      cl <- makeCluster(ncore)
      registerDoParallel(cl)

      edges_subset <- foreach(reciever = unodes,.packages=c("glmulti")) %dopar% {

        edger <- do.call('rbind',strsplit(names(within_r_corr),"_"))[,2]

        corr_r <- within_r_corr[which(edger==reciever)]

        corr_r <- corr_r[which(corr_r < 0)]
        if(length(corr_r) > 10) corr_r <- corr_r[order(corr_r)[1:10]]
        screenedr <- names(corr_r)
        if(length(screenedr)>0){
          #Xyr <- data_with_aest[which(data_with_aest[,node]==reciever),c(screenedr,event)]
          Xyr <- data_with_aest[,c(screenedr,event)]
          names(Xyr)[ncol(Xyr)] <- "y"
          offr <- off[which(data_with_aest[,node]==reciever)]

          sdxy <- apply(Xyr,2,sd)

          Xyr <- Xyr[,which(sdxy > 0)]

          increased <- T

          est0 <- glm(Xyr[,"y"]~1,family=binomial)

          BICfit <- BIC(est0)

          vars <- 1
          res.best.logistic <- NULL
          while(increased){

            res.best.logistic.v <-
              glmulti(y="y",xr=names(Xyr)[-ncol(Xyr)], data = Xyr,
                      level = 1,               # No interaction considered
                      method = "h",            # Exhaustive approach
                      crit = "bic",            # AIC as criteria
                      confsetsize = 1,         # Keep 5 best models
                      plotty = F, report = F,  # No plot or interim reports
                      fitfunction = "glm",     # glm function
                      maxsize = vars,
                      minsize = ifelse(vars==1,0,vars),
                      family = binomial,offset=off)
            BICvars <- BIC(attributes(res.best.logistic.v)$objects[[1]])
            increased <- BICvars < BICfit

            if(increased){
              res.best.logistic <- res.best.logistic.v
              BICfit <- BICvars
              vars <- vars + 1
            }

          }


          best_vars <- names(coef(res.best.logistic))

          best_edges <- best_vars[which(is.element(best_vars,screenedr))]

          best_edges

        }

      }

      stopCluster(cl)

      edges_subset <- unique(c(unlist(edges_subset)))

      converged <- all(is.element(edges_subset,old_edges)) & all(is.element(old_edges,edges_subset))

    }

    edges_inferred <- c(edges_subset_pos,edges_subset)

  }

  edge_cols <- data_with_aest[,edges_inferred]
  edges_combined <- apply(edge_cols,1,sum)
  edge_dat <- data.frame(edge_cols,edges_combined)
  data_for_neha <- data.frame(data_with_aest[,c(node,time,event,cascade,covariates)],edge_dat)
  combined_formula <- as.formula(paste(event,"~",paste(c(covariates,"edges_combined"),collapse="+"),sep=""))
  separate_formula <- as.formula(paste(event,"~",paste(c(covariates,edges_inferred),collapse="+"),sep=""))


  # results
  list(a_est = a_est,edges=edges_inferred,data_for_neha = data_for_neha,combined_formula=combined_formula,separate_formula=separate_formula)

}
