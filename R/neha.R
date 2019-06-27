
#' A function to create data for NEHA with discrete time EHA data
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param node A character string name of the variable that gives the node id.
#' @param time A character string name of the variable that gives the time, in integers.
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param cascade A character string name of the variable that gives the cascade id
#' @param threshold An integer such that an edge variable for node pair i/j will not be constructed if j did not experience more than 'threshold' events after i experienced them.
#' @return A data frame in which, in addition to all of the variables in 'eha_data', there is one column for each directed dyad that meets the 'threshold' condition, named 'i_j' in which the value indicates the number of time points before 'time' that 'i' experienced the event in the respective 'cascade'. If 'i' did not experience the cascade before 'time', the value is 0.
#' @export
data_neha_discrete <- function(eha_data,node,time,event,cascade,covariates){

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

#' A function to estimate NEHA parameters. This is used internally by bolasso.neha.
#' @import glmnet
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param A character string name of the variable that gives the node id
#' @param time A character string name of the variable that gives the time, in integers
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param A character string name of the variable that gives the cascade id
#' @param threshold An integer such that an edge variable for node pair i/j will not be constructed if j did not experience more than 'threshold' events after i experienced them.
#' @param covariates character vector of covariate names to include in the neha, excluding the intercept.
#' @param a A non-negative numeric value that models the exponential decay of sender influence. The effect of a previous event experienced by a diffusion tie sender on the log odds of an event at time t is gamma*exp(-exp(a)*(t-source_time)), where source_time is the time the sender experienced the event.
# note1: this function does not currently permit dyadic covariates
neha <- function(data_for_neha,node,time,event,cascade,covariates,threshold=0,a=-8){
  diffusion_effects_variables <- data_for_neha[,(length(covariates)+5):ncol(data_for_neha)]
  diffusion_effects_variables <- (diffusion_effects_variables>0)*exp(-exp(a)*(diffusion_effects_variables))
  if(length(covariates) == 0){
    x_for_glmnet <- as.matrix(diffusion_effects_variables)
    colnames(x_for_glmnet) <- colnames(diffusion_effects_variables)
  }
  if(length(covariates) > 0){
    covariate_variables <- cbind(data_for_neha[,covariates])
    x_for_glmnet <- cbind(as.matrix(covariate_variables),as.matrix(diffusion_effects_variables))
    colnames(x_for_glmnet) <- c(covariates,colnames(diffusion_effects_variables))
  }

  y_for_glmnet <- data_for_neha[,event]

  penalty.factors <- c(rep(0,length(covariates)),rep(1,ncol(diffusion_effects_variables)))

  lower.vals <- c(rep(-Inf,length(covariates)),rep(0,ncol(diffusion_effects_variables)))

  neha_estimate <- glmnet::cv.glmnet(x_for_glmnet,y_for_glmnet,family="binomial",penalty.factor=penalty.factors,lower.limits=lower.vals,nfolds=10,alpha=1,type.logistic="modified.Newton")

  neha_estimate

}


#' A function for NEHA estimation with bolasso (Bach 2008) selection.
#'
#' @import doParallel
#' @import parallel
#' @import speedglm
#' @param eha_data A dataframe that includes one observation for each node at risk of experiencing the event during each at-risk time point in each cascade. Note, it is assumed that each node can experience an event in each cascade once, at most.
#' @param node A character string name of the variable that gives the node id
#' @param time A character string name of the variable that gives the time, in integers
#' @param event A character string name of the variable that gives the binary 0/1 indicator of event occurrence.
#' @param cascade A character string name of the variable that gives the cascade id
#' @param covariates character vector of covariate names to include in the neha, excluding the intercept.
#' @param a A non-negative numeric value that models the exponential decay of sender influence. The effect of a previous event experienced by a diffusion tie sender on the log odds of an event at time t is gamma*exp(-exp(a)*(t-source_time)), where source_time is the time the sender experienced the event.
#' @param estimate.a A logical value indicating whether or not to estimate the value of 'a'.
#' @param n_jobs Integer, number of jobs to run in parallel. -1 means using all processors.
#'
#' @references Bach, Francis R. "Bolasso: model consistent lasso estimation through the bootstrap." In Proceedings of the 25th international conference on Machine learning, pp. 33-40. ACM, 2008.
#' @return A list with the following objects.
#' \itemize{
#'   \item neha.estimate - glm object giving the final neha estimates.
#'   \item a.estimate - estimate of a.
#'   \item formula.neha - formula needed to run final neha specification.
#'   \item data.for.neha - data frame that includes the edge variables required to run final neha specification.
#' }
#' @examples
#' \dontrun{
#' # Simulation study of the precision and recall of NEHA
#' # basic data parameters
#' cascades <- 100
#' nodes <- 10
#' times <- 20
#' nties <- 25
#'
#' # generate dataframe
#' time <- sort(rep(1:times,nodes))
#' node <- as.character(rep(1:nodes,times))
#' intercept <- rep(1,length(time))
#' covariate <- runif(length(time))-2
#' data_for_sim <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
#'
#' # regression parameters
#' beta <- cbind(c(-2,.25))
#' rownames(beta) <- c("intercept","covariate")

#' # generate network effects
#' possible_ties <- rbind(t(combn(1:nodes,2)),t(combn(1:nodes,2))[,c(2,1)])
#' possible_ties <- paste(possible_ties[,1],possible_ties[,2],sep="_")
#' ties <- sample(possible_ties,nties)
#' gamma <- cbind(exp(rnorm(nties)/2))
#' rownames(gamma) <- ties

#' # initiate simulated data object
#' simulated_data <- NULL

#' # generate the data one cascade at a time
#' for(c in 1:cascades){
#'   simulated_cascade <- simulate_neha_discrete(x=data_for_sim,node="node",time="time",beta=beta,gamma=gamma,a=-8)
#'   simulated_cascade <- data.frame(simulated_cascade,cascade=c,stringsAsFactors=F)
#'   simulated_data <- rbind(simulated_data,simulated_cascade)
#' }
#'
#' bolasso.results <- bolasso.neha(simulated_data,node="node",time="time",event="event",cascade="cascade",covariates="covariate",a=-8)
#'
#' bolasso.coefs <- names(coef(bolasso.results[[1]]))[-(1:length(beta))]
#' recall <- mean(is.element(paste('e',ties,sep=""),bolasso.coefs))
#' precision <- mean(is.element(bolasso.coefs,paste('e',ties,sep="")))
#' recall
#' precision
#' summary(bolasso.results)
#' }
#' @export
bolasso.neha <- function(eha_data,node,time,event,cascade,covariates=NULL,a=-8,estimate.a = T, n_jobs=1,data_only=F){
		
  # find out if node is numeric
  n1c <- substr(as.character(eha_data[,node]),1,1)
  
  if(!all(is.element(tolower(n1c),letters))){
  	print("appending n_ to beginning of all node id's since at least one seems numeric")
  	eha_data[,node] <- paste("n_",eha_data[,node],sep="")
  }


  eha_data <- eha_data[,c(node,time,event,cascade,covariates)]

  if(n_jobs == -1) n_jobs <- detectCores()

  data_for_neha <- data_neha_discrete(eha_data,node=node,time=time,event=event,cascade=cascade,covariates=covariates)

  print("data organization complete")

  if(data_only) return(data_for_neha)
  diffusion_effects_variables <- data_for_neha[,(ncol(eha_data)+1):ncol(data_for_neha)]

  if(length(covariates)==0){
    x_for_glmnet <- cbind(as.matrix(diffusion_effects_variables>0)*exp(-exp(a)*(diffusion_effects_variables)))
    colnames(x_for_glmnet) <- colnames(diffusion_effects_variables)
  }

  if(length(covariates) > 0){
    covariate_variables <- cbind(eha_data[,covariates])
    x_for_glmnet <- cbind(as.matrix(covariate_variables),as.matrix(diffusion_effects_variables>0)*exp(-exp(a)*(diffusion_effects_variables)))
    colnames(x_for_glmnet) <- c(covariates,colnames(diffusion_effects_variables))
  }

  y_for_glmnet <- eha_data[,event]

  a.likelihood <- function(a.par,diffusion_effects_variables,y){

    x <- cbind(as.matrix(diffusion_effects_variables>0)*exp(-exp(a.par)*(diffusion_effects_variables)))

    x <- as.matrix(x)

    x <- data.frame(x)

    form <- paste("y~",paste(names(x),collapse="+"))
    form <- as.formula(form)

    x$y <- y

    logLik(speedglm::speedglm(form,family=binomial(),data=x))
  }

  if(estimate.a){
   find.a <- optim(a,a.likelihood,method="BFGS",control=list(fnscale=-1),diffusion_effects_variables=diffusion_effects_variables,y=y_for_glmnet)
   a <- find.a$par
   print("a estimation complete")
  }

  nboot <- 15


  seeds <- as.integer(round(1+2147483640*runif(nboot)))


  boot.nonzero.est <- NULL
  for(boot.iter in 1:nboot){
    boot.ind <- NULL
    set.seed(seeds[boot.iter])
    unique.cascades <- unique(data_for_neha[,cascade])
    boot.samp <- sample(unique.cascades, length(unique.cascades), rep=T)
    #boot.data <- lapply(boot.samp, function(x){
    #  eha_data[eha_data[, cascade] == boot.samp[x], ]
    #})
    # boot.data <- data.frame(do.call('rbind', boot.data),stringsAsFactors=F)
    boot.data <- NULL
    for(c in boot.samp){
      boot.ind <- c(boot.ind,which(data_for_neha[, cascade] == c))
    }
    boot.data <- data_for_neha[boot.ind,]



    if(boot.iter >= 1){
      boot.neha_estimate <- neha(data_for_neha = boot.data,
                                 node = node,
                                 time = time,
                                 event = event,
                                 covariates = covariates,
                                 cascade = cascade,
                                 a = a)
      betas <- boot.neha_estimate$glmnet.fit$beta
      sel <- which(boot.neha_estimate$lambda == boot.neha_estimate$lambda.min)
      boot.neha.coef <- betas[, sel][betas[, sel] != 0]
      coef_names <- names(boot.neha.coef)
      not_sel <- c(1:length(covariates))
      boot.recovered_edges <- substr(coef_names[-not_sel], 1,
                                     nchar(coef_names[-not_sel]))
      boot.nonzero.est <- c(boot.nonzero.est,boot.recovered_edges)
      print(paste("bootstrap iteration",boot.iter))
    }


  }

  freq.boot.est <- table(boot.nonzero.est)
  bolasso.eff <- names(freq.boot.est)[which(freq.boot.est == nboot)]

  bolasso.eff <- bolasso.eff[is.element(bolasso.eff,colnames(diffusion_effects_variables))]

  bolasso.x <- cbind(cbind(as.matrix(diffusion_effects_variables>0)*exp(-exp(a)*(diffusion_effects_variables)))[,match(bolasso.eff,colnames(diffusion_effects_variables))])
  colnames(bolasso.x) <- bolasso.eff
  if(length(covariates) == 0){
    bolasso.x <- bolasso.x
    colnames(bolasso.x) <- switch((length(bolasso.eff)>0) + 1,NULL,paste("e_",bolasso.eff,sep=""))
  }
  if(length(covariates) > 0){
    bolasso.x <- cbind(covariate_variables,bolasso.x)
    colnames(bolasso.x) <- c(covariates,switch((length(bolasso.eff)>0) + 1,NULL,paste("e_",bolasso.eff,sep="")))
  }
  data.for.neha <- data.frame(y_for_glmnet,bolasso.x)
  names(data.for.neha) <- c("y_for_glmnet",colnames(bolasso.x))
  if(length(colnames(bolasso.x)) ==0) formula.neha <- as.formula("y_for_glmnet~1")
  if(length(colnames(bolasso.x)) > 0){
    formula.neha <- as.formula(paste("y_for_glmnet~",paste(colnames(bolasso.x),collapse="+"),collapse=""))
  }
  bolasso.est <- glm(formula.neha,family="binomial",data=data.for.neha,x=T,y=T)

  return(list(neha.estimate = bolasso.est,a.estimate=a,formula.neha=formula.neha,data.for.neha=data.frame(data.for.neha,eha_data[,c(node,time,event,cascade)])))
}




