library(testthat)
library(neha)

test_that("data_neha_discrete works", {
  # basic data parameters
  nodes <- 3
  times <- 3
  # generate dataframes
  time <- rep(1:times,nodes)
  node <- paste("n",as.character(rep(1:nodes,times)),sep="")
  intercept <- rep(1,length(time))
  covariate <- runif(length(time))-2
  data_for_sim1 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim2 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim3 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim1$cascade = 1
  data_for_sim2$cascade = 2
  data_for_sim3$cascade = 3

  data_for_sim1$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim1 <- data_for_sim1[-c(2,3,6),]
  data_for_sim2$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim2 <- data_for_sim1[-c(2,3,6),]
  data_for_sim3$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim3 <- data_for_sim1[-c(2,3,6),]

  data_for_sim <- rbind(data_for_sim1,data_for_sim2,data_for_sim3)

  test_dat <- data_neha_discrete(data_for_sim,
                                 node='node',
                                 time='time',
                                 event='event',
                                 cascade='cascade')

  expect_equal(all(c(nrow(test_dat)==12,ncol(test_dat)==9)), TRUE)
})


test_that("simulate_neha_discrete works", {
  # Simulate data for NEHA
  # basic data parameters
  set.seed(9202011)
  nodes <- 3
  times <- 3
  nties <- 2
  # generate dataframe
  time <- sort(rep(1:times,nodes))
  node <- paste("n",as.character(rep(1:nodes,times)),sep="")
  intercept <- rep(1,length(time))
  covariate <- runif(length(time))-2
  data_for_sim <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)

  # regression parameters
  beta <- cbind(c(-1.5,.25))
  rownames(beta) <- c("intercept","covariate")

  # generate network effects
  possible_ties <- rbind(t(combn(1:nodes,2)),t(combn(1:nodes,2))[,c(2,1)])
  possible_ties <- paste(paste("n",possible_ties[,1],sep=""),paste("n",possible_ties[,2],sep=""),sep="_")
  ties <- sample(possible_ties,nties)
  gamma <- cbind(rep(1.5,length(ties)))
  rownames(gamma) <- ties

  # initiate simulated data object
  simulated_cascade <- simulate_neha_discrete(x=data_for_sim,node="node",time="time",beta=beta,gamma=gamma,a=-6)

  expect_equal(all(c(nrow(simulated_cascade)==6,
                     ncol(simulated_cascade)==5,
                     sum(simulated_cascade$event) == 2)), TRUE)
})


test_that("neha_geta works", {
  # basic data parameters
  nodes <- 3
  times <- 3
  # generate dataframes
  time <- rep(1:times,nodes)
  node <- paste("n",as.character(rep(1:nodes,times)),sep="")
  intercept <- rep(1,length(time))
  covariate <- runif(length(time))-2
  data_for_sim1 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim2 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim3 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim1$cascade = 1
  data_for_sim2$cascade = 2
  data_for_sim3$cascade = 3

  data_for_sim1$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim1 <- data_for_sim1[-c(2,3,6),]
  data_for_sim2$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim2 <- data_for_sim1[-c(2,3,6),]
  data_for_sim3$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim3 <- data_for_sim1[-c(2,3,6),]

  data_for_sim <- rbind(data_for_sim1,data_for_sim2,data_for_sim3)

  test_dat <- neha_geta(data_for_sim,
                        node='node',
                        time='time',
                        event='event',
                        cascade='cascade')

  expect_equal(log(log(0.9)/(1-2)) == test_dat$a_est, TRUE)
})



test_that("update_a works", {
  # basic data parameters
  nodes <- 3
  times <- 3
  # generate dataframes
  time <- rep(1:times,nodes)
  node <- paste("n",as.character(rep(1:nodes,times)),sep="")
  intercept <- rep(1,length(time))
  covariate <- rep(1,length(time))
  data_for_sim1 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim2 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim3 <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)
  data_for_sim1$cascade = 1
  data_for_sim2$cascade = 2
  data_for_sim3$cascade = 3

  data_for_sim1$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim1 <- data_for_sim1[-c(2,3,6),]
  data_for_sim2$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim2 <- data_for_sim1[-c(2,3,6),]
  data_for_sim3$event <- c(1,0,0,0,1,rep(0,3),1)
  data_for_sim3 <- data_for_sim1[-c(2,3,6),]

  data_for_sim <- rbind(data_for_sim1,data_for_sim2,data_for_sim3)

  test_dat <- data_neha_discrete(data_for_sim,
                                 node='node',
                                 time='time',
                                 event='event',
                                 cascade='cascade')

  ua <- neha:::update_a(test_dat,covariates="covariate",
                 edges_subset=c("n1_n2","n1_n3","n2_n3"),
                 edge_vars=c("n1_n2","n1_n3","n2_n3"),
                 old_a =2,
                 event="event",
                 ncore=2)

  expect_equal(ua[[1]],2.2)
})



test_that("neha works", {
  set.seed(9202011)
  cascades <- 10
  nodes <- 5
  times <- 10
  nties <- 4
  # generate dataframe
  time <- sort(rep(1:times,nodes))
  node <- paste("n",as.character(rep(1:nodes,times)),sep="")
  intercept <- rep(1,length(time))
  covariate <- runif(length(time))-2
  data_for_sim <- data.frame(time,node,intercept,covariate,stringsAsFactors=F)

  # regression parameters
  beta <- cbind(c(-2.5,.25))
  rownames(beta) <- c("intercept","covariate")

  # generate network effects
  possible_ties <- rbind(t(combn(1:nodes,2)),t(combn(1:nodes,2))[,c(2,1)])
  possible_ties <- paste(paste("n",possible_ties[,1],sep=""),paste("n",possible_ties[,2],sep=""),sep="_")
  ties <- sample(possible_ties,nties)
  gamma <- cbind(rep(1.5,length(ties)))
  rownames(gamma) <- ties

  # initiate simulated data object
  simulated_data <- NULL

  # generate the data one cascade at a time
  for(c in 1:cascades){
    simulated_cascade <- simulate_neha_discrete(x=data_for_sim,node="node",time="time",beta=beta,gamma=gamma,a=-6)
    simulated_cascade <- data.frame(simulated_cascade,cascade=c,stringsAsFactors=F)
    simulated_data <- rbind(simulated_data,simulated_cascade)
  }
  # estimate NEHA
  neha_results <- neha(simulated_data,node="node",time="time",event="event",cascade="cascade",covariates="covariate",ncore=3)

  expect_equal(all(length(neha_results[[2]]) == 2,
                   sum(is.element(c("n5_n3","n4_n5") ,neha_results[[2]]))),TRUE)
})










