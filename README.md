## Model Overview 

An R package which implements network event history analysis (NEHA).

* Jeffrey J. Harden, Bruce A. Desmarais, Mark Brockway, Frederick J. Boehmke, Scott J. LaCombe, Fridolin Linder, and Hanna Wallach. A Diffusion Network Event History Estimator. The Journal of Politics, Forthcoming

## Installation

The package can be installed from GitHub via devtools.

    install.packages("devtools")
    
Now we can install from Github using the following line:

    devtools::install_github("desmarais-lab/neha")

## Example

```
library(neha)
# Simulate data for NEHA
# basic data parameters
cascades <- 50
nodes <- 20
times <- 30
nties <- 25
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

```
