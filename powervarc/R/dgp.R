
rm(list = ls())

library('plm')
library('lmtest')
library('multiwayvcov')
library('miceadds')
library('parallel')
library('clubSandwich')
source('C:/Users/asanc/OneDrive/Documents/GitHub/Econometrics-IV/bindToEnv.R')

#' Randomly assign units to treatment and control groups
#'
#' @param data The dataset to be used
#' @param clustervar The cluster identifier
#' @param p_vec A data frame that contains p_I, p_K
#' @param theta A data frame that conntaine a0,a1
#'
#' @return A dataset with the treatment variables
#' @export
#'
#' @examples
sample_randomize <- function(data,clustervar,p_vec,theta) {

  # Extract parameters of interest
  p_K   <- p_vec$p_K
  p_I   <- p_vec$p_I
  a0    <- theta$a0
  a1    <- theta$a1

  # Identify unique clusters
  cluster_id <- data[,c(clustervar)]
  unique_cluster_id <- unique(cluster_id)

  # Identify number of observations
  num_obs      <- dim(data)[1]
  num_clusters <- length(unique_cluster_id)

  # Randomize
  rand_obs      <- rbinom(num_obs,1,p_I)
  rand_clusters <- rbinom(num_clusters,1,p_K)

  treatment_cluster <- as.data.frame(cbind(unique_cluster_id,rand_clusters))
  data <- merge(data,treatment_cluster,by.x=clustervar,by.y="unique_cluster_id")

  data$rand_obs  <- rand_obs
  data$treatment <- (data$rand_obs)*(data$rand_clusters)

  rand_spillover     <- rbinom(num_obs,1,a0)
  rand_noncompliance <- rbinom(num_obs,1,1-a1)

  data$treatment_tilde <- (data$rand_obs)*(data$rand_clusters)*rand_noncompliance +
                          (1-data$rand_obs)*(data$rand_clusters)*rand_spillover

  return(data)

}


#' Generate data for the outcome variable
#'
#' @param data  A dataframe that contain the "treatment_tilde" variable.
#' @param treatment_tilde A string with the name of the treatment_tilde variable
#' @param clustervar A string with the name of the clustervar variable
#' @param theta A dataframe that contains the variables beta,sigma_u,sigma_e
#'
#' @return
#' @export
#'
#' @examples
sample_outcome <- function(data,treatment_tilde,clustervar,theta) {

  cluster_id <- data[,c(clustervar)]
  unique_cluster_id <- unique(cluster_id)

  num_obs      <- dim(data)[1]
  num_clusters <- length(unique_cluster_id)

  beta    <- theta$beta
  sigma_u <- theta$sigma_u
  sigma_e <- theta$sigma_e

  rand_e <- rnorm(num_obs,0,sigma_e)
  rand_u <- rnorm(num_clusters,0,sigma_u)

  rand_cluster <- as.data.frame(cbind(unique_cluster_id,rand_u))

  data <- merge(data,rand_cluster,by.x=clustervar,by.y="unique_cluster_id")

  data$outcome <- data$rand_u + rand_e + data[,c(treatment_tilde)]*beta 

  return(data)

}


sample_cluster <- function(clustervar,theta) {

  cluster_num  <- theta$cluster_num
  cluster_size <- theta$cluster_size
  
  unique_cluster_id <- (1:cluster_num)
  # data <- as.data.frame(kronecker(unique_cluster_id,seq(1,1,length.out=cluster_size)))
  #  data <- expand.grid(cluster_id = 1:cluster_num,clustervar = 1:cluster_num )
  data <- data.frame(clustervar = rep(1:cluster_num,cluster_size))
  # names(data) <- c(clustervar)

  return(data)

}

simulate_test <- function(p_vec,theta) {

  a_vec         <- c(theta$a0,theta$a1)
  cluster_num   <- theta$cluster_num
  cluster_size  <- theta$cluster_size
  param         <- c(theta$beta,theta$sigma_u,theta$sigma_e)
  
  # Create simulate
  sample <- sample_cluster('clustervar',theta)
  sample <- sample_randomize(sample,'clustervar',p_vec,theta)
  sample <- sample_outcome(sample,'treatment_tilde','clustervar',theta)
  # sample <- invisible(plm::pdata.frame(sample,index = 'clustervar'))

  # Carry out hypothesis test
  # results_reg <- miceadds::lm.cluster(sample,outcome ~ treatment_tilde,'clustervar')
  # tstat       <- coeftest(results_reg)[2,3]

  # sample$treatment_tilde <- rep(1,cluster_num*cluster_size)
  
  # results_reg <- plm(outcome ~ treatment_tilde, data = sample, model = "pooling")
  results_reg <- lm(outcome ~ treatment_tilde ,sample)
  # cluster_vcov <- plm::vcovHC(results_reg, type = "HC0", cluster = "group", adjust = T)
  
  if(sum(!is.na(results_reg$coefficients)) > 1) {

    # pvalue <- clubSandwich::coef_test(results_reg, vcov = "CR1",
    #                                  cluster = sample$clustervar, test = "naive-t")[2,3]

    pvalue <- clubSandwich::coef_test(results_reg, vcov = "CR2",
                                      cluster = sample$clustervar, test = "Satterthwaite")[2,4]
    
    summary(results_reg)  
    summary(results_reg,cluster = c("clustervar"))
    
    # results_reg <- lm(outcome ~ treatment_tilde,sample)
    # tstat       <- coeftest(results_reg)[2,3]
    
    reject      <- 1*(pvalue < 0.05)
    
    varCov <- as.matrix(sample$outcome) %*% as.matrix(t(sample$outcome))
    
  } else {

    pvalue <- NA
    reject <- NA

  }
  
  return(reject)

}


test_fn <- function(i,p_vec,theta) {
  
  bindToEnv(objNames=c('sample_randomize','sample_outcome','sample_cluster',
                       'simulate_test','coeftest','plm'))
  
  function(x) {
    simulate_test(p_vec[i,],theta[i,])
  }
  
}




# Input parameters
# p_vec <- data.frame(p_I = c(0.6),p_K = 0.2)
theta <- data.frame(a0 = 1,a1 = 0, cluster_num = 100,cluster_size = 50,
                    beta=0,sigma_u = 1,sigma_e = 1)

p_vec <- expand.grid(p_I = seq(0.001,1,length.out = 8), p_K=seq(0.001,1,length.out = 8))
theta <- theta[rep(row.names(theta), nrow(p_vec)),]

# Simulation with the first set of parameters
# simulate_test(p_vec[1,],theta[1,])

parallelCluster <- parallel::makeCluster(parallel::detectCores())
results <- data.frame(p_vec,power = NA)

for(i in 1:nrow(p_vec)) {
  
  print(i)
  
  ptm <- proc.time()
  
  simulate_power <- test_fn(i,p_vec,theta)
  simulate_power()
  hola <- parallel::parLapply(parallelCluster,rep(0,800),simulate_power)
  
  results[i,c('power')] <- mean(unlist(hola),na.rm=TRUE)
  
  proc.time() - ptm
  
}

sum(!is.na(unlist(hola)))

on.exit(stopCluster(cl))

