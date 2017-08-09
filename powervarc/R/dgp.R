

sample_randomize <- function(data,treatment_name, clustervar,
                             p_cluster, p_individual) {

  # data <- sample_data
  # p_cluster    <- 0.5
  # p_individual <- 0.2
  # clustervar <- 'variable1'

  cluster_id <- data[,c(clustervar)]
  unique_cluster_id <- unique(cluster_id)

  num_obs      <- dim(sample_data)[1]
  num_clusters <- length(unique_cluster_id)

  rand_obs      <- rbinom(num_obs,1,p_individual)
  rand_clusters <- rbinom(num_clusters,1,p_cluster)

  treatment_cluster <- as.data.frame(cbind(unique_cluster_id,rand_clusters))
  data <- merge(data,treatment_cluster,by.x=clustervar,by.y="unique_cluster_id")

  data$rand_obs  <- rand_obs
  data$treatment <- (data$rand_obs)*(data$rand_clusters)

  return(data)

}


empty_matrix <- matrix(c(rbinom(1000,20,0.5),rnorm(1000)),1000,2)

sample_data <- as.data.frame(empty_matrix)
names(sample_data) <- c('variable1','variable2')

sample_data_transform <- sample_randomize(sample_data,'treatment','variable1',0.5,0.5)

sum(sample_data_transform$treatment)
