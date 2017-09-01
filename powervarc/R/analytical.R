
setwd(paste('D:/Dropbox/UPenn/Research Projects/',
            'Spillover effects power calculation/figures/simulations',
            sep = ""))

library(ggplot2)
library(plot3D)

power_uc <- function(p_vec,theta) {
  
  p_K  <- p_vec$p_K
  p_I  <- p_vec$p_I
  a0   <- theta$a0
  a1   <- theta$a1
  K    <- theta$cluster_num
  M_K  <- theta$cluster_size
  
  beta <- theta$beta
  sigma_u <- theta$sigma_u
  sigma_e <- theta$sigma_e

  E_theta <- (1-a1-((1-p_I)*p_K/(1-p_I*p_K))*a0)*beta
  V_theta <- (1/(K*M_K))*((1/((p_I*p_K)*(1-p_I*p_K)))*(sigma_e^2) +
                          (M_K*(sigma_u^2))*((1-p_K)/((1-p_I*p_K)^2)) +
                          (beta^2)*(((a0)*(1-a0)/((1-p_I*p_K)^2)) + 
                                    ((a1)*(1-a1)/((p_I*p_K)^2))))
  
  V_theta_miss <- (1/(K*M_K))*((1/((p_I*p_K)*(1-p_I*p_K)))*(sigma_e^2) +
                               (M_K*(sigma_u^2))*((1-p_K)/((1-p_I*p_K)^2)))
  
  factor_miss <- 1 # sqrt(V_theta_miss/V_theta) # 1
  
  power <- 1-(pnorm(1.96*factor_miss-E_theta/sqrt(V_theta))-
              pnorm(-1.96*factor_miss-E_theta/sqrt(V_theta)))
                          
  return(power)

}


example_plot <- function(){
  
  theta <- data.frame(a0 = 0.33,a1 = 0.17, cluster_num = 22,cluster_size = 90,
                      beta=0.86,sigma_u = 0.84,sigma_e = 4.95)
  
  p_vec <- expand.grid(p_I = seq(0.001,0.999,length.out = 100),
                       p_K=seq(0.001,0.999,length.out = 100))
  
  power_vec <- power_uc(p_vec, theta)
  data <- data.frame(p_vec,power_vec <- power_vec)
  
  ggplot(data) + 
    aes(x = p_I, y = p_K, z = power_vec, fill = power_vec) + 
    geom_tile() + 
    coord_equal() +
    geom_contour(color = "white", alpha = 0.5) + 
    scale_fill_distiller(palette="Spectral", na.value="white") + 
    theme_bw() +
    ggtitle(paste('a0 = ',theta$a0,', a1 = ',theta$a0,
                  ', K = ', theta$cluster_num,', M_K',theta$cluster_size),
            paste('b = ',theta$beta,', s_u = ',theta$sigma_u,
                  ', s_e = ', theta$sigma_e))
}


# Generate data

example_random <- function(S) {

  results <- matrix(NA,S,10)
  
  set.seed(12345)
  theta_df <- data.frame(a0 = 0.5*runif(S),a1 = 0.5*runif(S),
                         cluster_num = 100*runif(S),cluster_size = 100*runif(S),
                         beta=runif(S),sigma_u = 5*runif(S),sigma_e = 5*runif(S))
  
  p_vec <- expand.grid(p_I = seq(0.001,0.999,length.out = 100),
                       p_K=seq(0.001,0.999,length.out = 100))
  
  for(i in 1:S) {
    
    print(i)
  }
  
  for(i in 1:S) {
    
    print(i)
    theta <- theta_df[i,]
    
    power_vec <- power_uc(p_vec, theta)
    data <- data.frame(p_vec,theta[rep(rownames(theta),nrow(p_vec)),],
                       power_vec = power_vec)                      
    
    # plot_power <- ggplot(data) + 
    #   aes(x = p_I, y = p_K, z = power_vec, fill = power_vec) + 
    #   geom_tile() + 
    #   coord_equal() +
    #   geom_contour(color = "white", alpha = 0.5) + 
    #   scale_fill_distiller(palette="Spectral", na.value="white",limits = c(0,1)) + 
    #   theme_bw() +
    #   ggtitle(paste('a0 = ',theta$a0,', a1 = ',theta$a1,
    #                 ', K = ', theta$cluster_num,', M_K',theta$cluster_size),
    #           paste('b = ',theta$beta,', s_u = ',theta$sigma_u,
    #                 ', s_e = ', theta$sigma_e))
    # 
    # 
    # print(plot_power)
    # dev.copy(png,paste('rplot_',i,'.png',sep=""))
    # dev.off()
    
    results_tentative <- as.matrix(data[power_vec == max(power_vec),])
    
    if(dim(results_tentative)[1] == 1) {
      
      results[i,] <- results_tentative
      
    }
    
    colnames(results) <- c('p_I','p_K','a0','a1','K','M_K','beta','s_u','s_e','pow')

  }
  
  return(round(results,2))
  
}

example_random(1000)


