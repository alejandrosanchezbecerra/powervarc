
library(ggplot2)

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
                          (beta^2)*(((a0)*(1-a0)/((1-p_I*p_K)^2)) + 
                                    ((a1)*(1-a1)/((p_I*p_K)^2))) +
                          (M_K*(sigma_u^2))*((1-p_K)/((1-p_I*p_K)^2)))
  
  power <- 1-(pnorm(1.96-E_theta/sqrt(V_theta))-pnorm(-1.96-E_theta/sqrt(V_theta)))
                          
  return(power)

}

theta <- data.frame(a0 = 1,a1 = 0, cluster_num = 100,cluster_size = 50,
                    beta=0,sigma_u = 1,sigma_e = 1)

p_vec <- expand.grid(p_I = seq(0,1,length.out = 10), p_K=seq(0,1,length.out = 10))

power_uc(p_vec, theta)

