



N <- 100

theta_vec <-  seq(0,2*pi,length.out = N)
t_vec <-  seq(0,1,length.out = N)

results <- matrix(NA,N,N)

target_fn <- function(a0,a1) {

  a1*(1-a1)*((1-a0)^2)-((1-a0-a1)^2)*(1-a1)*(a1)-((1-a0-a1)^2)*(1-a0

}

target_fn(0,0)


i <- 0

for( theta in theta_vec ) {

  i <- i + 1
  j <- 0

  for( t in t_vec) {

    j <- j + 1

    results[i,j] = abs( target_fn(t*cos(theta),t*sin(theta))) / t

  }


}

sum(is.na(results[,-1]))

heatmap(results[,-1],Colv = NULL)
