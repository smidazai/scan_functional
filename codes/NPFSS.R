###########################################
norm <- function(x) sqrt(sum(x^2))

# compute NPFSS method
compute_np <- function(c1, c2, my_mat) {
  # number of combinaison
  nb_combi <- length(c1)
  npoints <- nrow(my_mat)
  
  my_dist <- as(dist(t(my_mat)), "matrix")
  
  stat_np <- numeric(nb_combi)
  
  for (i in (1:nb_combi)) {
    vecindin <- c1[[i]]
    vecindout <- c2[[i]]
    
    nx <- length(vecindin)
    ny <- length(vecindout)
    
    if(nx == 1) {
      myX <- matrix(my_mat[,vecindin], nrow(my_mat), 1)
    } else {
      myX <- my_mat[, vecindin]
    }
    
    if(ny == 1) {
      myY <- matrix(my_mat[,vecindout], nrow(my_mat), 1)
    } else {
      myY <- my_mat[, vecindout]
    }
    
    ##### Non Parametric
    vec_int <- numeric(npoints)
    
    for(j in 1:nx)
      vec_int <- vec_int + my_mat[, vecindin[j]] * sum(1/my_dist[vecindin[j], vecindout])
    
    for(j in 1:ny)
      vec_int <- vec_int - my_mat[, vecindout[j]] * sum(1/my_dist[vecindout[j], vecindin])
    
    vec_int <- vec_int / (nx * ny)
    temp_stat_np <- norm(vec_int)
    
    stat_np[i] <-  sqrt((nx * ny) / (nx + ny)) * temp_stat_np
  }
  
  my_max_np <- max(stat_np)
  vecclus_np <-  c1[[which.max(stat_np)]]
  
  list(
    stat=my_max_np, vec=vecclus_np
  )
}
