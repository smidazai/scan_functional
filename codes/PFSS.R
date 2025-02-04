# compute the PFSS method
compute_p <- function(c1, c2, my_mat) {
  
  # number of combinaison
  nb_combi <- length(c1)
  stat_p <- numeric(nb_combi) 
  npoints <- nrow(my_mat)
  
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
    
    ##### Initi
    
    X_bar = rowMeans(myX)
    Y_bar = rowMeans(myY)
    
    ##### Parametric
    moy_glob <- ((nx * X_bar) + (ny * Y_bar))/(nx + ny)
    
    diff1 <- X_bar - moy_glob
    diff2 <- Y_bar - moy_glob
    nor1 <- norm(diff1)^2
    nor2 <- norm(diff2)^2
    numerateur <- (nx * nor1) + (ny * nor2)
    
    sum1 <- 0
    for (j in 1:nx){
      temp <- myX[, j] - X_bar
      sum1 <- sum1 + crossprod(temp)
    }
    
    sum2 <- 0
    for (j in 1:ny){
      temp <- myY[, j] - Y_bar
      sum2 <- sum2 + crossprod(temp)
    }
    
    denominateur <- (1/((nx + ny)-2))*(sum1+sum2)
    
    stat_p[i] <- numerateur / denominateur
    
  }
  
  my_max_p <- max(stat_p)
  vecclus_p <-  c1[[which.max(stat_p)]]
  
  # secondary cluster 
  order_stat <- order(stat_p, decreasing = T)
  c1_order <- c1[order_stat]
  c2_order <- c2[order_stat]
  
  list(
    stat=my_max_p, vec=vecclus_p
  )
}