# compute DFFSS method # 
compute_dffss <- function(c1, c2, my_mat) {
  
  # number of combinaison
  nb_combi <- length(c1)
  dffss <- numeric(nb_combi)
  npoints <- nrow(my_mat)
  
  for (i in (1:nb_combi)) {
    vecindin <- c1[[i]]
    vecindout <- c2[[i]]
    
    nx <- length(vecindin)
    ny <- length(vecindout)
    
    if(nx == 1) {
      myX <- matrix(my_mat[, vecindin], nrow(my_mat), 1)
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
    Xn_bar_matrix = matrix(X_bar, npoints, nx)
    
    Y_bar = rowMeans(myY)
    Xm_bar_matrix = matrix(Y_bar, npoints, ny)
    
    sd_x <- rowSums((myX - Xn_bar_matrix)^2)
    sd_y <- rowSums((myY - Xm_bar_matrix)^2)
    sd_den <- (sd_x + sd_y) / (nx + ny - 2) * (1/nx + 1/ny)
    dffss[i] <- max(abs(X_bar - Y_bar) / sqrt(sd_den), na.rm = T) 
  }
  
  my_max_dffss <- max(dffss)
  vecclus_dffss <-  c1[[which.max(dffss)]]
  
  list(
    stat=my_max_dffss, 
    vec=vecclus_dffss
  )
}