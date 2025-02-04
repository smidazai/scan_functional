# compute HFSS method
compute_h <- function(c1, c2, my_mat, d = 5,
                      plot_eigen = F, ...) {
  # number of combinaison
  nb_combi <- length(c1)
  npoints <- nrow(my_mat)
  hotelling_k <- numeric(nb_combi)
  cpv <- numeric(npoints)
  
  for (i in (1:nb_combi)) {
    vecindin <- c1[[i]]
    vecindout <- c2[[i]]
    
    nx <- length(vecindin)
    ny <- length(vecindout)
    
    if(nx == 1) {
      myX <- matrix(my_mat[,vecindin], nrow(my_mat), 1)
      cov_1 <- matrix(0, npoints, npoints)
    } else {
      myX <- my_mat[, vecindin]
      cov_1 <- cov(t(myX)) # ((myX - Xn_bar_matrix) %*% t(myX - Xn_bar_matrix)) / nx
    }
    
    if(ny == 1) {
      myY <- matrix(my_mat[,vecindout], nrow(my_mat), 1)
      cov_2 <- matrix(0, npoints, npoints)
    } else {
      myY <- my_mat[, vecindout]
      cov_2 <- cov(t(myY)) # ((myY - Xm_bar_matrix) %*% t(myY - Xm_bar_matrix)) / ny
    }
    
    ##### Initi
    
    X_bar = rowMeans(myX)
    Y_bar = rowMeans(myY)
    
    D <- 1 / (nx + ny - 2) * ((nx - 1) * cov_1  +  (ny - 1) * cov_2)
    D_hotelling <- (nx + ny) / (nx * ny) * D
    
    # Hotelling
    if (plot_eigen){
      temp_svd <- eigen(D_hotelling, symmetric = T)
    } else {
      temp_svd <- rARPACK::eigs_sym(D_hotelling, d)
    }
    
    tau_k <- temp_svd$values
    
    if (plot_eigen){
      cpv <- cpv + c(cumsum(tau_k) / sum(tau_k))
    }
    
    eigen_vec <- temp_svd$vectors
    a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
    rate <- a_k^2 / tau_k
    
    hotelling_k[i] <- sum(rate[1:d])
  }
  
  my_max_h <- max(hotelling_k)
  K <- which.max(hotelling_k)
  vecclus_h <-  c1[[K]]
  
  
  if(plot_eigen) {
    
    #par(oma = c(0, 0, 0, 0), mar = c(3.5, 3, 1, 1), 
    #    mgp = c(2, 1, 0), las = 1)
    plot(1:length(cpv), cpv / nb_combi, type = "o", ...)
    cat("Variance explained in % by the 10 first components: ", 
        round(cpv[1:10] / nb_combi * 100, 2), "\n")
    
  }
  
  
  
  list(stat=my_max_h, vec=vecclus_h)
}
