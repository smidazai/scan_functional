###########################################
simulvec <- function(npoints, shape = "gauss") {
  # initi
  vect <- (0:npoints)/npoints
  k <- 1
  sigma <- 1/((k - 0.5) * pi)
  
  vecyphi <- sqrt(2) * sin(vect/sigma)
  
  if(shape == "gauss") {
    y <- rnorm(1, mean = 0, sd = sigma)
  } else {
    if(shape == "student") { 
      y <- sigma * rt(1, 4)
    } else  {
      if(shape == "chisq") { 
        y <- sigma * rchisq(1, 4)
      } else {
        y <- sigma * rexp(1, 0.5)
      }
    }
  }
  vecY <- y * vecyphi
  exvecY <- vecY
  
  flag <- TRUE
  k <- 2
  while (flag) {
    
    sigma <- 1/((k - 0.5) * pi)
    if(shape == "gauss") {
      y <- rnorm(1, mean = 0, sd = sigma)
    } else {
      if(shape == "student") { 
        y <- sigma * rt(1, 4)
      } else  {
        if(shape == "chisq") { 
        y <- sigma * rchisq(1, 4)
        } else {
          y <- sigma * rexp(1, 4)  
        }
      }
    }
    
    vecyphi <- sqrt(2) * sin(vect/sigma)
    
    exvecY <- vecY
    vecY <- vecY + y * vecyphi
    flag <- (sum((vecY - exvecY) ^ 2) / sum((vecY) ^ 2) > 0.001) #citère d'arrêt. (il dépend de la convergence de notre processus)
    k <- k + 1
    
  }
  vecY
}

###########################################
norm <- function(x) sqrt(sum(x^2))

###########################################
# install.packages("progress")
require("progress")

compute_all <- function(c1, c2, my_mat, my_thresolhd, plot_res = T) {
  # number of combinaison
  nb_combi <- length(c1)
  my_dist <- as(dist(t(my_mat)), "matrix")
  npoints <- nrow(my_mat)
  
  camille <- numeric(nb_combi)
  # dm <- numeric(nb_combi)
  stat_p <- numeric(nb_combi) 
  stat_np <- numeric(nb_combi)
  hotelling_1 <- numeric(nb_combi)
  hotelling_2 <- numeric(nb_combi)
  hotelling_3 <- numeric(nb_combi)
  hotelling_4 <- numeric(nb_combi)
  hotelling_5 <- numeric(nb_combi)
  var_1 <- var_2 <- var_3 <- var_4 <- var_5 <- var_10 <- 
    numeric(nb_combi)
  #   hotelling_10 <- numeric(nb_combi)
  #   hotelling_15 <- numeric(nb_combi)
  hotelling_threshold <- numeric(nb_combi)
  hotelling_all <- numeric(nb_combi)
  
  tau_k_vec <- numeric(nb_combi)
  pb <- progress_bar$new(total = nb_combi)
  
  for(i in 1:nb_combi) {
    pb$tick()
    
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
    Xn_bar_matrix = matrix(X_bar, npoints, nx)
    
    Y_bar = rowMeans(myY)
    Xm_bar_matrix = matrix(Y_bar, npoints, ny)
    
    ##### Horvath
    # Statistic of test
    #U_NM <- sum((X_bar - Y_bar) ^ 2)
    #U_NM = nx * ny / (nx + ny) * U_NM
    
    #den <- ny  / (nx + ny) * sum((myX - Xn_bar_matrix) ^ 2) / nx +
    #  nx  / (nx + ny) * sum((myY - Xm_bar_matrix) ^ 2) / ny 
    
    #horvath_stat_1[i] <- U_NM / den
    #horvath_stat_2[i] <- U_NM / sqrt(den)
    #horvath_stat_3[i] <- U_NM 
    
    # Horvath
    #cov_1 <- cov(myX) # ((myX - Xn_bar_matrix) %*% t(myX - Xn_bar_matrix)) / nx
    #cov_2 <- cov(myY) # ((myY - Xm_bar_matrix) %*% t(myY - Xm_bar_matrix)) / ny
    
    # D <- 1 / (nx + ny) * (ny * cov_1  +  nx * cov_2)
    D <- 1 / (nx + ny - 2) * ((nx - 1) * cov_1  +  (ny - 1) * cov_2)
    D_hotelling <- (nx + ny) / (nx * ny) * D
    
    #temp_svd <- eigen(D)
    #tau_k <- temp_svd$values
    #eigen_vec <- temp_svd$vectors
    
    #tol <- tau_k[1] * sqrt(.Machine$double.eps)
    # select_k <- apply(A_lambda, 2, function(x) all(abs(x)>tol))
    #select_k <- tau_k >= tol
    # tol <- sqrt(.Machine$double.eps) 
    #eigen_vec <- eigen_vec[, select_k]
    #tau_k <- tau_k[select_k]
    
    #a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
    
    
    #horvath_stat_4[i] <- U_NM / max(tau_k)   
    #horvath_stat_5[i] <- nx * ny / (nx + ny) * sum((a_k^2 / tau_k))
    #horvath_stat_6[i] <- nx * ny / (nx + ny) * sum(a_k^2)
    
    # camille
    sd_x <- rowSums((myX - Xn_bar_matrix)^2)
    sd_y <- rowSums((myY - Xm_bar_matrix)^2)
    sd_den <- (sd_x + sd_y) / (nx + ny - 2) * (1/nx + 1/ny)
    camille[i] <- max(abs(X_bar - Y_bar) / sqrt(sd_den), na.rm = T) 
    
    # Hotelling
    temp_svd <- eigen(D_hotelling, symmetric = T)
    
    # rARPACK::eigs_sym(D_hotelling, 10)
    tau_k <- temp_svd$values
    eigen_vec <- temp_svd$vectors
    
    #eigen_vec <- eigen_vec[, cumsum(tau_k)/sum(tau_k) <= 0.9999]
    #tau_k <- tau_k[cumsum(tau_k)/sum(tau_k) <= 0.9999]
    #length_tau_k[i] <- length(tau_k)
    
    a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
    rate <- a_k^2 / tau_k
    hotelling_1[i] <- sum(rate[1])
    hotelling_2[i] <- sum(rate[1:2])
    hotelling_3[i] <- sum(rate[1:3])
    hotelling_4[i] <- sum(rate[1:4])
    hotelling_5[i] <- sum(rate[1:5])
    
    var_1[i] <- sum(tau_k[1:1]) / sum(tau_k)
    var_2[i] <- sum(tau_k[1:2]) / sum(tau_k)
    var_3[i] <- sum(tau_k[1:3]) / sum(tau_k)
    var_4[i] <- sum(tau_k[1:4]) / sum(tau_k)
    var_5[i] <- sum(tau_k[1:5]) / sum(tau_k)
    #  hotelling_10[i] <- sum(rate[1:10])
    var_10[i] <- sum(tau_k[1:10]) / sum(tau_k)     
    #    hotelling_15[i] <- sum(rate[1:15])
    
    eigen_vec <- eigen_vec[, tau_k >= sqrt(.Machine$double.eps)]
    tau_k <- tau_k[tau_k >= sqrt(.Machine$double.eps)]
    hotelling_all[i] <- sum(rate[1:length(tau_k)])
    
    tau_k_vec[i] <- length(which(cumsum(tau_k) / sum(tau_k) < my_thresolhd)) + 1
    eigen_vec <- eigen_vec[, 1:tau_k_vec[i]]
    tau_k <- tau_k[1:tau_k_vec[i]]  
    hotelling_threshold[i] <- sum(rate[1:tau_k_vec[i]])
    
    # Mahalanobis distance
    # dm[i] <- t((X_bar - Y_bar)) %*% MASS::ginv(D) %*% (X_bar - Y_bar)
    #####
    
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
  
  # Among the smallest p-values, we choose the statistic which is the highest
  # ind_p_value_min <- which(p_value == min(p_value))
  # my_max_1 <- max(horvath_stat_1)
  # vecclus_1 <-  cluster_g1[[which(horvath_stat_1 == my_max_1)]]
  # 
  # my_max_2 <- max(horvath_stat_2)
  # vecclus_2 <-  cluster_g1[[which(horvath_stat_2 == my_max_2)]]
  # 
  # my_max_3 <- max(horvath_stat_3)
  # vecclus_3 <-  cluster_g1[[which(horvath_stat_3 == my_max_3)]]
  # 
  # my_max_4 <- max(horvath_stat_4)
  # vecclus_4 <-  cluster_g1[[which(horvath_stat_4 == my_max_4)]]
  # 
  # my_max_5 <- max(horvath_stat_5)
  # vecclus_5 <-  cluster_g1[[which(horvath_stat_5 == my_max_5)]]
  # 
  # my_max_6 <- max(horvath_stat_6)
  # vecclus_6 <-  cluster_g1[[which(horvath_stat_6 == my_max_6)]]
  
  my_max_camille <- max(camille)
  vecclus_camille <-  c1[[which.max(camille)]]
  
  # my_max_dm <- max(dm)
  # vecclus_dm <-  cluster_g1[[which(dm == my_max_dm)]]
  
  my_max_p <- max(stat_p)
  vecclus_p <-  c1[[which.max(stat_p)]]
  
  my_max_np <- max(stat_np)
  vecclus_np <-  c1[[which.max(stat_np)]]
  
  my_max_h_1 <- max(hotelling_1)
  vecclus_h_1 <-  c1[[which.max(hotelling_1)]]
  
  my_max_h_2 <- max(hotelling_2)
  vecclus_h_2 <-  c1[[which.max(hotelling_2)]]
  
  my_max_h_3 <- max(hotelling_3)
  vecclus_h_3 <-  c1[[which.max(hotelling_3)]]
  
  my_max_h_4 <- max(hotelling_4)
  vecclus_h_4 <-  c1[[which.max(hotelling_4)]]
  
  my_max_h_5 <- max(hotelling_5)
  vecclus_h_5 <-  c1[[which.max(hotelling_5)]]
  
  #my_max_h_10 <- max(hotelling_10)
  #vecclus_h_10 <-  c1[[which.max(hotelling_10)]]
  
  my_max_h_threshold <- max(hotelling_threshold)
  vecclus_h_threshold <-  c1[[which.max(hotelling_threshold)]]
  
  my_max_h_all <- max(hotelling_all)
  vecclus_h_all <-  c1[[which.max(hotelling_all)]]
  
  if(plot_res) {
    par(mfrow = c(4, 3), oma = c(0, 0, 0, 0),
        mar = c(3, 3, 1, 1))
    plot(camille, type = "h", main = "camille")
    plot(stat_p, type = "h", main = "p")
    plot(stat_np, type = "h", main = "np")
    plot(hotelling_1, type = "h", main = paste0("hotelling_1: \n",
                                                var_1[which.max(hotelling_1)]))
    plot(hotelling_2, type = "h", main = paste0("hotelling_2: \n",
                                                var_2[which.max(hotelling_2)]))
    plot(hotelling_3, type = "h", main = paste0("hotelling_3: \n",
                                                var_3[which.max(hotelling_3)]))
    plot(hotelling_4, type = "h", main = paste0("hotelling_4: \n",
                                                var_4[which.max(hotelling_4)]))
    plot(hotelling_5, type = "h", main = paste0("hotelling_5: \n",
                                                var_5[which.max(hotelling_5)]))
    # plot(hotelling_10, type = "h", main = paste0("hotelling_10: \n",
    #                                              var_10[which.max(hotelling_10)]))
    plot(hotelling_all, type = "h", main = "All")
    plot(hotelling_threshold, type = "h", 
         main = paste0("hotelling best = \n",
                       tau_k_vec[which.max(hotelling_threshold)]))
    
  }
  
  list(
    camille = list(stat=my_max_camille, vec=vecclus_camille),
    #       dm = list(stat=my_max_dm, vec=vecclus_dm),
    stat_p = list(stat=my_max_p, vec=vecclus_p),
    stat_np = list(stat=my_max_np, vec=vecclus_np),
    hotelling_1 = list(stat=my_max_h_1, vec=vecclus_h_1),
    hotelling_2 = list(stat=my_max_h_2, vec=vecclus_h_2),
    hotelling_3 = list(stat=my_max_h_3, vec=vecclus_h_3),
    hotelling_4 = list(stat=my_max_h_4, vec=vecclus_h_4),
    hotelling_5 = list(stat=my_max_h_5, vec=vecclus_h_5),
    hotelling_threshold = list(stat=my_max_h_threshold, 
                               vec=vecclus_h_threshold,
                               k = tau_k_vec[which.max(hotelling_threshold)]),
    #     hotelling_15 = list(stat=my_max_h_15, vec=vecclus_h_15),
    hotelling_all = list(stat=my_max_h_all, vec=vecclus_h_all)
  )
}

##### DFFSS ####### 
compute_dffss <- function(c1, c2, my_mat) {
  
  # number of combinaison
  nb_combi <- length(c1)
  camille <- numeric(nb_combi)
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
    camille[i] <- max(abs(X_bar - Y_bar) / sqrt(sd_den), na.rm = T) 
  }
  
  my_max_camille <- max(camille)
  vecclus_camille <-  c1[[which.max(camille)]]
  
  list(
    stat=my_max_camille, 
    vec=vecclus_camille
  )
}

##### parametric ###

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
  
  list(
    stat=my_max_p, vec=vecclus_p
  )
}

######. Non parametric

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

######. Hotelling

compute_h <- function(c1, c2, my_mat, d = 5,
                              plot_eigen = F) {
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
    eigen_vec <- temp_svd$vectors
    
    a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
    rate <- a_k^2 / tau_k
    
    if (plot_eigen){
      cpv <- cpv + c(cumsum(tau_k) / sum(tau_k))
    }
    hotelling_k[i] <- sum(rate[1:d])
  }

  my_max_h <- max(hotelling_k)
  K <- which.max(hotelling_k)
  vecclus_h <-  c1[[K]]
  
  
  if(plot_eigen) {
    
    #par(oma = c(0, 0, 0, 0), mar = c(3.5, 3, 1, 1), 
    #    mgp = c(2, 1, 0), las = 1)
    plot(1:length(cpv), cpv / nb_combi, type = "o",
         xlab = "k", 
         ylab = TeX("$\\bar{CPV}$"), ylim = c(0, 1), cex = 0.5)
    cat("Variance explained in % by the 10 first components: ", 
        round(cpv[1:10] / nb_combi * 100, 2), "\n")
    
  }
    

  
  list(stat=my_max_h, vec=vecclus_h)
}

###########################################
statscan_all_tibo <- function(X, my_pairs, mini = 2, 
                               maxi = trunc(nrow(X) / 2),
                               my_thresolhd = 0.995,
                               second_cluster = F) {
  

  # initialization
  cluster_g1 <- my_pairs[[1]] 
  cluster_g2 <- my_pairs[[2]] 
  npoints <- nrow(X)
  
  # restriction on the size of the cluster 
  size_group <- sapply(cluster_g1, length)
  ind_restriction <- size_group >= mini & size_group <= maxi
  cluster_g1 <- cluster_g1[ind_restriction]
  cluster_g2 <- cluster_g2[ind_restriction] 
  
  # First Clusters
  res_1 <- compute_all(cluster_g1, cluster_g2, X, my_thresolhd)
  
  vecclus_camille <-  res_1$camille$vec
  vecclus_p <-  res_1$stat_p$vec
  vecclus_np <- res_1$stat_np$vec
  vecclus_h_1 <- res_1$hotelling_1$vec  
  vecclus_h_2 <- res_1$hotelling_2$vec  
  vecclus_h_3 <- res_1$hotelling_3$vec  
  vecclus_h_4 <-  res_1$hotelling_4$vec  
  vecclus_h_5 <-  res_1$hotelling_5$vec  
  vecclus_h_threshold <-  res_1$hotelling_threshold$vec  
#  vecclus_h_15 <-  res_1$hotelling_15$vec  
  vecclus_h_all <- res_1$hotelling_all$vec   
    
  if(second_cluster) {
    # camille 
    cluster_g1_temp <- sapply(cluster_g1, function(x) setdiff(x, vecclus_camille))
    cluster_g2_temp <- sapply(cluster_g2, function(x) setdiff(x, vecclus_camille))
    id_pos <- union(which(sapply(cluster_g1_temp, function(x) length(x) == 0)),
          which(sapply(cluster_g2_temp, function(x) length(x) == 0)))
    temp <- compute_camille(cluster_g1_temp[-id_pos], 
                           cluster_g2_temp[-id_pos], X)
    vecclus_camille_2 <- temp$vec 
    stat_camille_2 <- temp$stat 
    
    # p
    cluster_g1_temp <- sapply(cluster_g1, function(x) setdiff(x, vecclus_p))
    cluster_g2_temp <- sapply(cluster_g2, function(x) setdiff(x, vecclus_p))
    id_pos <- union(which(sapply(cluster_g1_temp, function(x) length(x) == 0)),
                    which(sapply(cluster_g2_temp, function(x) length(x) == 0)))
    temp <- compute_p(cluster_g1_temp[-id_pos], 
              cluster_g2_temp[-id_pos],
              X)
    vecclus_p_2 <- temp$vec 
    stat_p_2 <- temp$stat
    
    # np
    cluster_g1_temp <- sapply(cluster_g1, function(x) setdiff(x, vecclus_np))
    cluster_g2_temp <- sapply(cluster_g2, function(x) setdiff(x, vecclus_np))
    id_pos <- union(which(sapply(cluster_g1_temp, function(x) length(x) == 0)),
                    which(sapply(cluster_g2_temp, function(x) length(x) == 0)))
    temp <- compute_p(cluster_g1_temp[-id_pos], cluster_g2_temp[-id_pos],
                      X)
    vecclus_np_2 <- temp$vec 
    stat_np_2 <- temp$stat 
    
    # hotelling Threshold
    cluster_g1_temp <- sapply(cluster_g1, function(x) setdiff(x, vecclus_h_threshold))
    cluster_g2_temp <- sapply(cluster_g2, function(x) setdiff(x, vecclus_h_threshold))
    id_pos <- union(which(sapply(cluster_g1_temp, function(x) length(x) == 0)),
                    which(sapply(cluster_g2_temp, function(x) length(x) == 0)))
    temp <- compute_hotelling(cluster_g1_temp[-id_pos], cluster_g2_temp[-id_pos],
                              X, k = res_1$hotelling_threshold$k)
    vecclus_h_threshold_2 <- temp$vec
    stat_h_threshold_2 <- temp$stat
  } else {
    stat_camille_2 <- NULL
    stat_p_2 <- NULL
    stat_np_2 <- NULL
    stat_h_threshold_2 <- NULL
    
    vecclus_camille_2 <- NULL
    vecclus_p_2 <- NULL
    vecclus_np_2 <- NULL
    vecclus_h_threshold_2 <- NULL
  }
  
  list(
    camille = list(stat=res_1$camille$stat, vec=vecclus_camille, 
                   stat2 = stat_camille_2, vec2 = vecclus_camille_2),
    #       dm = list(stat=my_max_dm, vec=vecclus_dm),
    stat_p = list(stat=res_1$stat_p$stat, vec=vecclus_p, 
                  stat2 = stat_p_2, vec2 = vecclus_p_2),
    stat_np = list(stat=res_1$stat_np$stat, vec=vecclus_np, 
                   stat2 = stat_np_2, vec2 = vecclus_np_2),
    # hotelling_1 = list(stat=my_max_h_1, vec=vecclus_h_1),
    # hotelling_2 = list(stat=my_max_h_2, vec=vecclus_h_2),
    # hotelling_3 = list(stat=my_max_h_3, vec=vecclus_h_3),
    # hotelling_4 = list(stat=my_max_h_4, vec=vecclus_h_4),
    hotelling_t = list(stat=res_1$hotelling_threshold$stat, 
                       vec=vecclus_h_threshold, 
                       stat2 = stat_h_threshold_2,
                       vec2 = vecclus_h_threshold_2)
    # hotelling_10 = list(stat=my_max_h_10, vec=vecclus_h_10),
    # hotelling_15 = list(stat=my_max_h_15, vec=vecclus_h_15)
  )
  
}


##################################################

find_all_cluster <- function (Matcoord) {
  n <- nrow(Matcoord)
  Matdist <- as.matrix(dist(Matcoord, upper = TRUE))
  vecord_list <- vector("list", n)
  for (k in 1:n) {
    vecord_list[[k]] <- order(Matdist[, k])
  }
  res_cluster_g1 <- vector("list", 0)
  res_cluster_g2 <- vector("list", 0)
  
  matrix_g1 <- vector("list", n)
#  matrix_g2 <- vector("list", n)
  for(k in 1:n) {
    matrix_g1[[k]] <- rep(0, n)
 #   matrix_g2[[k]] <- rep(0, n)     
  }
  
  nb_combi <- 0
  for (k in 1:(n-1)) {
    for (j in 1:n) {
      temp_1 <- vecord_list[[j]][1:k]
      temp_2 <- vecord_list[[j]][(k+1):n]
      #cond_1 <- sapply(res_cluster_g1, function(x) all(temp_1 %in% x) & 
      #                   length(temp_1) == length(x))
      #cond_2 <- sapply(res_cluster_g1, function(x) all(temp_2 %in% x) & 
      #                   length(temp_2) == length(x))
      my_vec <- my_vec_2 <-numeric(n)
      my_vec[temp_1] <- 1
      my_vec_2[temp_2] <- 1
      # my_length <- k
      cond_1 <-  any(matrix_g1[[k]] %*% my_vec == k)
      cond_2 <-  any(matrix_g1[[n-k]] %*% my_vec_2 == n-k)
        
      if (!(any(cond_1) | any(cond_2))) {
        nb_combi <- nb_combi + 1
        res_cluster_g1[[nb_combi]] <- temp_1
        res_cluster_g2[[nb_combi]] <- temp_2
        matrix_g1[[k]] <- rbind(matrix_g1[[k]], my_vec)
    #    matrix_g2[[n-k]] <- rbind(matrix_g2[[n-k]], rep(1, n) - my_vec)
      }
    }
  }
  cat("Number of possible combinaison: ", nb_combi, "\n")
  return(list(vec_g1 = res_cluster_g1,
              vec_g2 = res_cluster_g2))
}

##################################################

power_simu <- function(input, niter = 99, nsimu = 100) {
  # init
  sizeclust <- parms_df[input, "sizeclust"] 
  alpha <- parms_df[input, "alpha"] 
  shape <- parms_df[input, "shape"] 
  type_shift <- parms_df[input, "type_shift"] 
  
  nb_est <- 10
  power <- nTP <- nFP <- c(
    camille = 0,
    p = 0,
    np = 0,
    hotelling_1 = 0,
    hotelling_2 = 0,
    hotelling_3 = 0,
    hotelling_4 = 0,
    hotelling_5 = 0,
    hotelling_10 = 0,
    hotelling_15 = 0
   ) 
  

  # parameters of simulation 
  nobs <- nrow(Matcoord)
  npoints <- 100
  ndrop <- 100
  t.disc <- (1:(npoints)) / (npoints) # (1:(npoints+ndrop)) / (npoints+ndrop)
  
  vecclus <- c(74, 92, 91, 93, 77, 90, 94, 76, 59, 27, 79, 44, 26, 75, 2, 61)[1:sizeclust]
  nclus <- length(vecclus)
  veccluster <- numeric(nobs)
  veccluster[vecclus] <- 1
  
  # restriction on the size of the cluster 
  mini <- 2 
  maxi <- trunc(nobs / 2)
  my_thresolhd <- 0.995
  cluster_g1 <- my_pairs[[1]] 
  cluster_g2 <- my_pairs[[2]] 
  size_group <- sapply(cluster_g1, length)
  ind_restriction <- size_group >= mini & size_group <= maxi
  cluster_g1 <- cluster_g1[ind_restriction]
  cluster_g2 <- cluster_g2[ind_restriction] 
  # number of combinaison
  nb_combi <- length(cluster_g1)
  
  # simulations 
  i_simu <- 1
  while (i_simu <= nsimu) {
    print(i_simu)
    X <- matrix(0, npoints+ndrop, nobs)
    for (k in 1:nobs) {
      X[, k] <- simulvec(npoints+(ndrop-1), shape = shape) 
      if(type_shift == 1) {
        X[(ndrop+1):nrow(X), k] <- X[(ndrop+1):nrow(X), k] + alpha * (veccluster[k] == 1)
        } else {
          if (type_shift == 2) {
            X[(ndrop+1):nrow(X), k] <- X[(ndrop+1):nrow(X), k] + alpha * t.disc * (veccluster[k] == 1)
            } else {
              if (type_shift == 3) {
                X[(ndrop+1):nrow(X), k] <- X[(ndrop+1):nrow(X), k] + alpha * t.disc * (1 - t.disc) * (veccluster[k] == 1)
                } else {
                  X[(ndrop+1):nrow(X), k] <- X[(ndrop+1):nrow(X), k] + alpha * exp(-100 * (t.disc - 0.5) ^ 2) / 3 * (veccluster[k] == 1)
                }
            }
        }
    }
    # drop first observation 
    X <- X[-(1:ndrop), ]
  
    # initialization
    npoints <- nrow(X)
    my_dist <- as(dist(t(X)), "matrix")
  
    # observed stat
    camille <- numeric(nb_combi)
    stat_p <- numeric(nb_combi) 
    stat_np <- numeric(nb_combi)
    hotelling_1 <- numeric(nb_combi)
    hotelling_2 <- numeric(nb_combi)
    hotelling_3 <- numeric(nb_combi)
    hotelling_4 <- numeric(nb_combi)
    hotelling_5 <- numeric(nb_combi)
    hotelling_10 <- numeric(nb_combi)
    hotelling_15 <- numeric(nb_combi)
    
    cpv <- numeric(npoints)
    
    for (i in (1:nb_combi)) {
      
      vecindin <- cluster_g1[[i]]
      vecindout <- cluster_g2[[i]]
    
      nx <- length(vecindin)
      ny <- length(vecindout)
    
      if(nx == 1) {
        myX <- matrix(X[,vecindin], nrow(X), 1)
        cov_1 <- matrix(0, npoints, npoints)
        } else {
          myX <- X[, vecindin]
          cov_1 <- cov(t(myX)) 
        }
      if(ny == 1) {
        myY <- matrix(X[,vecindout], nrow(X), 1)
        cov_2 <- matrix(0, npoints, npoints)
        } else {
          myY <- X[, vecindout]
          cov_2 <- cov(t(myY)) # 
        }
      ##### Initi
      X_bar = rowMeans(myX)
      Xn_bar_matrix = matrix(X_bar, npoints, nx)
    
      Y_bar = rowMeans(myY)
      Xm_bar_matrix = matrix(Y_bar, npoints, ny)
    
      D <- 1 / (nx + ny - 2) * ((nx - 1) * cov_1  +  (ny - 1) * cov_2)
      D_hotelling <- (nx + ny) / (nx * ny) * D
    
      # camille
      sd_x <- rowSums((myX - Xn_bar_matrix)^2)
      sd_y <- rowSums((myY - Xm_bar_matrix)^2)
      sd_den <- (sd_x + sd_y) / (nx + ny - 2) * (1/nx + 1/ny)
      camille[i] <- max(abs(X_bar - Y_bar) / sqrt(sd_den)) 
    
      # Hotelling
      #temp_svd <- eigen(D_hotelling, symmetric = T)
      #tau_k <- temp_svd$values
      #eigen_vec <- temp_svd$vectors
      #tol <- tau_k[1] * sqrt(.Machine$double.eps)
      #eigen_vec <- eigen_vec[, tau_k >= tol]
      #tau_k <- tau_k[tau_k >= tol]
      
      #eigen_vec <- eigen_vec[, cumsum(tau_k) / sum(tau_k) < my_thresolhd]
      #tau_k <- tau_k[cumsum(tau_k) / sum(tau_k) < my_thresolhd]  
      #vec_tau[i] = length(tau_k)
      temp_svd <- RSpectra::eigs_sym(D_hotelling, 15)
      tau_k <- temp_svd$values
      eigen_vec <- temp_svd$vectors
      
      a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
      rate_value <- a_k^2 / tau_k 
      hotelling_1[i] <- sum(rate_value[1])
      hotelling_2[i] <- sum(rate_value[1:2])
      hotelling_3[i] <- sum(rate_value[1:3])
      hotelling_4[i] <- sum(rate_value[1:4])
      hotelling_5[i] <- sum(rate_value[1:5])
      hotelling_10[i] <- sum(rate_value[1:10])
      hotelling_15[i] <- sum(rate_value[1:15])
      
      ##### Parametric
      moy_glob <- ((nx * X_bar) + (ny * Y_bar))/(nx + ny)
    
      diff1 <- X_bar - moy_glob
      diff2 <- Y_bar - moy_glob
      nor1 <- norm(diff1)^2
      nor2 <- norm(diff2)^2
      numerateur <- (nx * nor1) + (ny * nor2)
      denominateur <- (1/((nx + ny)-2))*(sum(sd_x) + sum(sd_y))
      stat_p[i] <- numerateur / denominateur
    
      ##### Non Parametric
      vec_int <- numeric(npoints)
       
      for(j in 1:nx)
        vec_int <- vec_int + X[, vecindin[j]] * sum(1/my_dist[vecindin[j], vecindout])
       
      for(j in 1:ny)
         vec_int <- vec_int - X[, vecindout[j]] * sum(1/my_dist[vecindout[j], vecindin])
       
      vec_int <- vec_int / (nx * ny)
      temp_stat_np <- norm(vec_int)
       
      stat_np[i] <-  sqrt((nx * ny) / (nx + ny)) * temp_stat_np
    }
    
    my_stat <- c(
      camille = max(camille),
      p = max(stat_p),
      np = max(stat_np),
      hotelling_1 = max(hotelling_1),
      hotelling_2 = max(hotelling_2),
      hotelling_3 = max(hotelling_3),
      hotelling_4 = max(hotelling_4),
      hotelling_5 = max(hotelling_5),
      hotelling_10 = max(hotelling_10),
      hotelling_15 = max(hotelling_15)
      )
    
    
    clus_found <- list(
      camille =  cluster_g1[[which.max(camille)]],
      p = cluster_g1[[which.max(stat_p)]],
      np = cluster_g1[[which.max(stat_np)]],
      hotelling_1 = cluster_g1[[which.max(hotelling_1)]],
      hotelling_2 = cluster_g1[[which.max(hotelling_2)]],
      hotelling_3 = cluster_g1[[which.max(hotelling_3)]],
      hotelling_4 = cluster_g1[[which.max(hotelling_4)]],
      hotelling_5 = cluster_g1[[which.max(hotelling_5)]],
      hotelling_10 = cluster_g1[[which.max(hotelling_10)]],
      hotelling_15 = cluster_g1[[which.max(hotelling_15)]]
      )
    
    # correct estimation  
    intersect_t <- sapply(clus_found, function(x) length(intersect(x, vecclus)))  
    setdiff_t <- sapply(clus_found, function(x) length(setdiff(x, vecclus))) 
  
    ################################. Replication by permutation 
    # initialisation result
    my_pvalue <- integer(nb_est) 
    stat_max <- numeric(niter)

    for (b in 1:niter) {
      print(b)
      #set.seed(b)
      perm <- sample(nobs)
      X_sim <- X[, perm]
      my_dist_sim <- my_dist[perm, perm]

      camille <- numeric(nb_combi)
      stat_p <- numeric(nb_combi) 
      stat_np <- numeric(nb_combi)
      hotelling_1 <- numeric(nb_combi)
      hotelling_2 <- numeric(nb_combi)
      hotelling_3 <- numeric(nb_combi)
      hotelling_4 <- numeric(nb_combi)
      hotelling_5 <- numeric(nb_combi)
      hotelling_10 <- numeric(nb_combi)
      hotelling_15 <- numeric(nb_combi)
    
      for (i in (1:nb_combi)) {
        vecindin <- cluster_g1[[i]]
        vecindout <- cluster_g2[[i]]
        nx <- length(vecindin)
        ny <- length(vecindout)
        
        if (nx == 1) {
          myX <- matrix(X_sim[,vecindin], nrow(X_sim), 1)
         cov_1 <- matrix(0, npoints, npoints)
          } else {
            myX <- X_sim[, vecindin]
            cov_1 <- cov(t(myX)) 
          }
        
        if(ny == 1) {
          myY <- matrix(X_sim[,vecindout], nrow(X_sim), 1)
          cov_2 <- matrix(0, npoints, npoints)
          } else {
            myY <- X_sim[, vecindout]
            cov_2 <- cov(t(myY)) # 
          }
        
        ##### Initi
        X_bar = rowMeans(myX)
        Xn_bar_matrix = matrix(X_bar, npoints, nx)
      
        Y_bar = rowMeans(myY)
        Xm_bar_matrix = matrix(Y_bar, npoints, ny)
      
        D <- 1 / (nx + ny - 2) * ((nx - 1) * cov_1  +  (ny - 1) * cov_2)
        D_hotelling <- (nx + ny) / (nx * ny) * D
  
        # camille
        sd_x <- rowSums((myX - Xn_bar_matrix)^2)
        sd_y <- rowSums((myY - Xm_bar_matrix)^2)
        sd_den <- (sd_x + sd_y) / (nx + ny - 2) * (1/nx + 1/ny)
        camille[i] <- max(abs(X_bar - Y_bar) / sqrt(sd_den)) 
      
        # Hotelling
        #temp_svd <- eigen(D_hotelling, symmetric = T)
      
        #tau_k <- temp_svd$values
        #eigen_vec <- temp_svd$vectors
      
        #tol <- tau_k[1] * sqrt(.Machine$double.eps)
        #eigen_vec <- eigen_vec[, tau_k >= tol]
        #tau_k <- tau_k[tau_k >= tol]

        #eigen_vec <- eigen_vec[, cumsum(tau_k) / sum(tau_k) < my_thresolhd]
        #tau_k <- tau_k[cumsum(tau_k) / sum(tau_k) < my_thresolhd] 
        #vec_tau[i] = length(tau_k)
        
        temp_svd <- RSpectra::eigs_sym(D_hotelling, 15)
        tau_k <- temp_svd$values
        eigen_vec <- temp_svd$vectors
        
        a_k <- t(eigen_vec) %*% (X_bar - Y_bar)
        rate_value <- a_k^2 / tau_k 
        hotelling_1[i] <- sum(rate_value[1])
        hotelling_2[i] <- sum(rate_value[1:2])
        hotelling_3[i] <- sum(rate_value[1:3])
        hotelling_4[i] <- sum(rate_value[1:4])
        hotelling_5[i] <- sum(rate_value[1:5])
        hotelling_10[i] <- sum(rate_value[1:10])
        hotelling_15[i] <- sum(rate_value[1:15])
      
        ##### Parametric
        moy_glob <- ((nx * X_bar) + (ny * Y_bar))/(nx + ny)
        diff1 <- X_bar - moy_glob
        diff2 <- Y_bar - moy_glob
        nor1 <- norm(diff1)^2
        nor2 <- norm(diff2)^2
        numerateur <- (nx * nor1) + (ny * nor2)
        denominateur <- (1/((nx + ny)-2))*(sum(sd_x) + sum(sd_y))
        stat_p[i] <- numerateur / denominateur
      
        ##### Non Parametric
        vec_int <- numeric(npoints)
      
        for(j in 1:nx)
          vec_int <- vec_int + X_sim[, vecindin[j]] * sum(1/my_dist_sim[vecindin[j], vecindout])
      
        for(j in 1:ny)
          vec_int <- vec_int - X_sim[, vecindout[j]] * sum(1/my_dist_sim[vecindout[j], vecindin])
      
        vec_int <- vec_int / (nx * ny)
        temp_stat_np <- norm(vec_int)
      
        stat_np[i] <-  sqrt((nx * ny) / (nx + ny)) * temp_stat_np
      }
      
      my_stat_perm <- c(
        camille = max(camille),
        p = max(stat_p),
        np = max(stat_np),
        hotelling_1 = max(hotelling_1),
        hotelling_2 = max(hotelling_2),
        hotelling_3 = max(hotelling_3),
        hotelling_4 = max(hotelling_4),
        hotelling_5 = max(hotelling_5),
        hotelling_10 = max(hotelling_10),
        hotelling_15 = max(hotelling_15)
        )
      
      # par(mfrow = c(2, 2))
      # plot(hotelling_1, type = "h")
      # abline(h = my_stat["hotelling_1"], col = "red")
      # plot(hotelling_5, type = "h")
      # abline(h = my_stat["hotelling_5"], col = "red")
      # plot(hotelling_10, type = "h")
      # abline(h = my_stat["hotelling_10"], col = "red")
      # plot(hotelling_15, type = "h")
      # abline(h = my_stat["hotelling_15"], col = "red")
      my_pvalue <- my_pvalue + (my_stat < my_stat_perm)
    }
    
    my_pvalue <- my_pvalue / (niter + 1)
  
    power[my_pvalue <= 0.05] <- power[my_pvalue <= 0.05] + 1
    nTP[my_pvalue <= 0.05] <- nTP[my_pvalue <= 0.05] + intersect_t[my_pvalue <= 0.05]
    nFP[my_pvalue <= 0.05] <- nFP[my_pvalue <= 0.05] + setdiff_t[my_pvalue <= 0.05]
  
    i_simu <- i_simu + 1
  }
  
  return(list(power = power,
            nTP = nTP,
            nFP = nFP))
}


###########################

signscan <- function(Matcoord, my_pairs, sizeclust = 8, 
                     alpha = 1, 
                     type_shift = 1, 
                     shape = "gauss",
                     nsim = 1000, niter = 999, 
                     nb_cluster = 10) {
  
  nb_est <- 8
  power <- integer(nb_est)
  nTP <- integer(nb_est)
  nFP <- integer(nb_est)
  
  # parameters of simulation 
  nobs <- nrow(Matcoord)
  npoints <- 100

  vecclus <- c(74, 92, 91, 93, 77, 90, 94, 76, 59, 27, 79, 44, 26, 75, 2, 61)[1:sizeclust]
  nclus <- length(vecclus)
  veccluster <- numeric(nobs)
  veccluster[vecclus] <- 1
  
  #####
  
  for (sim in (1:nsim)) {
    
    print(paste("## Simulation =", sim))
    
    X <- matrix(0, npoints+1, nobs)
    for (k in 1:nobs) {
      X[, k] <- simulvec(npoints, shape = shape) 
      if(type_shift == 1) {
        X[, k] <- X[, k] + alpha * (veccluster[k] == 1)
      } else {
        if (type_shift == 2) {
          t.disc <- (1:(npoints+1)) / npoints
          X[, k] <- X[, k] + alpha * t.disc * (veccluster[k] == 1)
        } else {
           if (type_shift == 3) {
             t.disc <- (1:(npoints+1)) / npoints
             X[, k] <- X[, k] + alpha * t.disc * (1 - t.disc) * (veccluster[k] == 1)
           } else {
             t.disc <- (1:(npoints+1)) / npoints
             X[, k] <- X[, k] + alpha * exp(-100 * (t.disc - 0.5) ^ 2) / 3 * (veccluster[k] == 1)
           }
        }
      }
    }
    # drop first observation  
    X <- X[-1, ]
    
    # observed stat
    temp <- statscan_all_tibo(X, my_pairs)
    my_stat <- sapply(temp, function(x) x$stat)
    clus_found <- lapply(temp, function(x) x$vec)
    
    # correct estimation  
    intersect_t <- sapply(clus_found, function(x) length(intersect(x, vecclus)))  
    setdiff_t <- sapply(clus_found, function(x) length(setdiff(x, vecclus))) 
    
    # initialisation result
    my_pvalue <- integer(nb_est) 
    my_stat_iter <- vector("list", nb_est)
    
    statscan_all_par <- function(vec_b, true_value, nb_est) {
      
      nobs <- ncol(X)
      res_loc <- numeric(nb_est)
      
      for(b in vec_b) {
        set.seed(b)
        perm <- sample(nobs)
        MatXsim <- X[, perm]
        temp <- statscan_all_tibo(MatXsim, my_pairs)
        res_loc <- res_loc + (true_value < sapply(temp, function(x) x$stat))
      }
      
      return(res_loc)
    } 
    
    compute_per_cluster <- niter / nb_cluster
      
    split_vector <- split(1:niter, ceiling(seq_along(1:niter) / compute_per_cluster))
    require(parallel)
    cl <- makeCluster(nb_cluster)
    clusterExport(cl, c("statscan_all_tibo", "norm", "X", "my_pairs"))
    system.time({
      res_par <- clusterApplyLB(cl, split_vector, statscan_all_par, 
                                true_value = my_stat, nb_est = nb_est) 
    })
    stopCluster(cl)
    
    # compute p_values
    vec_p_value <- integer(nb_est)
    for(b in 1:nb_cluster)
      vec_p_value <- vec_p_value + res_par[[b]]
    
    vec_p_value <- vec_p_value / (niter + 1)
    
    power[vec_p_value <= 0.05] <- power[vec_p_value <= 0.05] + 1 
    nTP[vec_p_value <= 0.05] <- nTP[vec_p_value <= 0.05] + intersect_t[vec_p_value <= 0.05]
    nFP[vec_p_value <= 0.05] <- nFP[vec_p_value <= 0.05] + setdiff_t[vec_p_value <= 0.05]
    
  }
  
  
  return(list(power = power,
              nTP = nTP,
              nFP = nFP))
  
}


