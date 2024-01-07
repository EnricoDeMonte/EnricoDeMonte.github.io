# Functions to conduct for "Nonparametric instrumental regression with two-way fixed effects"
# De Monte (2023), Journal of Econometric Methods.  

# October 2023

require(MASS)
LMU.CVCMB = function( h, 
                      x, 
                      y, 
                      N,
                      T, 
                      effects = c("individual", "time", "indiv-time"),
                      type = c("ll", "lq", "l_cubic") , 
                      jump,
                      c1, c2 ){
  NT = N*T
  mx = matrix(0, nrow = NT, ncol = 1)   # mx = conditional mean vector
  y = as.vector(y)
  x = as.vector(x)
  
  if(type == "ll"){
    Dx = matrix(0, nrow = NT, ncol = 1)  # Dx = 1st gradients vector 
    colnames(Dx) <- c("gradient.1")
    P = 1
  }
  
  if(type == "lq"){
    Dx = matrix(0, nrow = NT, ncol = 2)  # Dx = 2nd gradients vector 
    colnames(Dx) <- c("gradient.1", "gradient.2")
    P = 2 
  }
  
  if(type == "l_cubic"){
    Dx = matrix(0, nrow = NT, ncol = 3)  # Dx = 3rd gradients vecor 
    colnames(Dx) <- c("gradient.1", "gradient.2", "gradient.3")
    P = 3
  }
  
  
  i <- 1 # leave-one out while-loop 
  while(i <= NT) {
    
    if(effects =="individual"){  
      kxw_sum_0 = matrix(NA, nrow = NT, ncol = 1)
      kx_sum_0 = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_0 = matrix(NA, nrow = NT, ncol = 1)
      
      kxw_sum = matrix(NA, nrow = NT, ncol = 1)
      kx_sum = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_x = matrix(NA, nrow = NT, ncol = 1)
    }
    
    if(effects == "time") {
      kxw_sum_0 = matrix(NA, nrow = NT, ncol = 1)
      kx_sum_0 = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_0 = matrix(NA, nrow = NT, ncol = 1)
      
      kxw_sum = matrix(NA, nrow = NT, ncol = 1)
      kx_sum = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_x = matrix(NA, nrow = NT, ncol = 1)
    } 
    
    if(effects == "indiv-time"){
      
      kxw_sum_0_a = matrix(NA, nrow = NT, ncol = 1)
      kxw_sum_0_b = matrix(NA, nrow = NT, ncol = 1)
      kxw_sum_0_c = matrix(NA, nrow = NT, ncol = 1)
      
      kxw_sum_a = matrix(NA, nrow = NT, ncol = 1)
      kxw_sum_b = matrix(NA, nrow = NT, ncol = 1)
      kxw_sum_c = matrix(NA, nrow = NT, ncol = 1)
      
      gamma_it_0_a = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_0_b = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_0_c = matrix(NA, nrow = NT, ncol = 1)
      
      gamma_it_x_a = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_x_b = matrix(NA, nrow = NT, ncol = 1)
      gamma_it_x_c = matrix(NA, nrow = NT, ncol = 1)
      
      
    } 
    
    X_star_x = matrix(NA, nrow = NT, ncol = P)
    X_star_0 = matrix(NA, nrow = NT, ncol = P)
    
    Y_star_x = matrix(NA, nrow = NT, ncol = 1)
    Y_star_0 = matrix(NA, nrow = NT, ncol = 1) 
    
    # Calculates weights and transformed data 
    # 1.) Kernel weights for x and 0 
    
    # leave-one out, by putting NA for the i th left-out observation
    
    kx = dnorm((x - x[i])/h, 0, 1 )
    kx_0 = dnorm((x - 0)/h, 0, 1 )
    
    kxw = kx
    kxw[i] = NA
    
    kxw_0 = kx_0
    kxw_0[i] = NA
    
    # loop to sum of given i over all each t's    
    if(effects == "individual"){ 
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        kxw_sum[n_index,1] = sum(kxw[n_index_loo])
        kxw_sum_0[n_index,1] = sum(kxw_0[n_index_loo])  # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
      }
    }
    if(effects == "time"){
      for(t in 1:T){
        # t_index holds the indeces for of all i's for a given t
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        if(i %in% t_index){ # remove the index of the left-out observation
          t_index_loo = t_index[-which(t_index == i)]
        }else{
          t_index_loo = t_index
        }
        kxw_sum[t_index] = sum(kxw[t_index_loo]) # calculate sum of given t over all each i's 
        kxw_sum_0[t_index] = sum(kxw_0[t_index_loo]) 
      }
    }  
    if(effects == "indiv-time"){
      # individual 
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        kxw_sum_a[n_index,1] = sum(kxw[n_index_loo])
        kxw_sum_0_a[n_index,1] = sum(kxw_0[n_index_loo])  # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
      }
      # time 
      for(t in 1:T){
        # t_index holds the indeces for of all i's for a given t
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        if(i %in% t_index){ # remove the index of the left-out observation
          t_index_loo = t_index[-which(t_index == i)]
        }else{
          t_index_loo = t_index
        }
        kxw_sum_b[t_index] = sum(kxw[t_index_loo]) # calculate sum of given t over all each i's 
        kxw_sum_0_b[t_index] = sum(kxw_0[t_index_loo]) 
      }
      
      kxw_sum_c[,1] = sum(kxw[-i]) # calculate sum of given t over all each i's 
      kxw_sum_0_c[,1] = sum(kxw_0[-i]) 
      
    }  
    sum(is.na(kxw_sum_b))
    if(effects == "indiv-time"){
      # Add constants c1 and C2
      kxw_sum_a = ifelse((is.na(kxw_sum_a) == FALSE), (kxw_sum_a + c1), NA)
      kxw_sum_b = ifelse((is.na(kxw_sum_b) == FALSE), (kxw_sum_b + c1), NA)
      kxw_sum_c = ifelse((is.na(kxw_sum_c) == FALSE), (kxw_sum_c + c1), NA)
      
      kxw_sum_0_a = ifelse((is.na(kxw_sum_0_a) == FALSE), (kxw_sum_0_a + c1), NA)
      kxw_sum_0_b = ifelse((is.na(kxw_sum_0_b) == FALSE), (kxw_sum_0_b + c1), NA)
      kxw_sum_0_c = ifelse((is.na(kxw_sum_0_c) == FALSE), (kxw_sum_0_c + c1), NA)
      
      wx_a = kxw/(kxw_sum_a )    # weights (Lee (2018) p. 4 equation (5))
      wx_b = kxw/(kxw_sum_b )
      wx_c = kxw/(kxw_sum_c )
      # sum(is.na(kxw)) ; sum(is.na(kxw_sum_a)) ; sum(is.na(wx_b))
      wx_0_a = kxw_0/(kxw_sum_0_a )
      wx_0_b = kxw_0/(kxw_sum_0_b )
      wx_0_c = kxw_0/(kxw_sum_0_c )
      
    }else{
      kxw_sum = ifelse((is.na(kxw_sum) == FALSE), (kxw_sum + c1), NA)
      kxw_sum_0 = ifelse((is.na(kxw_sum_0) == FALSE), (kxw_sum_0 + c1), NA)
      
      wx = kxw/(kxw_sum )           # weights (Lee (2018) p. 4 equation (5))
      wx_0 = kxw_0/(kxw_sum_0 )  
    }
   
    # sum(wx[n_index])  # important condition == 1
    # loop for calculating local within transformation

    if(effects == "individual"){ 
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        # substract from y and x the local weighted sums of given t over all each i's 
        for(p in 1:P){  
          X_star_x[n_index, p] = (x[n_index] - x[i])^p - sum( (x[n_index_loo] - x[i])^p * wx[n_index_loo])
          X_star_0[n_index, p] = (x[n_index] - x[i])^p - sum( (x[n_index_loo] - x[i])^p * wx_0[n_index_loo])
        }  
        Y_star_x[n_index,1] = y[n_index] - sum(y[n_index_loo]*wx[n_index_loo])
        Y_star_0[n_index,1] = y[n_index] - sum(y[n_index_loo]*wx_0[n_index_loo])
      }
    }
    if(effects == "time"){
      for(t in 1:T){
        # t_index holds the indeces for of all i's for a given t
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        if(i %in% t_index){ # remove the index of the left-out observation
          t_index_loo = t_index[-which(t_index == i)]
        }else{
          t_index_loo = t_index 
        }
        # substract from y and x the local weighted sums of given t over all each i's
        for(p in 1:P){  
          X_star_x[t_index, p] = (x[t_index] - x[i])^p - sum( (x[t_index_loo] - x[i])^p * wx[t_index_loo])
          X_star_0[t_index, p] = (x[t_index] - x[i])^p - sum( (x[t_index_loo] - x[i])^p * wx_0[t_index_loo])
        }  
        Y_star_x[t_index,1] = y[t_index] - sum(y[t_index_loo]*wx[t_index_loo])
        Y_star_0[t_index,1] = y[t_index] - sum(y[t_index_loo]*wx_0[t_index_loo])
      }
    }
    if(effects == "indiv-time"){
      #individual
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        # substract from y and x the local weighted sums of given t over all each i's 
        for(p in 1:P){  
          X_star_x[n_index, p] = (x[n_index] - x[i])^p - sum( (x[n_index_loo] - x[i])^p * wx_a[n_index_loo])
          X_star_0[n_index, p] = (x[n_index] - x[i])^p - sum( (x[n_index_loo] - x[i])^p * wx_0_a[n_index_loo])
        }  
        Y_star_x[n_index,1] = y[n_index] - sum(y[n_index_loo]*wx_a[n_index_loo])
        Y_star_0[n_index,1] = y[n_index] - sum(y[n_index_loo]*wx_0_a[n_index_loo])
        
        # time sum(is.na(x_star))
        for(t in 1:T){
          # t_index holds the indeces for of all i's for a given t
          t_index = matrix(0, nrow = N, ncol = 1)
          for(j in 1:N){
            t_index[j,1] = (T*(j-1)+t)
          }
          if(i %in% t_index){ # remove the index of the left-out observation
            t_index_loo = t_index[-which(t_index == i)]
          }else{
            t_index_loo = t_index 
          }
          # substract from y and x the local weighted sums of given t over all each i's 
          for(p in 1:P){  
            X_star_x[t_index, p] = X_star_x[t_index, p] - sum( (x[t_index_loo] - x[i])^p * wx_b[t_index_loo])
            X_star_0[t_index, p] = X_star_0[t_index, p] - sum( (x[t_index_loo] - x[i])^p * wx_0_b[t_index_loo])
          }  
          Y_star_x[t_index,1] = Y_star_x[t_index,1] - sum(y[t_index_loo]*wx_b[t_index_loo])
          Y_star_0[t_index,1] = Y_star_0[t_index,1] - sum(y[t_index_loo]*wx_0_b[t_index_loo])
        }
        
        #indiv-time 
        for(p in 1:P){    
          X_star_x[, p] = X_star_x[, p] + sum( (x[-i]-x[i])^p * wx_c[-i])
          X_star_0[, p] = X_star_0[, p] + sum( (x[-i]-x[i])^p * wx_0_c[-i])
        }
        Y_star_x[,1] = Y_star_x[,1] + sum(y[-i]*wx_c[-i])
        Y_star_0[,1] = Y_star_0[,1] + sum(y[-i]*wx_0_c[-i])
      }
    }
    
    # Leave-one-out
    X_star_x = X_star_x[-i , ]
    Y_star_x = Y_star_x[-i , ]
    
    X_star_0 = X_star_0[-i , ]
    Y_star_0 = Y_star_0[-i , ]
    
    kx = kx[-i]
    kx_0 = kx_0[-i]
    
    beta_x = ginv(t(X_star_x)%*%diag(kx)%*%X_star_x + c2 )%*%t(X_star_x)%*%diag(kx)%*%Y_star_x
    beta_0 = ginv(t(X_star_0)%*%diag(kx_0)%*%X_star_0 + c2)%*%t(X_star_0)%*%diag(kx_0)%*%Y_star_0
    
    # Calculation of gamma_i(x) and gamma_i(0); p.5 with out equation number  
    
    if(effects =="indiv-time"){
      
      if(type == "ll"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_c
        
        gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_a 
        gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_b 
        gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_c
        
        gamma_it_x_a[i,1] = NA; gamma_it_x_b[i,1] = NA; gamma_it_x_c[i,1] = NA;
        gamma_it_0_a[i,1] = NA; gamma_it_0_b[i,1] = NA; gamma_it_0_c[i,1] = NA;
      }
      
      if(type == "lq"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_c
        
        gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_a 
        gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_b 
        gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_c
        
        gamma_it_x_a[i,1] = NA; gamma_it_x_b[i,1] = NA; gamma_it_x_c[i,1] = NA;
        gamma_it_0_a[i,1] = NA; gamma_it_0_b[i,1] = NA; gamma_it_0_c[i,1] = NA;
      }
      if(type == "l_cubic"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_c
        
        gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_a 
        gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_b 
        gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_c
        
        gamma_it_x_a[i,1] = NA; gamma_it_x_b[i,1] = NA; gamma_it_x_c[i,1] = NA;
        gamma_it_0_a[i,1] = NA; gamma_it_0_b[i,1] = NA; gamma_it_0_c[i,1] = NA;
      }
    }else{
      if(type == "ll"){
        gamma_it_x[,1] = (y - as.vector(beta_x)*(x - x[i]))*wx
        gamma_it_0[,1] = (y - as.vector(beta_0)*(x - 0))*wx_0
        gamma_it_x[i,1] = NA ; gamma_it_0[i,1] = NA
      }
      if(type == "lq"){
        gamma_it_x[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2  )*wx
        gamma_it_0[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0
        gamma_it_x[i,1] = NA ; gamma_it_0[i,1] = NA
      }
      if(type == "l_cubic"){
        gamma_it_x[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3 )*wx
        gamma_it_0[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3  )*wx_0
        gamma_it_x[i,1] = NA ; gamma_it_0[i,1] = NA
      }
    }
    
    if(effects == "individual"){
      # loop to sum of given i over all each t's
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        
        gamma_it_x[n_index,1] = sum(gamma_it_x[n_index_loo,1])
        gamma_it_0[n_index,1] = sum(gamma_it_0[n_index_loo,1])
      }
      n_single_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        n_single_index[j,1] = (T*(j-1)+1)
      }
      gamma_i_x = gamma_it_x[n_single_index,1]
      gamma_i_0 = gamma_it_0[n_single_index,1] # only keep first observation ; gamma_i_0 is of dimenstion (Nx1) 
    }    
    
    if(effects == "time"){
      # loop to sum of given t over all each i's
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        if(i %in% t_index){
          t_index_loo = t_index[-which(t_index == i)]
        }else{
          t_index_loo = t_index
        } 
        gamma_it_x[t_index,1] = sum(gamma_it_x[t_index_loo,1])
        gamma_it_0[t_index,1] = sum(gamma_it_0[t_index_loo,1])
      }
      
      gamma_t_x = as.matrix(na.omit(gamma_it_x[1:T,1]))
      gamma_t_0 = as.matrix(na.omit(gamma_it_0[1:T,1]))
    }
    
    if(effects == "indiv-time"){
      # individual 
      # loop to sum of given i over all each t's
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        if(i %in% n_index){ # remove the index of the left-out observation
          n_index_loo = n_index[-which(n_index == i)]
        }else{
          n_index_loo = n_index
        }
        
        gamma_it_x_a[n_index,1] = sum(gamma_it_x_a[n_index_loo,1])
        gamma_it_0_a[n_index,1] = sum(gamma_it_0_a[n_index_loo,1])
      }
      n_single_index = matrix(0, nrow = N, ncol = 1)
      
      for(j in 1:N){
        n_single_index[j,1] = (T*(j-1)+1)
      }
      gamma_i_x_a = gamma_it_x_a[n_single_index,1]
      gamma_i_0_a = gamma_it_0_a[n_single_index,1] # only keep first observation ; gamma_i_0 is of dimenstion (Nx1) 
      
      # time 
      # loop to sum of given t over all each i's
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        if(i %in% t_index){
          t_index_loo = t_index[-which(t_index == i)]
        }else{
          t_index_loo = t_index
        } 
        
        gamma_it_x_b[t_index,1] = sum(gamma_it_x_b[t_index_loo,1])
        gamma_it_0_b[t_index,1] = sum(gamma_it_0_b[t_index_loo,1])
      }
      
      gamma_t_x_b = as.matrix(na.omit(gamma_it_x_b[1:T,1]))
      gamma_t_0_b = as.matrix(na.omit(gamma_it_0_b[1:T,1]))
      
      gamma_x_c = sum(gamma_it_x_c[-i,1])
      gamma_0_c = sum(gamma_it_0_c[-i,1])    
      
    }
    
    if(effects == "individual"){
      mx[i,1] = mean(gamma_i_x - gamma_i_0) # mhat(x)
    }
    if(effects == "time"){
      mx[i,1] = mean(gamma_t_x - gamma_t_0) # mhat(x)
    }
    if(effects == "indiv-time"){
      mx[i,1] = (mean(gamma_i_x_a - gamma_i_0_a) + 
                   mean(gamma_t_x_b - gamma_t_0_b) - 
                   (gamma_x_c - gamma_0_c))
    }
    
    i = i+jump
    
  } # end while loop   
  
  CV_err = mean((y - mx[,1])^2)
  return(CV_err)
}


LMU.estimator = function( h, 
                          y, 
                          x,
                          N, 
                          T,
                          effects = c("individual", "time", "indiv-time"),
                          type = c("ll", "lq", "l_cubic"),
                          c1, c2){ 
  
  NT =  N*T 
  y = as.vector(y)
  x = as.vector(x)
  
  if(type == "ll"){
    Dx = matrix(0, nrow = NT, ncol = 1)  # Dx = 1st gradients vector 
    colnames(Dx) <- c("gradient.1")
    P = 1
    
  }
  
  if(type == "lq"){
    Dx = matrix(0, nrow = NT, ncol = 2)  # Dx = 2nd gradients vector
    colnames(Dx) <- c("gradient.1", "gradient.2")
    P = 2
  }
  if(type == "l_cubic"){
    Dx = matrix(0, nrow = NT, ncol = 3)  # Dx = 3rd gradients vector 
    colnames(Dx) <- c("gradient.1", "gradient.2", "gradient.3")
    P = 3
  }
  
  X_star_0 = matrix(0, nrow = length(x), ncol = P )
  X_star_x = matrix(0, nrow = length(x), ncol = P )
  
  Y_star_0 = matrix(0, nrow = length(x), ncol = 1 )
  Y_star_x = matrix(0, nrow = length(x), ncol = 1 )
  
  mx = matrix(0, nrow = NT, ncol = 1)  # mx = conditional mean vector
  colnames(mx) <- "mhat"
  
  
  if(effects == "individual"){
    kx_sum_0 = matrix(0, nrow = NT, ncol = 1)
    gamma_it_0 = matrix(0, nrow = NT, ncol = 1)
  }
  
  if(effects == "time"){
    kx_sum_0 = matrix(0, nrow = NT, ncol = 1)
    gamma_it_0 = matrix(0, nrow = NT, ncol = 1)
  }
  
  if(effects == "indiv-time"){
    kx_sum_0_a = matrix(0, nrow = NT, ncol = 1) 
    kx_sum_0_b = matrix(0, nrow = NT, ncol = 1)
    kx_sum_0_c = matrix(0, nrow = NT, ncol = 1)
    
    gamma_it_0_a = matrix(0, nrow = NT, ncol = 1)
    gamma_it_0_b = matrix(0, nrow = NT, ncol = 1)
    gamma_it_0_c = matrix(0, nrow = NT, ncol = 1)
    
  }
  
  
  kx_0 = dnorm(((x - 0)/h), 0, 1 )  # kernel weights
  # loop to sum of given i over all each t's  
  if(effects == "individual"){
    for(j in 1:N){
      n_index = ((1+T*(j-1)):(j*T))
      # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
      kx_sum_0[n_index,1] = sum(kx_0[n_index])  
    }
  }
  # loop to sum of given t over all each i's  
  if(effects == "time"){
    for(t in 1:T){
      # t_index holds the indeces for of all i's for a given t
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t) 
      }
      kx_sum_0[t_index, 1] = sum(kx_0[t_index]) # sum of t over all i's 
    }
  } 
  
  if(effects == "indiv-time"){
    # individual 
    for(j in 1:N){
      n_index = ((1+T*(j-1)):(j*T))
      kx_sum_0_a[n_index,1] = sum(kx_0[n_index])  # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
    }
    #time 
    for(t in 1:T){
      # t_index holds the indeces for of all i's for a given t
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t) 
      }
      kx_sum_0_b[t_index, 1] = sum(kx_0[t_index]) # sum of t over all i's 
    }
    
    # individual and time 
    kx_sum_0_c[, 1] = sum(kx_0)
    
  }
  
  if(effects == "indiv-time"){
    wx_0_a = kx_0 / (kx_sum_0_a + c1)  
    wx_0_b = kx_0 / (kx_sum_0_b + c1)  
    wx_0_c = kx_0 / (kx_sum_0_c + c1)  
  }else{
    wx_0 = kx_0 / (kx_sum_0 + c1)  #  weight matrix for x = 0 
  }
  
  # calculation of y* and x* for x = 0 
  
  if(effects == "individual"){
    
    # local within transformation (demeaning for all i's their individual means)
    
    for(j in 1:N){
      n_index = ((1+T*(j-1)):(j*T))
      for(p in 1:P){ 
        X_star_0[n_index,p] = (x[n_index] - 0)^p - sum( (x[n_index] - 0)^p * wx_0[n_index,1])  
      }
      Y_star_0[n_index,1] =  y[n_index]  - sum( y[n_index] * wx_0[n_index,1])
    }
  }
  if(effects == "time"){
    #local within transformation (demeaning for all t's their time means over all i's)
    
    for(t in 1:T){
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t) 
      }
      for(p in 1:P){  
        X_star_0[t_index,p] = (x[t_index] - 0)^p - sum( (x[t_index] - 0)^p * wx_0[t_index,1])
      }
      Y_star_0[t_index,1] = y[t_index] - sum( y[t_index] * wx_0[t_index,1])
    }
  }
  if(effects == "indiv-time"){
    # local within transformation (demeaning for all i's their individual means)
    
    for(j in 1:N){
      n_index = ((1+T*(j-1)):(j*T))
      for(p in 1:P){ 
        X_star_0[n_index,p] = (x[n_index] - 0)^p - sum( (x[n_index] - 0)^p * wx_0_a[n_index,1])
      }
      Y_star_0[n_index,1] = y[n_index] - sum( y[n_index] * wx_0_a[n_index,1])
    }
    
    #local within tranformation (demeaning for all t's their time means over all i's)
    for(t in 1:T){
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t) 
      }
      for(p in 1:P){ 
        X_star_0[t_index,p] = X_star_0[t_index,p] - sum( (x[t_index] - 0)^p * wx_0_b[t_index,1])
      }
      Y_star_0[t_index,1] = Y_star_0[t_index,1] - sum( y[t_index] * wx_0_b[t_index,1])
    }
    for(p in 1:P){ 
      X_star_0[,p] = X_star_0[,p] + sum( (x - 0)^p * wx_0_c)
    }
    Y_star_0[,1] = Y_star_0[,1] + sum( y * wx_0_c)
    
  }
  
  
  beta_0 = ginv(t(X_star_0)%*%diag(kx_0)%*%X_star_0 + c2)%*%t(X_star_0)%*%diag(kx_0)%*%Y_star_0 # regression with transformed data
  
  if(effects == "indiv-time"){
    
    if(type == "ll"){ 
      
      gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_a 
      gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_b 
      gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0_c
    }
    
    if(type == "lq"){ 
      
      gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_a 
      gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_b 
      gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 )*wx_0_c
    }
    
    
    if(type == "l_cubic"){ 
      
      gamma_it_0_a[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_a 
      gamma_it_0_b[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_b 
      gamma_it_0_c[,1] = (y - as.vector(beta_0[1])*(x - 0) - as.vector(beta_0[2])*(x - 0)^2 - as.vector(beta_0[3])*(x - 0)^3)*wx_0_c
    }
    
  }else{
    
    if(type == "ll"){ gamma_it_0[,1] = (y - as.vector(beta_0[1])*(x - 0))*wx_0   }
    if(type == "lq"){ gamma_it_0[,1] = (y - as.vector(beta_0[1])*(x - 0) -  as.vector(beta_0[2])*(x - 0)^2 )*wx_0  }
    if(type == "l_cubic"){  gamma_it_0[,1] = (y - as.vector(beta_0[1])*(x - 0) -  as.vector(beta_0[2])*(x - 0)^2 -  as.vector(beta_0[3])*(x - 0)^3 )*wx_0 }
    
  }
  if(effects == "individual"){
    for(j in 1:N){
      n_index =  ((1+T*(j-1)):(j*T))
      gamma_it_0[n_index, 1] = sum(gamma_it_0[n_index,1])  # take sums over t for each individual. 
    }                                                      # for each t of an indiv we have the same sum
    n_single_index = matrix(0, nrow = N, ncol = 1)
    for(j in 1:N){
      n_single_index[j,1] = (T*(j-1)+1)
    }
    gamma_i_0 = as.matrix(gamma_it_0[n_single_index,1]) # only keep first observation ; gamma_i_0 is of dimenstion (Nx1)
   
  }
  
  if(effects == "time"){
    for(t in 1:T){
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t)
      }
      gamma_it_0[t_index, 1] = sum(gamma_it_0[t_index,1])  
    }
    gamma_t_0 = as.matrix(gamma_it_0[1:T,1]) 

  }
  
  if(effects == "indiv-time"){
    # individual : gamma_a 
    for(j in 1:N){
      n_index =  ((1+T*(j-1)):(j*T))
      gamma_it_0_a[n_index, 1] = sum(gamma_it_0_a[n_index,1])  
    }                                                   
    n_single_index = matrix(0, nrow = N, ncol = 1)
    for(j in 1:N){
      n_single_index[j,1] = (T*(j-1)+1)
    }
    gamma_i_0_a = as.matrix(gamma_it_0_a[n_single_index,1]) 
    
    # time : gamma_b
    for(t in 1:T){
      t_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        t_index[j,1] = (T*(j-1)+t)
      }
      gamma_it_0_b[t_index, 1] = sum(gamma_it_0_b[t_index,1])   
    }
    gamma_t_0_b = as.matrix(gamma_it_0_b[1:T,1]) 
  
    # indiv and time: gamma_b 
    gamma_0_c = sum(gamma_it_0_c[,1]) 
  }
  
  
 
  # Step 2: calculate gamma_i(x), Lee et al. (2019) p. 5
  
  # loop for local linear regression for each i
  for(i in 1:NT){  
    if(effects == "individual"){
      kx_sum = matrix(0, nrow = NT, ncol = 1)
      gamma_it_x = matrix(0, nrow = NT, ncol = 1)
    }
    
    if(effects == "time"){
      kx_sum = matrix(0, nrow = NT, ncol = 1)
      gamma_it_x = matrix(0, nrow = NT, ncol = 1)
    }
    
    if(effects == "indiv-time"){
      kx_sum_a = matrix(0, nrow = NT, ncol = 1) 
      kx_sum_b = matrix(0, nrow = NT, ncol = 1)
      kx_sum_c = matrix(0, nrow = NT, ncol = 1)
      
      gamma_it_x_a = matrix(0, nrow = NT, ncol = 1)
      gamma_it_x_b = matrix(0, nrow = NT, ncol = 1)
      gamma_it_x_c = matrix(0, nrow = NT, ncol = 1)
    }
    
    
    # Calculation of weighting matrix and transformed data to remove the fixed effect 
    # Lee et al. (2019) p.4 equation (5)
    
    kx = dnorm(((x - x[i])/h), 0, 1 )
    
    if(effects == "individual"){
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))          
        kx_sum[n_index,1] = sum(kx[n_index])  # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
      }
    }
    
    if(effects == "time"){
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t) 
        }
        kx_sum[t_index, 1] = sum(kx[t_index])  
      }
    }
    
    if(effects == "indiv-time"){
      # individual 
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        kx_sum_a[n_index,1] = sum(kx[n_index])  # Note: by (1+T*(j-1)):(j*T) we calculate for each indiv i its sum over all t. 
      }
      #time 
      for(t in 1:T){
        # t_index holds the indeces for of all i's for a given t
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t) 
        }
        kx_sum_b[t_index, 1] = sum(kx[t_index]) # sum of t over all i's 
      }
      
      # individual and time 
      kx_sum_c[, 1] = sum(kx)
      
    }
    
    # c1 =  1/T^2 
    if(effects == "indiv-time"){    # The constant leads to a strong impact if the kx values are very small
      wx_a = kx / (kx_sum_a + c1)   
      wx_b = kx / (kx_sum_b + c1)  
      wx_c = kx / (kx_sum_c + c1)
    }else{
      wx = kx / (kx_sum + c1)  
    }                          
    if(effects == "individual"){
      
      # local within tranformation (demeaning for all i's their individual means)
      
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        for(p in 1:P){ 
          X_star_x[n_index,p] = (x[n_index] - x[i])^p - sum( (x[n_index] - x[i])^p * wx[n_index,1])  
        }
        Y_star_x[n_index,1] =  y[n_index]  - sum( y[n_index] * wx[n_index,1])
      }
    }
    if(effects == "time"){
      #local within tranformation (demeaning for all t's their time means over all i's)
      
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t) 
        }
        for(p in 1:P){  
          X_star_x[t_index,p] = (x[t_index] - x[i])^p - sum( (x[t_index] - x[i])^p * wx[t_index,1])
        }
        Y_star_x[t_index,1] = y[t_index]  - sum( y[t_index] * wx[t_index,1])
      }
    }
    if(effects == "indiv-time"){
      # local within tranformation (demeaning for all i's their individual means)
      
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        for(p in 1:P){ 
          X_star_x[n_index,p] = (x[n_index] - x[i])^p - sum( (x[n_index] - x[i])^p * wx_a[n_index,1])
        }
        Y_star_x[n_index,1] = y[n_index] - sum( y[n_index] * wx_a[n_index,1])
      }
      
      #local within tranformation (demeaning for all t's their time means over all i's)
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t) 
        }
        for(p in 1:P){ 
          X_star_x[t_index,p] = X_star_x[t_index,p] - sum( (x[t_index] - x[i])^p * wx_b[t_index,1])
        }
        Y_star_x[t_index,1] = Y_star_x[t_index,1] - sum( y[t_index] * wx_b[t_index,1])
      }
      for(p in 1:P){ 
        X_star_x[,p] = X_star_x[,p] + sum( (x - x[i])^p * wx_c)
      }
      Y_star_x[,1] = Y_star_x[,1] + sum( y * wx_c)
      
    }
    
    beta_x = ginv(t(X_star_x)%*%diag(kx)%*%X_star_x + c2)%*%t(X_star_x)%*%diag(kx)%*%Y_star_x # regression with transformed data
    
    if(type == "ll"){
      Dx[i,1] = beta_x[1]
    }
    if(type == "lq"){
      Dx[i,1] = beta_x[1]
      Dx[i,2] = beta_x[2]
    }
    if(type == "l_cubic"){
      Dx[i,1] = beta_x[1]
      Dx[i,2] = beta_x[2]
      Dx[i,3] = beta_x[3]
    }
    if(effects=="indiv-time"){
      if(type == "ll"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx_c
      }
      if(type == "lq"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 )*wx_c
      }
      if(type == "l_cubic"){ 
        gamma_it_x_a[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_a 
        gamma_it_x_b[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_b 
        gamma_it_x_c[,1] = (y - as.vector(beta_x[1])*(x - x[i]) - as.vector(beta_x[2])*(x - x[i])^2 - as.vector(beta_x[3])*(x - x[i])^3)*wx_c
      }
    }else{
      if(type == "ll"){ gamma_it_x[,1] = (y - as.vector(beta_x[1])*(x - x[i]))*wx   }
      if(type == "lq"){ gamma_it_x[,1] = (y - as.vector(beta_x[1])*(x - x[i]) -  as.vector(beta_x[2])*(x - x[i])^2 )*wx  }
      if(type == "l_cubic"){  gamma_it_x[,1] = (y - as.vector(beta_x[1])*(x - x[i]) -  as.vector(beta_x[2])*(x - x[i])^2 -  as.vector(beta_x[3])*(x - x[i])^3 )*wx}
    }
    
    
    if(effects == "individual"){
      for(j in 1:N){
        n_index = ((1+T*(j-1)):(j*T))
        gamma_it_x[n_index, 1] = sum(gamma_it_x[n_index,1])  
      }
      n_single_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        n_single_index[j,1] = (T*(j-1)+1) 
      }
      gamma_i_x = as.matrix(gamma_it_x[n_single_index,1])
    }
    
    if(effects == "time"){
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        gamma_it_x[t_index, 1] = sum(gamma_it_x[t_index,1])   
      }
      gamma_t_x = as.matrix(gamma_it_x[1:T,1]) 
    }
    
    
    if(effects == "indiv-time"){
      # individual : gamma_a 
      for(j in 1:N){
        n_index =  ((1+T*(j-1)):(j*T))
        gamma_it_x_a[n_index, 1] = sum(gamma_it_x_a[n_index,1])  
      }                                                     
      n_single_index = matrix(0, nrow = N, ncol = 1)
      for(j in 1:N){
        n_single_index[j,1] = (T*(j-1)+1)
      }
      gamma_i_x_a = as.matrix(gamma_it_x_a[n_single_index,1]) 
      # time : gamma_b
      for(t in 1:T){
        t_index = matrix(0, nrow = N, ncol = 1)
        for(j in 1:N){
          t_index[j,1] = (T*(j-1)+t)
        }
        gamma_it_x_b[t_index, 1] = sum(gamma_it_x_b[t_index,1])   
      }
      gamma_t_x_b = as.matrix(gamma_it_x_b[1:T,1]) 
      # indiv qnd time: gamma_b 
      gamma_x_c = sum(gamma_it_x_c[,1]) # we only need to take the values of the first id 1:T for gamma_i (dimension (Tx1))
    }
    
    # print(i)
    
    if(effects == "individual"){
      mx[i,] = as.numeric(mean(gamma_i_x - gamma_i_0))  
    }
    
    if(effects == "time"){
      mx[i,] = as.numeric(mean(gamma_t_x - gamma_t_0))  
    }
    
    if(effects == "indiv-time"){
      mx[i,] = ( as.numeric(mean(gamma_i_x_a - gamma_i_0_a)) + 
                   as.numeric(mean(gamma_t_x_b - gamma_t_0_b)) -
                   as.numeric((gamma_x_c - gamma_0_c)) )
    }
  }
  estimates = cbind(mx, Dx)
  return( estimates )
}


npregivfe  = function(y,
                   x,
                   z, 
                   c,
                   tol,
                   N,
                   T, 
                   max.iter,
                   bw_method = c("optimal","plug_in", "sylverman"),
                   effects =c("individual", "time", "indiv-time"),
                   type = c("ll", "lq", "l_cubic"), 
                   bws 
){
  
  y = as.matrix(as.numeric(y))
  x = as.matrix(as.numeric(x))
  z = as.matrix(as.numeric(z)) 
  
  output = list()
  crit_k = list()
  phi_k_list = list()
  max.iter = max.iter
  
  NT=N*T 
  
  c1 =  1/T^2 ; c2 = 1/((N*T)^2)
  
  if(bw_method == "optimal"){
    print(paste("compute bandwidth E(Y|z)"))
    Ehat_yz_bw = optimize(f = LMU.CVCMB,
                          interval = c(0,5*sd(z)),
                          lower = 0, upper =  5*sd(z),
                          tol = 0.001, 
                          x = z, y = y ,N = N, T = T, 
                          effects = effects , 
                          type =  type, 
                          jump = 1, 
                          c1 = c1, c2 = c2)$minimum
  }
  if(bw_method == "plug_in"){ Ehat_yz_bw = bws[1] }
  if(bw_method == "sylverman"){Ehat_yz_bw = sqrt(var(z))*NT^(-1/7)}
  
  print(paste("Estimate E(Y|z)"))
  Ehat_yz = LMU.estimator(h=Ehat_yz_bw,
                           x = z, y = y, N = N, T = T,
                           effects = effects,
                           type = type, 
                           c1 = c1, c2 = c2)[,"mhat"]
  
  k=1
  while( k  < max.iter){
    if(k == 1){ 
      # 1) initial guess for phi(x) = ghat = E(Y|x)
      print(paste("Iteration", k, "compute bandwidth E(Y|x) "))
      #  phi_k1_bw = 1.002509
      if(bw_method == "optimal"){
        phi_k1_bw = optimize(f = LMU.CVCMB,
                             interval = c(0,5*sd(x)),
                             lower = 0, upper = 5*sd(x),
                             tol = 0.001,
                             x = x, y = y ,N = N,T = T,
                             effects = effects,
                             type = type,
                             jump = 1,
                             c1 = c1, c2 = c2)$minimum
        
      }
      if(bw_method == "plug_in"){ phi_k1_bw = bws[2] }
      if(bw_method == "sylverman"){phi_k1_bw = sqrt(var(x))*NT^(-1/7)}
      
      print(paste("Iteration", k, "Estimate E(Y|x) "))
      
      phi_k1 = LMU.estimator(h=phi_k1_bw,
                              x = x, y = y, N = N, T = T,
                              effects = effects,
                              type = type,
                              c1 = c1, c2 = c2)[,"mhat"]
    }else{ 
      phi_k1 = phi_k 
    }
    # 2) E(Y - ghat_k(x)|z)  
    y2_k = y - phi_k1  
    
    if(k == 1){ 
      print(paste("Iteration", k, "compute bandwidth E(Y - ghat_k(x)|z)"))
      
      # Ehat_y2z_bw = 0.1419789
      if(bw_method == "optimal"){
        Ehat_y2z_bw = optimize(f = LMU.CVCMB,
                               interval = c(0,5*sd(z)),
                               lower = 0, upper = 5*sd(z),
                               tol = 0.001,
                               x = z, y = y2_k ,N = N, T = T,
                               effects = effects,
                               type = type, 
                               jump = 1, 
                               c1 = c1, c2 = c2)$minimum
      }
      if(bw_method == "plug_in"){ Ehat_y2z_bw = bws[3] }
      if(bw_method == "sylverman"){Ehat_y2z_bw = sqrt(var(z))*NT^(-1/7)}
    }
    print(paste("Iteration", k, "Estimate E(Y - ghat_k(x)|z) "))
    Ehat_y2z = LMU.estimator(h=Ehat_y2z_bw, x=z, y=y2_k, N=N, T=T, effects = effects, type = type, c1 = c1, c2 = c2)[,"mhat"]
  
    # 3) E(E(Y - ghat_k(x)|z) | x)
    if(k == 1){ 
      print(paste("Iteration", k, "compute bandwidth E(E(Y - ghat_k(x)|z) | x)"))
      
      # EEhat_y2zx_bw  = 1.097986
      if(bw_method == "optimal"){
        h2 = sqrt(var(x))*NT^(-1/7) 
       
        EEhat_y2zx_bw = optimize(f = LMU.CVCMB,
                                 interval = c(0,5*sd(x)),
                                 lower = 0, upper = 5*sd(x),
                                 tol = 0.001,
                                 x=x, y=Ehat_y2z ,N=N,T=T,
                                 effects = effects,
                                 type = type,
                                 jump = 1,
                                 c1 = c1, c2 = c2)$minimum
      }
      if(bw_method == "plug_in"){ EEhat_y2zx_bw = bws[4] }
      if(bw_method == "sylverman"){EEhat_y2zx_bw = sqrt(var(x))*NT^(-1/7)}
    }
    print(paste("Iteration", k, "Estimate E(E(Y - ghat_k(x)|z) | x)"))
    
    EEhat_y2zx = LMU.estimator( h = EEhat_y2zx_bw,
                                 x = x, y = Ehat_y2z, N = N, T = T,
                                 effects = effects,
                                 type = type,
                                 c1 = c1, c2 = c2)[,"mhat"]
    
    # 4) Calculate phi_k
    
    phi_k =  phi_k1 + c*EEhat_y2zx
    
    phi_k_list[[k]] <- phi_k 
    
    # 5) Stropping criteria 
    if(k == 1){ 
      
      print(paste("Iteration", k, "compute bandwidth E(phi_hat_k(x) | z)"))
      # Ehat_phi_kz_bw = 0.17984
      
      if(bw_method == "optimal"){
        Ehat_phi_kz_bw = optimize(f = LMU.CVCMB,
                                  interval = c(0,5*sd(z)),
                                  lower = 0, upper = 5*sd(z),
                                  tol = 0.001,
                                  x=z, y=phi_k ,N=N,T=T,
                                  effects = effects,
                                  type = type,
                                  jump = 1,
                                  c1 = c1, c2 = c2)$minimum
      } 
      if(bw_method == "plug_in"){  Ehat_phi_kz_bw = bws[5] }
      if(bw_method == "sylverman"){Ehat_phi_kz_bw = sqrt(var(z))*NT^(-1/7)}
    }
    print(paste("Iteration", k, "Estimate E(phi_hat_k(x) | z)"))
    
    Ehat_phi_kz = LMU.estimator( h = Ehat_phi_kz_bw,
                                  x = z, y = phi_k,
                                  N = N, T = T,
                                  effects = effects,
                                  type = type, c1 = c1, c2 = c2)[,"mhat"]
    
    crit_k[[k]] = sum( ( ( Ehat_yz - Ehat_phi_kz ) / Ehat_yz )^2 )
    if( k == 1){
      k = k + 1
    }else{
      if(  abs( 100 - ( crit_k[[k]]*100 / crit_k[[k-1]] ) ) <= tol  ){
        break
      }else{
        k = k + 1}
    }
    print(k)
  }
  
  bws = c(Ehat_yz_bw, phi_k1_bw, Ehat_y2z_bw, EEhat_y2zx_bw, Ehat_phi_kz_bw)
  
  phi_k_df <- do.call(cbind, phi_k_list) 
  
  cols = c()  
  
  for(i in 1:ncol(phi_k_df)){ cols[i] <- as.character( paste( "k",i,sep = "" ) )}  
  colnames(phi_k_df) <- cols 
  
  crit_k_vec = do.call( rbind, crit_k )
  colnames( crit_k_vec ) <- "crit"
  
  output[[1]] <- phi_k_df
  output[[2]] <- crit_k_vec
  output[[3]] <- bws
  
  return(output)
}




npregiv2 = function(y,z,w, c,tol, max.iter, bw_method = c("optimal","plug_in", "sylverman"), bws ){
  
  
  
  y = as.vector(y)
  z = as.vector(z)
  w = as.vector(w)
  
  crit_k = list()
  phi_k_list = list()
  output = list()
  max.iter = max.iter
  
  NT = length(y)
  
  # Part of stopping criterion (only once to be computed)
  
  if(bw_method == "optimal"){ Ehat_yw_bw = npregbw(y~w, regtype = "ll", bwmethod = "cv.ls")$bw   }
  if(bw_method == "plug_in"){ Ehat_yw_bw = bws[1] }
  if(bw_method == "sylverman"){ Ehat_yw_bw = sqrt(var(w))*NT^(-1/7) }
  
  
  print(c("Ehat_yw_bw:", Ehat_yw_bw )) 
  
  Ehat_yw = fitted(npreg(bws = Ehat_yw_bw, tydat = y, txdat = w, regtype = "ll"))
  
  c = c   # constant to be specified: c<1 is requiered
  tol = tol # "stabilizing" tolerance threshold: when stopping criterion "stabilizes" i.e. changes in the criterion are less than "tol"
  k=1
  while ( k  < max.iter){
    if(k == 1){ 
      # 1) initial guess for phi(z) = ghat = E(Y|z)
      if(bw_method=="optimal"){phi_k1_bw = npregbw(y~z, regtype = "ll",  bwmethod = "cv.ls")$bw }
      if(bw_method=="plug_in"){phi_k1_bw = bws[2] }
      if(bw_method=="sylverman"){phi_k1_bw = sqrt(var(z))*NT^(-1/7) }
      
      print(c("phi_k1_bw:", phi_k1_bw)) 
      phi_k1 = fitted(npreg(bws = phi_k1_bw, tydat =  y, txdat = z, regtype = "ll" )) 
    }else{ 
      # if k >1, use phi_k-1(z)
      phi_k1 = phi_k 
    }
    # 2) E(Y - ghat_k(z)|w)  
    y2_k = y - phi_k1  # transform the data
    if(k==1){
      if(bw_method=="optimal"){ Ehat_y2w_bw = npregbw(y2_k~w, regtype = "ll", bwmethod = "cv.ls" )$bw }
      if(bw_method=="plug_in"){ Ehat_y2w_bw = bws[3] }
      if(bw_method=="sylverman"){ Ehat_y2w_bw = sqrt(var(w))*NT^(-1/7) }  
    }
    
    print(c("Ehat_y2w_bw:", Ehat_y2w_bw )) 
    Ehat_y2w = fitted(npreg(bws = Ehat_y2w_bw, tydat = y2_k, txdat = w, regtype = "ll" ))  
    
    
    
    # 3) E(E(Y - ghat_k(z)|w) | z)
    if(k==1){
      if(bw_method=="optimal"){ EEhat_y2wz_bw = npregbw(Ehat_y2w~z, regtype = "ll", bwmethod = "cv.ls" )$bw }
      if(bw_method=="plug_in"){ EEhat_y2wz_bw = bws[4]}
      if(bw_method=="sylverman"){EEhat_y2wz_bw = sqrt(var(z))*NT^(-1/7) }  
    }
    print(c("EEhat_y2wz_bw:", EEhat_y2wz_bw )) 
    EEhat_y2wz = fitted(npreg(bws = EEhat_y2wz_bw, tydat =  Ehat_y2w, txdat = z, regtype = "ll" ) )
    
    
    # 4) Calculate phi_k
    
    phi_k =  phi_k1 + c*EEhat_y2wz
    
    phi_k_list[[k]] <- phi_k 
    
    # 5) Stropping criteria 
    # || ( E(Y|w) - E(phi_k(z)|w) ) / E(Y|w) ||?  
  
    # E(phi_hat_k(z) | w)
    if(k==1){
      if(bw_method=="optimal"){Ehat_phi_kw_bw = npregbw(phi_k~w, regtype = "ll", bwmethod = "cv.ls" )$bw }
      if(bw_method=="plug_in"){ Ehat_phi_kw_bw = bws[5]}
      if(bw_method=="sylverman"){Ehat_phi_kw_bw = sqrt(var(w))*NT^(-1/7) }
      
    }  
    print(c("Ehat_phi_kw_bw:", Ehat_phi_kw_bw))
    Ehat_phi_kw = fitted(npreg(bws = Ehat_phi_kw_bw, tydat = phi_k, txdat = w, regtype = "ll" ))
    
    crit_k[[k]] = sum(((Ehat_yw - Ehat_phi_kw )/Ehat_yw)^2)
    if( k == 1){
      k = k + 1
    }else{
      if(    
        abs(100 - (crit_k[[k]]*100/crit_k[[k-1]])) <= tol # here tol change in % 
      ){
        # return(phi_k)
        break
      }else{
        k = k + 1}
    }
    print(k)
  }
  
  bws = c(Ehat_yw_bw, phi_k1_bw, Ehat_y2w_bw, EEhat_y2wz_bw, Ehat_phi_kw_bw)
  phi_k_df <- do.call("cbind", phi_k_list) 
  cols = c()  
  for(i in 1:ncol(phi_k_df)){ cols[i] <- as.character(paste("k",i,sep = ""))}  
  colnames(phi_k_df) <- cols 
  
  crit_k_vec = do.call("rbind", crit_k)
  colnames(crit_k_vec) <- "crit"
  
  output[[1]] <- phi_k_df
  output[[2]] <- crit_k_vec
  output[[3]] <- bws
  
  return(output)
  
}

# Data Generating process 

dgpivfe = function(N,T){ 
  
  library(data.table)
  NT=N*T
  iota_T <- rep(1,T)
  iota_NT <- rep(1,NT)
  iota_N <- rep(1,N)
  trend_T <- c(1:T)
  ind_N <- c(1:N)
  col_T <- iota_N %x% trend_T 
  col_N <- ind_N %x% iota_T # numerate the industries from 1 to N
  t = c(1:T)
  
  # Parts of DGP in Racine (2019, p. 288) 
  
  # Correlation parameters 
  rho.xz = 0.2  
  rho.ux = 0.8
  sigma.u = 0.05
  
  # DGP of the nonlinear model 
  dgp = function(x){x^2}
  
  # Generate endogeneity in x 
  
  v1 = rnorm(NT,0,0.8)
  v2 = rnorm(NT,0,0.15)
  eps = rnorm(NT,0,sigma.u)
  u = rho.ux*v1 + eps # E(u,x) != 0 

  # Individual effects
  alpha_i = runif(N, 0,1.5)
  alpha_i = as.matrix(alpha_i %x% iota_T)
  
  # Time fixed effects 
  delta_t = runif(T, 0,1.7)
  delta_t = as.matrix(rep(delta_t,N))
  
  # Regressor correlated with individual effects
  x = matrix(NT, nrow = NT, ncol = 1)
  
    for(i in 1:NT){
      x[i,1] = rnorm(1, mean = alpha_i[i]+delta_t[i], sd = 1 ) 
    }
  z = rho.xz*x + v2
  x = x + v1
  
  
  y =  as.matrix(dgp(x=x)) + alpha_i + delta_t + u 
  
  
  Data = cbind(col_N, col_T, y, x, z)
  colnames(Data) <- c("col_N", "col_T", "y", "x", "z")
  Data = data.table(Data)
  
  return(Data)
  
}
