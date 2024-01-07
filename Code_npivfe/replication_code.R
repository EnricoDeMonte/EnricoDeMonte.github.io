# Replication code 
# De Monte (2023). Nonparametric Instrumental Regression with Two-Way Fixed Effects.
# Journal of Econoteric Methods


library(latex2exp)
library(colorRamps)
library(data.table)
library(np)
library(xtable)
library(np)

source("C:/npivfe/npivfe.R")

# DGP ----
set.seed(49)

# Panel structure 
N = 100  
T= 20   
NT=N*T
iota_T <- rep(1,T)
iota_N <- rep(1,N)
trend_T <- c(1:T)
ind_N <- c(1:N)
col_T <- iota_N %x% trend_T 
col_N <- ind_N %x% iota_T 


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
xi_i = runif(N, 0,1.5)
xi_i = as.matrix(xi_i %x% iota_T)

# Time fixed effects 
delta_t = runif(T, 0,1.7)
delta_t = as.matrix(rep(delta_t,N))

# Regressor correlated with individual effects
x = matrix(NT, nrow = NT, ncol = 1)

for(i in 1:NT){
  x[i,1] = rnorm(1, mean = xi_i[i]+delta_t[i], sd = 1 ) 
}

z = rho.xz*x + v2
x = x + v1

# Dependent variable 
y = as.matrix(dgp(x=x)) + xi_i + delta_t + u 

covmat = cov(cbind(x,z,xi_i, delta_t, u) ) 
colnames(covmat) <- c("x","z","xi_i", "delta_t", "u")
rownames(covmat) <- c("x","z","xi_i", "delta_t", "u")
covmat[lower.tri(covmat, diag = TRUE)] = 0 
covmat = round(covmat,2)
print(xtable(covmat, caption = "Covariance matrix"  ), include.rownames = TRUE)


Data = data.table(cbind(col_N, col_T, y, x, z))
names(Data) <- c("col_N", "col_T","y","x", "z")



# NP IV with two-way fixed effects ----
# Landweber-Fridman with LMU estimator 

rStart <- Sys.time()
res.npivfe = npregivfe( y = y,
                        x = x,
                        z = z, 
                        c = 0.5,   # constant to be specified: c<1 is requiered
                        tol = 1,   # "stabilizing" tolerance threshold: when stopping criterion "stabilizes" i.e. changes in the criterion are less than "tol"
                        N = N,
                        T = T, 
                        max.iter = 100,
                        bw_method = "optimal",
                        effects = "indiv-time",
                        type = "ll" )
rEnd <- Sys.time()
print(paste0(round(rEnd - rStart, 2), attr(rEnd - rStart, "units")))

# Investigate regularized solution path 
phi_k_mat = res.npivfe[[1]]

write.csv(phi_k_mat, "C:/npivfe/phi_k_mat.csv")

writecsv(res.npivfe[[1]], )
# Total number of iterations
k_bar = ncol(phi_k_mat) 

# Fitted values of the function of interest
fit.npivfe = phi_k_mat[, k_bar]  

# Plot the initial guess phi_0
plot(x[order(x)], phi_k_mat[, 1][order(x)],
     col = blue2green(k_bar)[1],
     type = "l",
     lwd = 2,
     ylim = range(y),
     ylab = "Y", xlab = "X", main = "" ) 
# Plot the successive estimates phi_k
for(k in 1:k_bar){
  lines( x[order(x)], phi_k_mat[, k][order(x)], 
         type = "l",
         lwd = 2, 
         col =  blue2green(k_bar)[k] )
}
# Add the true curve  
curve(dgp, 
      min(x), max(x),
      add = TRUE, 
      col = "red", 
      lwd = 2, lty = 2)
legend(-2, 40, 
       legend = c("DGP", TeX( '$\\hat{\\varphi}(x)_{0}$'), 
                  TeX('$\\hat{\\varphi}(x)_{\\bar{k}}$')), 
       col = c("red",blue2green(k_bar)[1], blue2green(k_bar)[k_bar]),
       lty = c(2, 1, 1),
       lwd = c(2, 2, 2),
       cex = 0.65)


# Stopping criterion ----- 

val_crit_k = res.npivfe[[2]]
write.csv(val_crit_k, "C:/npivfe/val_crit_k.csv")
plot(1:k_bar,
     val_crit_k, type = "l", 
     col = "black", lwd = 2, 
     xlim = c(0,50),
     main = "", 
     xlab = "k-iterations", ylab= TeX('$Criterion_k $'))

points(k_bar, val_crit_k[k_bar] )

# Comparison with other nonparamteric estimators ----

# Local linear kernel regression (LL)
fit.np_bw = npregbw(xdat = as.vector(x), ydat = as.vector(y), 
                    ckertype="gaussian", 
                    regtype = "ll", 
                    bwmethod = "cv.ls")
# Obtain fitted values 
fit.np = fitted(npreg(fit.np_bw))

# LMU fixed effects estimator 
# Obtain optimal banwdth
h_LMU = optim(par = sqrt(var(x))*NT^(-1/7),
              LMU.CVCMB, 
              method = "Nelder-Mead",
              control = list(reltol=0.001), 
              x = x,
              y = y,
              N = N,
              T = T,
              effects = "indiv-time",
              type = "ll",
              jump = 1, 
              c1 =  1/T^2, c2 = 1/((N*T)^2))$par

# Apply the LMU estimator 
fit.LMU = LMU_estimator(h = h_LMU,
                        x = x,
                        y = y,
                        N = N,
                        T = T,
                        effects = "indiv-time",
                        type = "ll",
                        c1 =  1/T^2, c2 = 1/((N*T)^2))[,"mhat"]



# Nonparametric IV (without fixed effects)
# fit.npiv = npregiv(y = y, 
#                    z = x, 
#                    w = z, 
#                    method = "Landweber-Fridman")$phi.mat

# Nonparametric IV with Landweber-Fridman 
res.npivf2 = npregiv2( y = y,
                        z = x, 
                        w = z, 
                        c = 0.5,
                        tol=1,
                        max.iter = 100,
                        bw_method = "optimal" )

fit.npivf2 = res.npivf2[[1]][,ncol(res.npivf2[[1]])]

res.estimators = cbind(col_N, col_T, fit.np, fit.LMU, fit.npivf2, fit.npivfe )
write.csv(res.estimators, "C:/npivfe/res.estimators.csv")

# Plot the results 
plot(x, y, 
     pch = 21, 
     col = "grey", 
     ylab = "Y", 
     xlab = "X" , main = "")
lines(x[order(x)], fit.np[order(x)],
      col= "orange",
      type = "l",
      lty =1,
      lwd = 2)
lines(x[order(x)], fit.npivf2[order(x)],
      col= "black",
      type = "l",
      lty = 1,
      lwd = 2)
lines(x[order(x)], fit.LMU[order(x)],
      col= "darkblue",
      type = "l",
      lty = 1,
      lwd = 2)
lines(x[order(x)],
      fit.npivfe[order(x)],
      type = "l",
      col = blue2green(k_bar)[k_bar], lty =1  ,lwd = 2)
curve(dgp(x),
      col="red", 
      add=TRUE, 
      lty = 2, 
      lwd = 2)
legend(-1,40, 
       legend = c("DGP", "local-linear", "L-F", "LMU", "L-F/LMU"), 
       col = c("red", "orange", "black" , "darkblue", blue2green(k_bar)[k_bar]),
       lty = c(2,1,1,1,1), 
       lwd = c(2,2,2,2,2), 
       cex = 0.65 )


# Bootstrap ----


# Use bandwidths from initial estimation of the conditional mean. 
bws = res.npivfe[[3]]

# Follow Malikov et al. (2020)  
# Bootstrap iterations 
B = 400

# Use bandwidths from initial estimation of the conditional mean. 
bws = res.npivfe[[3]]


# Compute residuals u_hat 
v_hat = y - fit.npivfe 

# Center residuals u_hat
v_hat_c = as.numeric(scale(v_hat, center = TRUE))

# Create empty lists to store bootstrapp estimates
phi_hat_boot_list <- list()

for(b in 1:B){
  
  # Generate bootstrap weights 
  # each individual i keeps its weight for all t
  b_i = sample( c((1+sqrt(5))/2, (1-sqrt(5))/2 ), 
                n, 
                prob = c( ( sqrt( 5 ) - 1 ) / ( 2 * sqrt( 5 ) ), ( sqrt( 5 )+1 ) / ( 2*sqrt( 5 ) ) ),
                replace = TRUE )
  
  # Generate new disturbances
  v_hat_b = (b_i %x% iota_T) * v_hat_c
  
  # Generate new outcome variable
  y_b = phi_hat + v_hat_b
  
  # Re-estimate phi_hat using y_b 
  res.npivfe_b = npivfe(y = y_b, 
                        x = x, 
                        z = z, 
                        N = N, 
                        T = T,
                        bw_method = "plug_in", 
                        effects = "indiv-time", 
                        type = "ll", 
                        c = 0.5,
                        tol = 1, 
                        bws = bws,
                        max.iter = 100)
  
  phi_hat_b = res.npivfe_b[[1]][,length(res.npivfe_b[[2]])]
  
  phi_hat_boot_list[[b]] <- phi_hat_b
}

# Bind the results in a (NT x B) matrix.  
phi_boot = do.call(cbind, phi_hat_boot_list)

# Compute the pointwise 95% confidence intervals 
phi_025 = apply(phi_boot, 1, quantile, probs=0.025)
phi_975 = apply(phi_boot, 1, quantile, probs=0.975)

# Plot the estimates along with the confidence intervals 
plot(x[order(x)], fit.npivfe[order(x)],
     ylim = c(0,40),
     type = "l",
     lty = 1,
     lwd = 2,
     ylab = "Y",
     xlab = "X",
     main = "")
lines(x[order(x)], phi_025[order(x)],
      lty = 2, 
      col = "black")
lines(x[order(x)], phi_975[order(x)],
      lty = 2,
      col = "black")
legend(-1, 35, 
       legend = c(TeX( '$\\hat{\\varphi}(x)_{\\bar{k}}$'), "95 % CI"), 
       col = c("black", "black"),
       lty = c(1, 2, 2), lwd = c(2, 1, 1),
       cex = 0.8 )

# Monte Carlo simulation ----

# Vectors defining the length of the data sets
no_indiv = c(50, 50, 100, 100)
no_time =  c(10, 20, 10, 20)

# Generate matrices to store the values of the RMSE and the IMAE
Sim_res_RMSE_1 = Sim_res_RMSE_2 = Sim_res_RMSE_3 = Sim_res_RMSE_4 = matrix( 0, nrow = rep, ncol = 4 )
colnames(Sim_res_RMSE_1) = colnames(Sim_res_RMSE_2) = colnames(Sim_res_RMSE_3) = colnames(Sim_res_RMSE_4) <- c("LF/LMU", "LF", "LMU", "LL")

Sim_res_IMAE_1 = Sim_res_IMAE_2 = Sim_res_IMAE_3 = Sim_res_IMAE_4 = matrix( 0, nrow = rep, ncol = 4 )
colnames(Sim_res_IMAE_1) = colnames(Sim_res_IMAE_2) = colnames(Sim_res_IMAE_3) = colnames(Sim_res_IMAE_4) <- c("LF/LMU", "LF", "LMU", "LL")

dgp = function(x){ x^2 }

# Monte Carlo repetitions
rep = 400

for(d in 1:4){
  for(i in 1:rep){ 
    # Specify length of the data
    N = no_indiv[d]; T = no_time[d] ; NT = N*T
    # Draw the data 
    Data = dgpivfe(N=N, T=T)
    
    # 1) The LF/LMU estimator (IV with two-way fixed effects)
    if(i == 1){ # Bandwidths are fixed to those obtained of the first repetition 
      res.npivfe = npregivfe( y = Data$y,
                              x = Data$x, 
                              z = Data$z,
                              N = N,
                              T = T,
                              c = 0.8, 
                              tol = 1,
                              max.iter = 100,
                              bw_method = "optimal" )
      bws_LF_LMU = phi_hat_LF_LMU[[3]]
      fit.npivfe =res.npivfe[[1]][, ncol( res.npivfe[[1]] ) ]
    }else{
      res.npivfe = res.npivfe( y = Data$y,
                               x = Data$x,
                               z = Data$z,
                               N = N, 
                               T = T,
                               c = 0.8,
                               tol=1,
                               max.iter=100,
                               bw_method = "plug_in",
                               bws = bws_LF_LMU )
    }
    # 2) The LF estimator (nonparametric IV only)
    if(i == 1){
      res.npiv2 = npregiv2( y = Data$y,
                            z = Data$x, 
                            w = Data$z, 
                            c = 0.8,
                            tol=1,
                            max.iter = 100,
                            bw_method = "optimal" )
      bws_LF = res.npiv2[[3]]
      fit.npivf2 = res.npivf2[[1]][, ncol( res.npivf2[[1]] ) ]
    }else{
      res.npiv2 = npivf2( y = Data$y,
                          z = Data$x,
                          w = Data$z,
                          c = 0.8, 
                          tol=1,
                          max.iter = 100,
                          bw_method = "plug_in",
                          bws = bws_LF  )
      fit.npivf2 = res.npivf2[[1]][, ncol(res.npivf2[[1]])]
    }
    
    # The LUM estimator (nonparametric two-way fixed effects) 
    if(i==1){
      bw_LMU = optimize( f = Lee.CVCMB,
                         interval = c(0,5*sd(Data$x)),
                         lower = 0,
                         upper = 5*sd(Data$x),
                         tol = 0.001,
                         x = Data$x,
                         y = Data$y,
                         N = N,
                         T = T,
                         effects = "indiv-time",
                         type = "ll",
                         jump = 1,
                         c1 = 1/T^2, c2 = 1/((N*T)^2))$minimum
      
      fit.LMU = LMU.estimator( h = bw_LMU,
                               x = Data$x,
                               y = Data$y,
                               N = N,
                               T = T,
                               effects = "indiv-time",
                               type = "ll",
                               c1 = 1/T^2, c2 = 1/((N*T)^2))[,"mhat"]
    }else{
      fit.LMU = LMU.estimator( h = bw_LMU,
                               x = Data$x, 
                               y = Data$y,
                               N = N,
                               T = T,
                               effects = "indiv-time",
                               type = "ll",
                               c1 = 1/T^2, c2 = 1/((N*T)^2))[,"mhat"]
    }
    
    # The LL estimator (local-linear kernel regression) 
    if(i==1){
      bw_ll = npregbw( Data$y~Data$x,
                       ckertype="gaussian" ,
                       regtype = "ll",
                       bwmethod = "cv.ls")$bw 
      fit.LL = fitted( npreg( bws = bw_ll,
                              tydat = Data$y, 
                              txdat = Data$x ) ) 
    }else{
      fit.LL = fitted( npreg( bws = bw_ll,  
                              tydat = Data$y,
                              txdat = Data$x,
                              ckertype="gaussian"))
    }
    if(d == 1){ 
      # RMSE 
      Sim_res_RMSE_1[i,"LF/LMU"] = sqrt( mean( ( fit.npivfe - dgp(Data$x ) )^2 ) )
      Sim_res_RMSE_1[i,"LF"] = sqrt( mean( ( fit.npiv2 - dgp( Data$x ) )^2 ) )
      Sim_res_RMSE_1[i,"LMU"] = sqrt( mean( ( fit.LMU - dgp( Data$x ) )^2) )
      Sim_res_RMSE_1[i,"LL"] = sqrt( mean( ( fit.LL - dgp(Data$x) )^2) )
      # IMEA 
      Sim_res_IMAE_1[i,"LF/LMU"] = mean( abs( fit.npivfe - dgp( Data$x ) ) )
      Sim_res_IMAE_1[i,"LF"] = mean( abs( fit.npiv2 - dgp( Data$x ) ) )
      Sim_res_IMAE_1[i,"LMU"] = mean( abs( fit.LMU - dgp( Data$x ) ) )
      Sim_res_IMAE_1[i,"LL"] = mean( abs( fit.LL - dgp( Data$x ) ) )
    }
    if(d==2){
      # RMSE 
      Sim_res_RMSE_2[i,"LF/LMU"] = sqrt( mean( ( fit.npivfe - dgp(Data$x ) )^2 ) )
      Sim_res_RMSE_2[i,"LF"] = sqrt( mean( ( fit.npiv2 - dgp( Data$x ) )^2 ) )
      Sim_res_RMSE_2[i,"LMU"] = sqrt( mean( ( fit.LMU - dgp( Data$x ) )^2) )
      Sim_res_RMSE_2[i,"LL"] = sqrt( mean( ( fit.LL - dgp(Data$x) )^2) )
      # IMEA 
      Sim_res_IMAE_2[i,"LF/LMU"] = mean( abs( fit.npivfe - dgp( Data$x ) ) )
      Sim_res_IMAE_2[i,"LF"] = mean( abs( fit.npiv2 - dgp( Data$x ) ) )
      Sim_res_IMAE_2[i,"LMU"] = mean( abs( fit.LMU - dgp( Data$x ) ) )
      Sim_res_IMAE_2[i,"LL"] = mean( abs( fit.LL - dgp( Data$x ) ) )
    }
    if(d==3){
      # RMSE 
      Sim_res_RMSE_3[i,"LF/LMU"] = sqrt( mean( ( fit.npivfe - dgp(Data$x ) )^2 ) )
      Sim_res_RMSE_3[i,"LF"] = sqrt( mean( ( fit.npiv2 - dgp( Data$x ) )^2 ) )
      Sim_res_RMSE_3[i,"LMU"] = sqrt( mean( ( fit.LMU - dgp( Data$x ) )^2) )
      Sim_res_RMSE_3[i,"LL"] = sqrt( mean( ( fit.LL - dgp(Data$x) )^2) )
      # IMEA 
      Sim_res_IMAE_3[i,"LF/LMU"] = mean( abs( fit.npivfe - dgp( Data$x ) ) )
      Sim_res_IMAE_3[i,"LF"] = mean( abs( fit.npiv2 - dgp( Data$x ) ) )
      Sim_res_IMAE_3[i,"LMU"] = mean( abs( fit.LMU - dgp( Data$x ) ) )
      Sim_res_IMAE_3[i,"LL"] = mean( abs( fit.LL - dgp( Data$x ) ) )
    }
    if(d==4){
      # RMSE 
      Sim_res_RMSE_4[i,"LF/LMU"] = sqrt( mean( ( fit.npivfe - dgp(Data$x ) )^2 ) )
      Sim_res_RMSE_4[i,"LF"] = sqrt( mean( ( fit.npiv2 - dgp( Data$x ) )^2 ) )
      Sim_res_RMSE_4[i,"LMU"] = sqrt( mean( ( fit.LMU - dgp( Data$x ) )^2) )
      Sim_res_RMSE_4[i,"LL"] = sqrt( mean( ( fit.LL - dgp(Data$x) )^2) )
      # IMEA 
      Sim_res_IMAE_4[i,"LF/LMU"] = mean( abs( fit.npivfe - dgp( Data$x ) ) )
      Sim_res_IMAE_4[i,"LF"] = mean( abs( fit.npiv2 - dgp( Data$x ) ) )
      Sim_res_IMAE_4[i,"LMU"] = mean( abs( fit.LMU - dgp( Data$x ) ) )
      Sim_res_IMAE_4[i,"LL"] = mean( abs( fit.LL - dgp( Data$x ) ) )
    }
    
  }
}

# Output table of the MC simulation
tab.MC = cbind(no_indiv, no_time, rbind( round( colMeans( Sim_res_RMSE_1 ), 3 ), 
                                         round( colMeans( Sim_res_RMSE_2 ), 3 ),
                                         round( colMeans( Sim_res_RMSE_3 ), 3 ), 
                                         round( colMeans( Sim_res_RMSE_4 ), 3 ) ),
               rbind( round( colMeans( Sim_res_IMAE_1 ), 3 ), 
                      round( colMeans( Sim_res_IMAE_2 ), 3 ),
                      round( colMeans( Sim_res_IMAE_3 ), 3 ), 
                      round( colMeans( Sim_res_IMAE_4 ), 3 ) ) )

colnames(tab.MC)[1:2] <- c("N", "T")

