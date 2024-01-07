Readme file for the package "npivfe.R" to conduct nonparametric instrumental regression with two-way fixed effects.
See De Monte (2023), Nonparametric instrumental regression with two-way fixed effects, Journal of Econometric Methods. 


--------------------------------------------------------------------------------------------

1) 
# LMU.CVCMB( h, 
             x, 
             y, 
             N,
             T, 
             effects = c("individual", "time", "indiv-time"
             type = c("ll", "lq", "l_cubic") , 
             jump,   
             c1, c2 ) 


# The function LMU.CVCMB allows to compute conitional mean based (CMB) 
# leave-one-out cross-validated (CV) bandwidths (see Appendix A in De Monte (2023)).
# The function usually is used to be minimized w.r.t. h the bandwidh parameter

Notes: 
# The function requires the MASS package.
# Applicable only for balanced panel data. 

# Arguments 
# h := kernel bandwidth parameter (should be treated a free parameter for optimal bandwith slection) 
# x := (univariate) explanatory variable (NT x 1) vector
# y := dependent variable (NT x 1) vector
# N := number of individuals 
# T := number of time periods 
# effect: i) "individual" if panel model with only individual fixed effects 
#         ii) "time" if panel model with only temporal fixed effects 
#         ii) "indiv-time" if panel model with two-way fixed effects 
# type: "ll" local-linear method will be used
#        "lq" local-quadratic method will be used
#        "l-cubic" local-cubic method will be used
# jump := allows to speed up by leaving out more than one, scalar
# c1 and c2 : constants to be chooses 
# (usually:  c1 =  1/T^2 ; c2 = 1/((N*T)^2)

Output
# Scalar. That is, the function yields the error of the cross-validation function,
#  which usually should be minimized. 
--------------------------------------------------------------------------------------------

2) 
# LMU.estimator = function( h, 
                            y, 
                            x,
                            N, 
                            T,
                            effects = c("individual", "time", "indiv-time"),
                            type = c("ll", "lq", "l_cubic"),
                            c1, c2)

# The function implements the estimator described in 
# Lee, Y., Mukherjee, D., & Ullah, A. (2019). Nonparametric estimation of the marginal effect in fixed-effect panel data models. 
# Journal of Multivariate Analysis, 171, 53-67. and in De Monte (2023) (for obtaining the conditional mean function). 

Notes: 
# The function requires the MASS package.
# Applicable only for balanced panel data. 

# Arguments 
# h := kernel bandwidth parameter, scalar 
# x := (univariate) explanatory variable, (NT x 1) vector
# y := dependent variable, (NT x 1) vector
# N := number of individuals, scalar
# T := number of time periods, scalar 
# effect: i) "individual" if panel model with only individual fixed effects 
#         ii) "time" if panel model with only temporal fixed effects 
#         ii) "indiv-time" if panel model with two-way fixed effects 
# type: "ll" local-linear method will be used
#        "lq" local-quadratic method will be used
#        "l-cubic" local-cubic method will be used
# c1 and c2 : constants to be chooses 
# (usually:  c1 =  1/T^2 ; c2 = 1/((N*T)^2)

# Output 
# NT x (1+D) matrix, with D the number of computed derivatives of the conditional mean function.
# The first columnt "mhat" contains the pointwise estimates of the conditional mean function. 
# The subsequent columns "gradient.1" (if type == "ll"), "gradient.2" (if type == "lq"),
# gradient.3 (if type == "l_cubic) provide the gradient(s) of the conditional mean function. 
--------------------------------------------------------------------------------------------

3) 
# npregivfe(y,
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
            bws ) 

# The Landweber-Fridman / LMU estimator 
# Nonparametric instrumental regression with two-way fixed effects, 
# De Monte (2023), Journal of Econometric Methods 

# Notes: 
# The function requires the functions LMU.estimator() and LMU.CVCMB(). 
# Applicable for balanced panel data only.
 
# Arguments: 
# y := dependent variable, (NT x 1) vector
# x := endogenous explanatory variable, (NT x 1) vector
# z := instrument (NT x 1) vector
# c := tuning parameter (0  < c < 1), scalar
# tol := tolerance level stopping criterion (in %), scalar
# N := number of individuals, scalar 
# T := number of time periods, scalar 
# max.iter := maximum number of iterations, scalar  
# bw_method = c("optimal", "plug_in" , "silverman")
#  i) "optimal" uses leave-one-out cross validation 
#  ii) "plug_in" provide bandwidth vector of five bandwidths  
#  iii) "silverman" computes and uses Silverman rule of thump bandwidths 
# effect = c("individual" , "time", "indiv-time") 
# i) "individual" if panel model with only individual fixed effects 
# ii) "time" if panel model with only temporal fixed effects 
# ii) "indiv-time" if panel model with two-way fixed effects 
# type = c("ll", "lq", "l_cubic")  
# "ll" local-linear method will be used
#  "lq" local-quadratic method will be used
# "l_cubic" local-cubic method will be used
# bws := (5 x 1) vector of plug-in bandwiths (five reguired) if bw_method == "plug_in"

# Output: 
# list object "output"
# output[[1]] := contains the estimates of the regularized solution path
# i.e. the conditional mean estimates phi(x)_k_hat until convergence of the Landweber-Fridman procedure 
# The final estimate of the conditional mean function can be obtained by 
# phi_hat = output[[1]][,ncol(output[[1]])]
# output[[2]] := contains the values of the function of the stopping criterion 
# ouput[[3]] := contains the used/estimated bandwidths 

--------------------------------------------------------------------------------------------
4) 
# npregiv2( y,
            z,
            w, 
            c,
            tol, 
            max.iter, 
            bw_method = c("optimal","plug_in", "sylverman"),
            bws )

# Nonparametric iv regression

# Notes: 
# The function requires the np package.
# Applicable for pooled data. 
# Also see npregiv from the np package, which should yield similar results.

# Arguments: 
 y := dependent variable, (NT x 1) vector
# x := endogenous explanatory variable, (NT x 1) vector
# z := instrument (NT x 1) vector
# w := tuning parameter (0  < c < 1), scalar
# tol := tolerance level stopping criterion (in %), scalar 
# max.iter := maximum number of iterations of the Landweber-Fridman regularization, scalar
## bw_method = c("optimal", "plug_in" , "silverman")
#  i) "optimal" uses leave-one-out cross validation 
#  ii) "plug_in" provide bandwidth vector of five bandwidths  
#  iii) "silverman" computes and uses Silverman rule of thump bandwidths 
# bws := (5 x 1) vector of plug-in bandwiths (five reguired) if bw_method == "plug_in"

# Output: 
# list object "output"
# output[[1]] := contains the estimates of the regularized solution path
# i.e. the conditional mean estimates phi(x)_k_hat until convergence 
# of the Landweber-Fridman procedure 
# The final estimate of the conditional mean function can be obtained by 
# phi_hat = output[[1]][,ncol(output[[1]])]
# output[[2]] := contains the values of the function of the stopping criterion 
# ouput[[3]] := contains the used/estimated bandwidths

--------------------------------------------------------------------------------------------

5) 
# dgpnpivfe( N, 
             T )

# The function allows to generate panel data from the data generating process 
# described in the paper. 

Note: The function requires the package data.table.

Arguments:
N := Number of individuals 
# := Number of time periods

# Output
# data.table dataset with N*T number of rows and 5 columns: The individual and time identifiers col_N and col_T
# as well as the depenent variable y, the endogenous explanatory variable x, and the instrument z.
# Note there is no seed specified. To obtain the same data as that used in the paper, specify set.seed(49). 


--------------------------------------------------------------------------------------------
  