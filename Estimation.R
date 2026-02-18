###########Goodness of fit indicators (BIC, Wasserstein, qqplots)#####

#Mean empirical wasserstein distance between 100 simulated pops and the data for each model

#' Computes the mean empirical wassersein distance between simulated populations 
#' for a given model and the data
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param X1 the data for time spent non smurf, an array
#' @param X2 the data for time spent smurf, an array with same indexes as X1
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#'  @param n an integer, the populations are simulated for n*10 
#'  default is 150 
#'  @returns an array with the mean empirical wasserstein distance on 100 simulations
#'  and the associated standard deviation
emp_wasserstein <- function(p,i,X1,X2,paramS,paramD,n=150){
  w = 0 
  wsd = 0
  N = length(X1)
  init_size = N
  pop_init <- population(data.frame( birth = rep(0, init_size), 
                                     death = NA, smurf = rep(FALSE,init_size), 
                                     time_smurf = rep(0, init_size)),id = TRUE)
  
  smurf_param = get_code_ks(p,paramS) #smurf rate parameters and code
  code_smurf = smurf_param$code
  paramS1 = smurf_param$par
  smurf = mk_event_individual(type = 'swap', 
                              intensity_code = paste('if (I.smurf) {result = 0;} 
                                      else {result = ',code_smurf,';}'),
                              kernel_code = 'I.smurf = true; I.time_smurf = t;'
  ) #defining smurf event
  death_param = get_code_kd(i,paramD) #death rate parameters and code
  code_death = death_param$code
  paramD1 = death_param$par
  params = c(paramS1,paramD1)
  death = mk_event_individual(type = 'death', 
                              #kernel_code = 'I.death = t;',
                              intensity_code = code_death
  )#defining death event
  birth_death <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death,smurf),
                          parameters = params)
  for (k in 1:100){
    sim_out <- popsim(model = birth_death,
                      initial_population = pop_init,
                      events_bounds = c('death' = 1, 'swap' = 1),
                      parameters = params,
                      time = (0:((n-1)))*10) #simulation
    Test_pop_h = sim_out$population #simulate population
    TD = Test_pop_h[[n-1]]$death
    TNS=Test_pop_h[[n-1]]$time_smurf
    mass1 = wpp(matrix(c(X2,X1),N,2), rep(1,N)) #experimental empirical distribition
    mass2 = wpp(matrix(c(TD-TNS,TNS),N,2), rep(1,N)) #simulated empirical distribution 
    w = w +  transport::wasserstein(mass1,mass2) #distance
    wsd = wsd +  (transport::wasserstein(mass1,mass2))^2 #square distance for standard deviation
  }
  w = w/100
  sd= sqrt(-w^2 + wsd/100) #standard deviation
  return(c(w,sd))
}

#' Computes the mean empirical wassersein distance between simulated populations 
#' for a given model and the data
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param X1 the data for time spent non smurf, an array
#' @param X2 the data for time spent smurf, an array with same indexes as X1
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#'  @param n an integer, the populations are simulated for n*10 
#'  default is 150 
#'  @returns an array with the mean empirical wasserstein distance on 100 simulations
#'  and the associated standard deviation
emp_wasserstein1D <- function(p,X1,paramS,n=200,K=100){
  w = 0 
  wsd = 0
  N = length(X1)
  init_size = N
  pop_init <- population(data.frame( birth = rep(0, init_size), 
                                     death = NA))
  
  smurf_param = get_code_ks(p,paramS) #smurf rate parameters and code
  code_smurf = smurf_param$code
  paramS1 = smurf_param$par
  params = c(paramS1)
  death = mk_event_individual(type = 'death', 
                              intensity_code = paste('{result = ',code_smurf,';}')
  )#defining death event
  mod <- mk_model(characteristics = get_characteristics(pop_init),
                  events = list(death),
                  parameters = params)
  for (k in 1:K){
    sim_out <- popsim(model = mod,
                      initial_population = pop_init,
                      events_bounds = c('death' = 1),
                      parameters = params,
                      time = (0:((n-1)))*10) #simulation
    Test_pop_h = sim_out$population #simulate population
    TD = Test_pop_h[[n-1]]$death
    w = w +  transport::wasserstein1d(TD,X1) #distance
    wsd = wsd +  (transport::wasserstein1d(TD,X1))^2 #square distance for standard deviation
  }
  w = w/K
  sd= sqrt(-w^2 + wsd/K) #standard deviation
  return(c(w,sd))
}



##BIC

#'Computes the BIC of a smurf or death jumping time model for a given hazard rate
#'@param X the jumping times for which we want to compute the BIC
#'@param Y an eventual covariable (smurf jumping times when computing the BIC for death jumping times)
#'@param haz the hazard rate of the density we want to compute the BIC for, a function of 2 variables, 
#'the jumping time and a covariable (can be NA)
#'@param k the number of parameters of hazard rate haz
#'@returns a number, the BIC of model with hazard rate haz for jumping times X
BIC <- function(X,Y,haz,k){
  fun = function(t,i)(haz(t,Y[i]))
  log_v = sum(sapply(1:length(X), function(i)(log(haz(X[i], Y[i]))))) - sum(sapply((1:length(X)), function(i) (integrate(function(t)(fun(t,i)),0,X[i])$value)) )
  return(k*log(length(X)) - 2*log_v)
} 

#'Computes the BIC of a smurf or death jumping time model for a given density function
#'@param X the jumping times for which we want to compute the BIC
#'@param Y an eventual covariable (smurf jumping times when computing the BIC for death jumping times)
#'@param dens the density for which we want to compute the BIC for, a function of 2 variables, 
#'the jumping time and a covariable (can be NA)
#'@param k the number of parameters of density dens
#'@returns a number, the BIC of model with hazard rate haz for jumping times X
BIC_dens <- function(X,Y,dens,k){
  fun = function(t,i)(dens(t,Y[i]))
  log_v = sum(sapply(1:length(X), function(i)(log(dens(X[i], Y[i])))))
  return(k*log(length(X)) - 2*log_v)
} 


##QQ plots time spent non-smurf

#' Returns an array of death or smurf jumping times for a simulated population 
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#'  @param case s or d depending on whether the function returns the time spent 
#'  smurf or non smurf
#'  @param N an integer, the population size, default is 1159
#'  @param n an integer, the populations are simulated for n*10 
#'  default is 150 
#'  @returns an array of simulated death or smurf jumping times 
plot_qq_pop <- function(p,i,paramS,paramD,case,N=1159, n=150){
  Test_pop = sim_pop_n(p,i,paramS,paramD,N,n)#simulate population
  TD = Test_pop[[n-1]]$death
  TNS=Test_pop[[n-1]]$time_smurf
  if(case =='s'){
    return(TNS)
  }
  else{
    return(TD-TNS)
  }
}


###########Non-parametric associated kernel estimation#####
##Defining some kernels and kernel estimators##

#' Gamma kernel 'of y' of bandwidth h estimated at x
#' (parameter of gamma distribution depend on y)
#' @param y array of jumping times
#' @param x value at which it is estimated
#' @param h bandwidth
#' @returns value of y,h-dependent gamma kernel evaluated at point x
gamma_kernel  = function(y, x, h) {
  indices = (y >= 2*h)
  y_new = y/h*indices + (1/4*(y/h)^2 + 1)*(!indices)
  h*dgamma(x, y_new , scale = h)
} #gamma kernel for non-parametric estimation

gamma_ker = function(y,x,h)(gamma_kernel(y,x,h)/h)


#' Gamma kernel estimator of hazard rate of jumping times in X with fixed bandwidth
#' h at point t
#' @param t point at which to evaluate kernel estimation
#' @param X array of jumping times
#' @param h bandwidth
#' @returns gamma kernel estimator with fixed bandwidth h of hazard rate of jumping times X at t 
ker_est_gamma_c <- function(t,X,h){
  return(sum(gamma_kernel(t, sort(X),h)/(length(X):(1+length(X)-length(sort(X)))))/h)
} # gamma kernel constant bandwidth




##Minimax bandwidth estimation##


#'Function V_0 for local minimax bandwidth choice
#'@param t time of estimation 
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns the value of V_0 at t for the bandwidth b (estimation of the variance 
#'of the estimator) 
V_0 = function(t,b,B,max,m, kappa_0 = NULL, lambda = NULL){
  if(is.null(kappa_0)){kappa_0 = 0.03}
  if(is.null(lambda)){lambda = 4}
  return(kappa_0*max*log(m)/(m*sqrt(b))) 
}

#'Function A_0 for local minimax bandwidth choice
#'@param t time of estimation 
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param Ttest array of jumping times for hazard estimation
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns the value of A_0 at t for the bandwidth b (estimation of the bias
#'of the estimator) 
A_0 = function(t,b,B,max,m,Ttest,kappa_0 = NULL, lambda = NULL){
  A = c()
  for (i in 1:length(B)){
    inter = (ker_est_gamma_c(t,Ttest,B[i]) - ker_est_gamma_c(t,Ttest,max(b,B[i])))^2 - V_0(t,B[i],B, max,m,kappa_0 , lambda )
    A[i] = inter}
  return(max(A,0))   
}

#'Function V for global minimax bandwidth choice
#'@param Grid array on which to perform estimation
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns the value of V on Grid for the bandwidth b (estimation of the variance 
#'of the estimator) 
V = function(Grid,b,B,max,m, kappa_2 = NULL, lambda = NULL, epsilon = NULL) {
  if(is.null(kappa_2)){kappa_2 = 20}
  if(is.null(lambda)){lambda = 4}
  if(is.null(epsilon)){epsilon = 0.5}
  return((1+epsilon)^2*kappa_2*max/(m*sqrt(b)))
}

#'Function A for global minimax bandwidth choice
#'@param Grid array on which to perform estimation
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param Ttest array of jumping times for hazard estimation
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns the value of A on Grid for the bandwidth b (estimation of the bias
#'of the estimator) 
A = function(Grid,b,B,max,m,Ttest,kappa_2 = NULL, lambda = NULL, epsilon = NULL){
  A = c()
  for (i in 1:length(B)){
    integrand = function(t)((ker_est_gamma_c(t,Ttest,B[i]) - ker_est_gamma_c(t,Ttest,max(b,B[i])))^2)
    integrand = Vectorize(integrand)
    inter = integrate(integrand,0, max(Grid),subdivisions=2000)$value - V(Grid,B[i],B, max,m,kappa_2, lambda, epsilon )
    A[i] = inter}
  return(max(A,0))   
}



#'Returns a set of bandwidth depending for adaptive bandwidth choice
#'for a given sample size
#'@param m a sample size
#'@returns a set of bandwidths
Bandwidth_set <- function(m){
  B = c()
  B[1] = 400*(log(m)/m)^2
  i = 1
  while (i <= 20*log(m)){
    B[(i-1)/4+1] = log(m)^2*i/m
    i = i+4
  }
  B = B[!is.na(B)]
  B = B[sqrt(B) <= min(1,6/log(m))]
  B = B[sqrt(B) >= B[1]]
  return(B)
}



#'Returns a set of bandwidth depending for adaptive bandwidth choice
#'for a given sample size for the global bandwidth choice procedure
#'@param m a sample size
#'@returns a set of bandwidths
Bandwidth_set_global <- function(m){
  B = c()
  B[1] = 1/m^(2/3)
  i = 1
  while (i <= (10*sqrt(m))){
    B[(i-1)/10+1] = i/m^(2/3)
    i = i+10
  }
  B = B[!is.na(B)]
  B = B[sqrt(B) <= min(1,6/log(m))]
  B = B[sqrt(B) >= sqrt(B[1])]
  return(B)
}


#'Gives the adaptive bandwidth estimator of the hazard rate for a data sample on  
#'a Grid and the chosen bandwidth for each point of the Grid
#'@param Times the data sample 
#'@param Grid the grid for the estimation 
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@param B (optional) a bandwidth set
#'@returns a list with two arrays : K, the estimated hazard rate on Grid 
#', B, the chosen bandwidth for each point of Grid 
minimax_pointwise_data <- function(Times, Grid,kappa_0 = NULL,lambda = NULL,B = NULL){
  m = length(Times)
  if(is.null(B)){B = Bandwidth_set(m)}
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated
  #upper bound of hazard rate
  Kopt = c()
  Bopt=c()
  for (j in 1:length(Grid)){ 
    t = Grid[j]
    C = sapply(B, function(b)(A_0(t,b,B,max,m,Times,kappa_0, lambda) + V_0(t,b,B,max,m,kappa_0, lambda)))
    bopt = B[which.min(C)]
    Bopt[j] = bopt #local minimax bandwidth
    Kopt[j] = ker_est_gamma_c(t,Times,bopt) #local estimator value
  }
  
  return(list('K' = Kopt, 'B'= Bopt))
}

#'Gives the adaptive bandwidth estimator of the hazard rate for a data sample on  
#'a Grid and the chosen global bandwidth 
#'@param Times the data sample 
#'@param Grid the grid for the estimation 
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@param B (optional) a bandwidth set
#'@returns a list with two arrays : K, the estimated hazard rate on Grid 
#', B, the chosen bandwidth 
minimax_global_data <- function(Times, Grid,kappa_2=NULL, lambda=NULL,epsilon=NULL,B = NULL){
  m = length(Times)
  if(is.null(B)){B = Bandwidth_set_global(m)}
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated
  #upper bound of hazard rate
  Bopt=c()
  C = sapply(B, function(b)(A(Grid,b,B,max,m,Times,kappa_2, lambda,epsilon) + V(Grid,b,B,max,m,kappa_2, lambda,epsilon)))
  bopt = B[which.min(C)] #global minimax bandwidth
  Kopt = sapply(Grid, function(t)(ker_est_gamma_c(t,Times,bopt))) #global minmax bandwidth gamma estimator
  return(list('K' = Kopt, 'B'= bopt))
}


#'Monte-Carlo integral approximation
#'@param fun function to integrate (of one variable)
#'@param A upper bound of integration
int_MC<- function(fun,A,n){
  samp = A*runif(n,0,1)
  res = A*mean(sapply(samp, fun))
  return(res)
}

#'CDF of a random variable with hazard rate haz_est, computed with monte-carlo
#'@param t real number point at which the CDF is given
#'@param haz_est function of a single variable (time)
#'@param n integer, number of repetitions for the monte-carlo estimation
#'@returns a real number, the cdf associated to the hazard rate haz_est, evaluated
#'at t
CDF_M_C <- function(t,haz_est,n){
  k_int = int_MC(haz_est,t, n)
  return(1-exp(-k_int))
}


##Variance of estimator
#'Variance of the kernel estimator evaluated at point t
#'@param t real number point at which the variance is evaluated
#'@param m number of data points
#'@param ker kernel function of three variables, the point of evaluation, the 
#'point on which the shape of the kernel depends and the bandwidth 
#'@param haz hazard rate, function of a single variable (time)
#'@param n integer, number of repetitions for the monte-carlo estimation
#'@param b positive bandwidth for the kernel 
#'@param max bound to use for the monte-carlo simulation to compute the integral 
#'on [0,max] (max should be as big as possible)
#'@returns a real number, the variance of the estimator haz or m data points
#' with kernel ker, evaluated at t
Var_M_C <- function(t,m, ker, haz,n,b,max){
  CDF = function(y)(CDF_M_C(y,haz,n))
  integrand = function(y)(ker(t,y,b)^2)
  return(int_MC(integrand, max,n)*haz(t)/(1-CDF(t))*1/m)
}



###########Parametric estimation of Smurf transition rate######

#Max likelihood estimation for c*exp(d*t)

#' Compute derivative of log likelihood of time spent non smurf with gompertz smurf 
#' hazard rate c.exp(d.t), at a given value of d 
#' @param d a parameter value
#' @param X1 an array of time spent non-smurf
#' @returns the derivative at d 
Der_vrais_exp<- function(d,X1){#derivative of the likelihood wrt d
  N = length(X1)
  c = (N/d+sum(X1))*d/(sum(X1*exp(d*X1)))#parameter c for which derivative wrt
  #c is 0
  return(N- c/d*(sum(exp(d*X1)) -N))
}


#max likelihood of ks = f + g*exp(h*t)
#'Find max likelihood estimations of gompertz makeham smurf hazard rate f + g.exp(h.t)
#' @param min min for h value
#' @param max max for h value
#' @param X1 an array of time spent non-smurf
#' @returns an array with parameters f, g and h of gompertz makeham smurf hazard rate 
param_Gomp_Mak <- function(X1,min,max){
  N = length(X1)
  X1ord = sort(append(0,X1,0))#reorder times
  NS = (N-1):0 #empirical survival function
  
  der_vrais_h <- function(h){ #subfunction to find the optimal h
    int = sum((exp(h*X1ord[2:(N+1)]) - exp(h*X1ord[1:N]))*NS)/h
    xm = min(X1ord[2:N])
    find_0 <- function(j){#function to determine parameter g by computing the 
      #derivative wrt f and max lik f for j and h
      f = (N-j*int)/sum((X1ord[2:(N+1)] - X1ord[1:N])*NS) #max lik f for this value of j and h
      #f is a decreasing function of j
      sum1 = sum(1/(f+j*exp(h*X1))) - sum((X1ord[2:(N+1)] - X1ord[1:N])*NS) #derivative wrt f
      #return(sum(log(f+j*exp(h*X1ord))) - sum(f*X1ord+j/h*exp(h*X1ord) - j/h)) 
      #return log-likelihood to check that it is max, commented 
      return(sum1) 
    }
    glim = N/(int - sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)*exp(h*xm)) #rough bound to determine where to look for g
    A_test = (1:99)*(glim - glim/10000)/(99.5) +glim/10000 #grid from 0 to glim
    apply_test = sapply(A_test, find_0) #apply find_0 to A_test
    i_max =  max(which(apply_test < 0)) #rough estimation of where find_0 goes from neg to pos
    g = try(uniroot(find_0,c(A_test[i_max],A_test[i_max+1]),tol=0.0000001)$root,silent=TRUE) 
    #precise estimation of g where the derivative is 0
    if(inherits(g,'try-error')){#if there is no such g, not return an error 
      #because it may just be the wrong value of h
      return(NA)
    } 
    int2 = sum((X1ord[2:(N+1)]*exp(h*X1ord[2:(N+1)]) - X1ord[1:N]*exp(h*X1ord[1:N]))*NS)/h 
    inttot = int2 - int/h
    f = (N-g*int)/sum((X1ord[2:(N+1)] - X1ord[1:N])*NS) #max lik f computed directly with analytic result
    sum1 = sum(X1*exp(h*X1)/(f+g*exp(h*X1))) - inttot #derivative wrt h
    return(c(sum1,g)) #return derivative and max lik g 
  }
  
  h = try(uniroot(function(h)(der_vrais_h(h)[1]),c(min,max))$root,silent=TRUE)#max lik h where derivative is 0
  if(inherits(h,'try-error')){#raise error if derivative does not cross 0 on the specified interval
    print('Problem finding h for which derivative is 0, the derivatives at the bounds are:')
    print(der_vrais_h(min))
    print(der_vrais_h(max))
  } 
  
  int = sum((exp(h*X1ord[2:(N+1)]) - exp(h*X1ord[1:N]))*NS)/h
  g = der_vrais_h(h)[2] #max lik g for the max lik h
  f = (N-g*int)/sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)#max lik f for the max lik h and g
  return(c(f,g,h))#return max likelihood parameters
}


#max likelihood kS = at^2 + b 
#' Compute derivative of log likelihood of time spent non smurf with polynomial
#'  at^2+b smurf hazard rate 
#' @param a value for a in parametrization
#' @param X1 an array of time spent non-smurf
#' @returns derivative of log likelihood at a 
der_vrais_poly <- function(a,X1){
  X1ord = sort(append(0,X1,0))
  N = length(X1)
  NS = N-1:N
  int1 = sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)
  int2 = sum((X1ord[2:(N+1)]^3 - X1ord[1:N]^3)*NS/3)
  b = (N- int2*a)/int1 #parameter b for which derivative wrt b is 0
  sum = sum(1/(a*X1ord[2:(N+1)]^2 + b))
  return( c(sum - int1,b)) 
}

##Confidence interval computation
#' Compute 95% CI for a max likelihood estimation 
#' @param param an array of parameter values for a given parametrization
#' @param loglike a function taking as input a parameter array, an array of time
#' spent non smurf and an array of time spent smurf and returning the loglikelihood
#' for a given model with the parameter array
#' @param X1 times spent non smurf
#' @param X2 times spent smurf
#' @returns an array of 95% standard error for each parameter in param
confidence_interval <- function(param,loglik,X1,X2){
  N = length(X1)
  H = numDeriv::hessian(function(p)(loglik(p,X1,X2)/N),
                        param,method.args=list(eps=1e-10, d=0.1, 
                                               zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))
  #numerical hessian
  A = expm::sqrtm(solve(-H)) #matrix square root of -H^-1
  M = abs(1.96/sqrt(N)*A)
  return(sapply(1:sqrt(length(M)), function(i)(M[i,i])))
}

#'Sum of log likelihood for gompertz makeham parametrization of ks
#'@param param an array of parameters (f,g,h) f+ g.exp(ht)
#'@param X1 array of jumping times on which to compute the log likelihood 
#'@param X2 default NULL, an eventual covariable, not used here but the parameter
#'is here for compatibility purposes 
#'@returns sum of log likelihood of gm ks parametrization on jumping times in X1
log_lik_gm <- function(param,X1,X2=NULL){
  a = param[1]
  b = param[2]
  c = param[3]
  return(sum(log(a+b*exp(c*X1)) - a*X1- b/c*exp(X1*c) +b/c))
}

#'Sum of log likelihood for polynomial parametrization of ks a.t^2 + b
#'@param param an array of parameters (a,b)
#'@param X1 array of jumping times on which to compute the log likelihood 
#'@param X2 default NULL, an eventual covariable, not used here but the parameter
#'is here for compatibility purposes 
#'@returns sum of log likelihood of poly ks parametrization on jumping times in X1
log_lik_poly <- function(param,X1,X2=NULL){
  a = param[1]
  b = param[2]
  return(sum(log(a*X1^2+b) - a*X1^3/3- b*X1))
}

#'Sum of log likelihood for gompertz parametrization of ks c.exp(dt)
#'@param param an array of parameters (c,d)
#'@param X1 array of jumping times on which to compute the log likelihood 
#'@param X2 default NULL, an eventual covariable, not used here but the parameter
#'is here for compatibility purposes 
#'@returns sum of log likelihood of g ks parametrization on jumping times in X1
log_lik_g <- function(param,X1,X2=NULL){
  a = param[1]
  b = param[2]
  return(sum(log(a*exp(b*X1)) - a/b*exp(b*X1) + a/b))
}

#'Sum of log likelihood for weibull parametrization of ks with k shape and l scale
#'parameters
#'@param param an array of parameters (k,l) shape / scale 
#'@param X1 array of jumping times on which to compute the log likelihood 
#'@param X2 default NULL, an eventual covariable, not used here but the parameter
#'is here for compatibility purposes 
#'@returns sum of log likelihood of w ks parametrization on jumping times in X1
log_lik_w <- function(param,X1,X2=NULL){
  k = param[1]
  l = param[2]
  N = length(X1)
  return(N*log(k/l)+sum((k-1)*log(X1/l) - (X1/l)^k))
}


########### Parametric estimation of death rate once Smurf#####
#Max likelihood constant kD with and without dependence

#' Computes max likelihood estimations of kd
#' for constant parametrization of death hazard rate
#' with dependence
#' @param X2 array of time spent smurf 
#' @param X1 array of time spent non smurf
#' @param g cox model coefficient, can be obtained with
#' as.numeric(coxph(Surv(X2,rep(1,length(X2)))~X1,list(X1,X2))$coef)
#' @returns a number, the estimator for the constant parameter
const_ML_kd_dep <- function(X2,X1,g){
  N = length(X2)
  m = mean(X1)
  X2p = sort(X2,index.return = TRUE) #reorder X2
  index = X2p$ix #indexes or reorder
  X2p = X2p$x #reordered vector 
  N_S = N - ( 0:length(X1))
  intk =  sum(exp(g*(X1[index]-m)) *X2p)
  kD2 = N/intk #analytic parameter value for which derivative is 0
  return(kD2)
}

#' Computes max likelihood estimations of kd
#' for constant parametrization of death hazard rate
#' with dependence only on TNS>lim
#' @param X2 array of time spent smurf 
#' @param X1 array of time spent non smurf
#' @param lim limit such that there is dependence only on TNS >= lim
#' @param g cox coefficient on TNS >= lim 
#' @returns a number, the estimator for the constant parameter
const_ML_kd_piec <- function(X2,X1,lim,g){
  N = length(X2)
  m = mean(X1[X1 >lim])
  X2p = sort(X2[X1 >lim],index.return = TRUE)
  index = X2p$ix
  X2p = X2p$x
  X1g = X1[X1 >lim] 
  intk =  sum(exp(g*(X1g[index]-m)) *X2p) + sum(X2[X1 <= lim])
  kD2 = N/intk #analytic parameter value for which derivative is 0
  return(kD2)
}#conditional dep on X1>lim


#' Computes max likelihood estimations of kd
#' for piecewise constant parametrization of death hazard rate
#' @param X2 array of time spent smurf 
#' @returns an array with the estimator and the 95% standard error
const_ML_kd <- function(X2){ #constant kD
  int1 = sum(X2)
  kD = length(X2)/int1 #estimator
  return(c(kD,1.96*kD/sqrt(N))) #returns estimator and confidence interval
}#no dependence

#Two constants parametric estimation 

#' Computes max likelihood estimations of kd1 and kd2 
#' for  piecewise constant (before and after 24 hrs)
#' parametrization of death hazard rate with dependence 
#' @param X1 array of time spent smurf 
#' @param X2 array of time spent non smurf
#' @param g cox coefficient 
#' @returns an array with kd1 and kd2
two_const_dep <- function(X1,X2,g){
  m = mean(X1)
  N = length(X2)
  Xord = sort(append(0,X2,0),index=TRUE) #reorder X2
  Index=sort(X2,index=TRUE)$ix #reordered indexes
  Xord = Xord$x #ordered array 
  cutoff = max(which(Xord<=24))
  kd1 = (cutoff)/(sum(Xord[2:(cutoff)]*exp(g*(X1[Index[1:(cutoff-1)]]-m))) + sum(24*exp(g*(X1[Index[(cutoff):N]]-m)))) #analytic parameter value for which derivative is 0
  kd2 = (N-(cutoff)) /( sum((Xord[(cutoff+1):(N+1)]-24)*exp(g*(X1[Index[(cutoff):N]]-m)) ) )#analytic parameter value for which derivative is 0
  kd1*1.96/sqrt(cutoff)#95 CI
  kd2*1.96/sqrt(N-cutoff)#95 CI
  return(c(kd1,kd2))
} #dependence

#' Computes max likelihood estimations of kd1 and kd2
#' for piecewise constant parametrization of death hazard rate (before and after
#' 24 hours)
#' @param X2 array of time spent smurf 
#' @returns an array with kd1, kd2, 95% standard error for kd1 and for kd2
two_const_indep <- function(X2){
  N = length(X2)
  Xord=sort(append(0,X2,0))
  cutoff2 = max(which(Xord<=24))
  kd1i = (cutoff2-1)/sum((Xord[2:cutoff2] - Xord[1:(cutoff2-1)])* (N:(N-(cutoff2-2)))) #analytic parameter value for which derivative is 0
  kd2i = (N-(cutoff2-1)) /sum((Xord[(cutoff2+1):(N+1)] - Xord[cutoff2:N])* ((N-(cutoff2-1))):1) #analytic parameter value for which derivative is 0
  ci1 = kd1i*1.96/sqrt(cutoff2) #95 CI
  ci2 = kd2i*1.96/sqrt(N-cutoff2)#95 CI
  return(c(kd1i,kd2i,ci1,ci2))
} 

#' Computes max likelihood estimations of kd1 and kd2
#' for piecewise constant parametrization of death hazard rate
#' with dependence only on TNS>lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent  smurf
#' @param lim limit such that there is dependence only on TNS >= lim
#' @param g cox coefficient on TNS >= lim 
#' @returns an array with kd1 and kd2
two_const_piec <- function(X1,X2,lim,g){
  m = mean(X1[X1>lim]) #mean for cox model
  N = length(X2)
  X1g = X1[X1>lim]
  X2g = X2[X1>lim]
  X1i = X1[X1<=lim]
  X2i = X2[X1<=lim]
  Xord = sort(append(0,X2,0),index=TRUE)#sort array
  Index=sort(X2,index=TRUE)$ix#get initial indexes of sorted array
  Xord = Xord$x #get sorted array
  cutoff = max(which(Xord<=24))
  kd1 = (cutoff)/(sum(X2g[X2g <=24]*exp(g*(X1g[X2g <=24]-m)))+sum(X2i[X2i<=24]) + sum(24*exp(g*(X1g[X2g>24]-m))) + sum(24*(X2i>24))) #analytic parameter value for which derivative is 0
  kd2 = (N-(cutoff)) /(sum((X2g[X2g >24]-24)*exp(g*(X1g[X2g >24]-m)))+sum(X2i[X2i>24] - 24)) #analytic parameter value for which derivative is 0
  return(c(kd1,kd2))
}




#' Computes max likelihood parameter estimation for 
#' gompertz makeham parametrization of death hazard rate kd1+kd2*exp(-d*t) 
#' with dependence
#' @param X1 array of time spent  smurf 
#' @param A array of values of d to test
#' @param min min value for kd1
#' @param max max value for kd1
#' @returns an array with kd1,kd2, log-likelihood value, boolean indicating 
#' an error and d
g_m_death <- function(X1, A, min, max){
  aux = function(d){#auxiliary function to call param_Gomp_Mak_D
    g = param_Gomp_Mak_D(X1,d,min,max) 
    return(c(g[3], g[4]))
  }
  der_lik = sapply(A,aux)
  if(sum(der_lik[2,]) == length(A)){ #checking that we found a 0 of the derivative
    print('None of the values of d have a derivative of the log-likelihood wrt k1 
          that goes from positive to negative, change A or min/max')
    return(NA)
  }
  der_lik = der_lik[1,] #array of likelihoods wrt d
  ind_d = which(der_lik[!is.na(der_lik)]== max(der_lik[!is.na(der_lik)])) #finding d_ei for which likelihood is maximum
  d_ei = A[ind_d] #parameter in the exponential
  if(length(d_ei)>1){ #raise a warning if max is attained for several values of d
    warning('Maximum likelihood attained for several values of d, not enough 
            numerical precision.')
  }
  d_ei = d_ei[length(d_ei)] #in case there is a maximum for several values of d,
  #just choose one
  gomp_D = param_Gomp_Mak_D(X1,d_ei,min,max)
  return(c(gomp_D,d_ei))
}

#' Computes derivative of log likelihood  of
#' gompertz makeham  parametrization of death hazard rate k1+k2*exp(-d*t) at 
#' a given value of d
#' @param X1 array of time spent smurf 
#' @param d parameter value
#' @param min min value for k1
#' @param max max value for k1
#' @returns an array with k1 and k2 and the log likelihood for that value of d 
#' and error which is 1 if an error occured and 0 else
param_Gomp_Mak_D <- function(X1,d,min,max){
  N = length(X1)
  error = 0 #to check if there is no 0 of the derivative for this value of d
  der_k1 <- function(k){ #derivative wrt k1
    k2 = (N - k*sum(X1))/(-sum(exp(-d*X1)/d) + N/d)
    sum1 = sum(1/(k+k2*exp(-d*X1))) - sum(X1)
    return(sum1)
  }
  min2 = floor(min*3500)
  max2 = floor(max*3500)+1
  A=(min2:max2)/3500 #grid to search where derivative wrt k1 crosses 0 
  app = sapply(A,der_k1) #computing derivative of log likelihood wrt a on A
  app2=sign(app) #sign of the derivative of the log likelihood wrt k1
  indi = 1
  while (indi < length(app2) && app2[indi+1] == app2[indi]){ #first change of sign 
    #of the derivative of the log likelihood
    indi = indi+1
  }
  if(indi == length(app2)){ #if no change of sign, mark this value of d as an error
    error = 1
    return(c(NA,NA,-Inf, error))
  }
  k1 = try(uniroot(der_k1,c(A[indi],A[indi+1]),tol=0.000000001)$root, silent=TRUE)
  #uniroot to find where it is 0 exactly
  if(inherits(k1,'try-error')){#should never happen
    error = 1
    return(c(NA,NA,-Inf, error))
  } 
  k2 =(N - k1*sum(X1))/(-sum(exp(-d*X1)/d) + N/d) #max likelihood value of k2 computed analytically
  return(c(k1,k2,sum(log(k1+k2*exp(-d*X1))) - sum(k1*X1 - k2/d*exp(-d*X1) +k2/d), error))
}




##functions to find max likelihood parameters in case of dependence
#' Computes derivative of log likelihood for of
#' gompertz makeham  parametrization of death hazard rate k1+k2*exp(-d*t) at 
#' a given value of d with dependence
#' @param X2 array of time spent smurf 
#' @param X1 array of time spent non smurf 
#' @param d parameter value
#' @param j depending on the model number and strata, either 12_s1, 12_s2, 9 or 6
#' @param g cox coefficient
#' @param lim limit for which there is dependence only on TNS >= lim, default 0
#' @param min min value for kd1
#' @param max max value for kd1
#' @returns an array with kd1 and kd2 and the log likelihood for that value of d
#'  and error which is 1 if an error occured and 0 else
param_Gomp_Mak_D_dep <- function(X2,X1,d,j,g,lim=0,min,max){
  N = length(X1)
  error = 0
  mean_sup = mean(X1[X1>=lim])
  XX = (X1- mean_sup) * (X1>=lim) #above the limit
  XX2 = (X1- mean(X1[X1<lim])) *(X1<lim) #under the limit
  #defining the dependence coefficient array pH depending on the model
  if(j=='12_s1'){
    pH = rep(1,length(X1[X1<lim]))
  }
  else if (j=='12_s2'){
    pH = exp( g*(XX[X1>=lim])) 
  }
  else if (j=='9'){pH = exp(g*(XX)) }
  else if (j=='6') {pH = exp(g*(X1-mean(X1)))}
  N = length(pH)
  der_k1 <- function(k){ #returns the derivative of the log-likelihood
    #wrt k1 
    k2 = (N - k*sum(X2*pH))/(N/d - sum(exp(-d*X2)*pH/d)) #analytic max lik value of k2
    sum1 = sum(1/(k+k2*exp(-d*X2))) - sum(X2*pH)
    return(sum1)
  }
  min2=floor(min*1700) + 1*(floor(min*1700)==0)
  max2 = floor(max*1700)+1
  A=(min2:max2)/1700 #Grid to search for where the derivative is 0
  XM = max(X2)
  app = sapply(A,der_k1) #derivative on grid
  app2=sign(app) #sign of the derivative of the log likelihood wrt k1
  indi = 1
  while (indi < length(app2) && app2[indi+1] == app2[indi]){ #first change of sign 
    #of the derivative of the log likelihood
    indi = indi+1
  }
  if(indi == length(app2)){
    error = 1
  }
  k1 = uniroot(der_k1,c(A[indi],A[indi+1]),tol=0.000000001)$root #finding where the derivative
  #is 0 in the interval we identified where it changed sign 
  k2 =(N - k1*sum(X2*pH))/(N/d - sum(exp(-d*X2)*pH/d))
  return(c(k1,k2,sum(log(k1+k2*exp(-d*X2))) - sum(pH*(k1*X2 - k2/d*exp(-d*X2) +k2/d)),error))
}

#' Computes max likelihood parameter estimation for 
#' gompertz makeham parametrization of death hazard rate kd1+kd2*exp(-d*t) 
#' with dependence
#' @param X2 array of time spent smurf 
#' @param X1 array of time spent non smurf 
#' @param j depending on the model number and strata, either 12_s1, 12_s2, 9 or 6
#' @param g cox coefficient
#' @param lim limit for which there is dependence only on TNS >= lim, default 0
#' @param min1 min value for kd1
#' @param max1 max value for kd1
#' @param min2 min value for d
#' @param max2 max value for d
#' @returns an array with kd1,kd2 and d
max_lik_gomp_mak_D <- function(X2,X1,j,g,lim=0,min1,max1, min2,max2){
  A =min2+ (max2-min2)*(1:1000)/1000
  aux = function(d){
    g = param_Gomp_Mak_D_dep(X2, X1,d,j,g,lim, min1,max1) 
    return(c(g[3], g[4]))
  }
  der_lik_dep = sapply(A,aux)
  if(sum(der_lik_dep[2,]== length(A))){
    print('None of the values of d have a derivative of the log-likelihood wrt k1 
          that goes from positive to negative, change min/max')
    return(NA)
  }
  der_lik_dep = der_lik_dep[1,]
  ind_d= which(der_lik_dep[!is.na(der_lik_dep)]== max(der_lik_dep[!is.na(der_lik_dep)]))
  d = A[ind_d]
  if(length(d)>1){ #raise a warning if max is attained for several values of d
    warning('Maximum likelihood attained for several values of d, not enough 
            numerical precision.')
  }
  d = d[length(d)] #in case there are several maxima 
  gomp_D_dep = param_Gomp_Mak_D_dep(X2,X1,d,j,g,lim,min1,max1)
  g1 = gomp_D_dep[1]
  g2 = gomp_D_dep[2]
  return(c(g1,g2,d))
}



#'Sum of log likelihood for exponential parametrization of kd
#'@param param an array of parameters 
#'@param X1 array of jumping times on which to compute the log likelihood 
#'@param X2 default NULL, an eventual covariable, not used here but the parameter
#'is here for compatibility purposes 
#'@returns sum of log likelihood of exp kd parametrization on jumping times in X1
log_lik_gm_d <- function(param,X1,X2=NULL){
  a = param[1]
  b = param[2]
  c = param[3]
  return(sum(log(a+b*exp(-c*X1)) - a*X1+b/c*exp(-X1*c) -b/c))
}



###########Parametric and non-parametric estimation of dependence#######
#Cox model estimation

#' Wrapper function to return just the coex coefficient computed with survival
#'  package function coxph
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf with same indexes as X1
#' @returns cox proportional hazard coefficient of X2 with covariable X1
cox_coef<- function(X1,X2){
  #returns the cox coefficient if X1 is times spent non-smurf 
  #and X2 corresponding times spent smurf
  t = list(X1,X2)
  cox = coxph(Surv(X2,rep(1,length(X2)))~X1,t)
  return(as.numeric(cox$coef))
}

##Comparing correlation 

#' Correlation coefficient on a subset of one of the variables
#'@param lim correlation is computed only on X1>lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf with same indexes as X1 
#' @returns correlation coefficient between X1 and X2 on X1>lim
cor_lim_sup <- function(lim,X1,X2){
  return(cor(X1[X1>lim],X2[X1>lim]))
}

#' Correlation coefficient on a subset of one of the variables
#'@param lim correlation is computed only on X1<lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf with same indexes as X1
#' @returns correlation coefficient between X1 and X2 on X1<lim
cor_lim_inf <- function(lim,X1,X2){
  return(cor(X1[X1<lim],X2[X1<lim]))
}

#Wald test results for significance of cox coef under or above a time spent NS limit

#'Wald test p-value for significance of cox coefficient between X1 and X2 for 
#'X1 under a certain limit
#'@param lim cox coef is computed on X1<lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf with same indexes as X1
#' @returns Wald test p-value for significance of cox coefficient between X1 
#' and X2 for X1<lim
#' 
score_dep_inf <- function(lim,X1,X2){#Wald test p-value for significance of cox model coef
  #for time spent non-smurf under lim
  sc = coxph(Surv(X2[X1 <= lim])~X1[X1 <=lim])$wald.test
  return(log( 1-pchisq(sc, 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)))
  #return(coxph(Surv(X2e[X1e <= lim])~X1e[X1e <=lim])$coef)
  
}
#'Wald test p-value for significance of cox coefficient between X1 and X2 for 
#'X1 above a certain limit
#'@param lim cox coef is computed on X1>lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf with same indexes as X1 
#' @returns Wald test p-value for significance of cox coefficient between X1 
#' and X2 for X1>lim
score_dep_sup <- function(lim,X1,X2){#Wald test p-value for significance of cox model coef
  #for time spent non-smurf over lim
  sc = coxph(Surv(X2[X1 > lim])~X1[X1>lim])$wald.test
  return(log( 1-pchisq(sc, 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)))
}


#Cox model coefficient above or under a limit of time spent ns

#'Cox coefficient between X1 and X2 for 
#'X1 under a certain limit only computed with coxph of package survival
#'@param lim cox coef is computed on X1<=lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf  with same indexes as X1
#' @returns coex coefficient computed only on X1<=lim
g_dep_inf <- function(lim,X1,X2){
  sc = coxph(Surv(X2[X1 <= lim])~X1[X1 <=lim])$coef
  return(sc)
}

#'Cox coefficient between X1 and X2 for 
#'X1 above a certain limit only computed with coxph of package survival
#'@param lim cox coef is computed on X1>lim
#' @param X1 array of time spent non smurf 
#' @param X2 array of time spent smurf  with same indexes as X1
#' @returns coex coefficient computed only on X1>lim
g_dep_sup <- function(lim,X1,X2){
  sc = coxph(Surv(X2[X1 > lim])~X1[X1>lim])$coef
  return(sc)
}

###########Parametric estimation of apparent death rate#######

#'Computes likelihood of gamma gompertz model with 3rd parameters explicitely 
#'computed where derivative of likelihood is 0
#'@param betagg beta parameter of model
#'@param bgg b parameter of model
#' @param X array of jumping times 
#' @returns array with log likelihood of gamma gompertz model on X with parameters betagg,
#' bgg and s where derivative of log likelihood is 0, and s 
lik_gg <- function(betagg,bgg,X){
  N = length(X)
  sgg = N/(sum(log(betagg-1 + exp(bgg*X)))-N*log(betagg) )
  return(c(N*log(bgg)+N*log(sgg)+bgg*sum(X)+N*sgg*log(betagg) - (sgg+1)*sum(log(betagg-1+exp(bgg*X))),sgg))
}

#'Computes derivative of likelihood of generized gamma model 
#'@param pga  parameter of model
#' @param X array of jumping times 
#' @param min min value for parameter l
#' @param max max value for parameter l 
#' @returns array with derivative of log likelihood, max lik parameter d , 
#' max lik parameter l
max_lik_ga <- function(pga,X,min,max){
  N = length(X)
  find_0_ga <- function(aga,X){
    dga = (pga/N)*sum((X/aga)^pga)
    return(-N*log(aga) +sum(log(X)) -  digamma(dga/pga)/pga)
  }
  aga = uniroot( function(aga)(find_0_ga(aga,X)) ,c(min,max))$root
  dga = (pga/N)*sum((X/aga)^pga)
  return(c(N/pga - sum(log(X/aga)*(X/aga)^pga) +digamma(dga/pga)*dga/pga^2,dga,aga))
}



###########Grid search fit with population data#########


#'Log likelihood of model 6 with parameters thet for death times in x
#'@param x array of death times 
#'@param thet an array of parameters (f,g,h,kd1e,kd2e,de,gam) for model 6
#'@returns the sum log likelihood of model 6 with parameters thet for death times in x
#'or -Inf if there was an error in computation (most likely integration) along the way
LL <- function(x, thet){
  N <- length(x)
  r <- 0
  for (i in 1:N){
    xi <- x[i]
    aux <- try(density_int(xi,thet),silent = TRUE)
    if(inherits(aux,'try-error')){
      r <- -Inf
      return(r)
    } 
    if(aux > 0){
      r <- r + log(aux)
    }
  }
  return(r)
}


#'Max likelihood grid search for parameters g and h in Model 6
#'@param thet parameter values for other parameters 
#'@param X array of death times 
#'@param minb min value for g
#'@param maxb max value for g 
#'@param minc min value for h
#'@param maxc max value for h 
#'@returns an array with the optimal values of g and h according to the log likelihood
grid_search <- function(thet,X,minb,maxb,minc,maxc){
  grid_b = minb + (0:10)/10*(maxb)
  grid_c = minc + (0:10)/10*(maxc)
  lb = length(grid_b)
  thet_res = thet
  lc = length(grid_c)
  res = matrix(0,lb,lc)
  for(i in 1:lb){
    thet_res[2] = grid_b[i]
    for (j in 1:lc){
      thet_res[3] <- grid_c[j]
      res[i,j] = LL(X, thet_res) #likelihood for set of parameters 
    }
  }
  ind_c = floor(which.max(res)/lb)+1*as.numeric( floor(which.max(res)/lb) != which.max(res)/lb)
  ind_b = which.max(res)-floor(which.max(res)/lb)*lb 
  ind_b = ind_b*(ind_b>0) + lb*(ind_b==0) #
  b_best = grid_b[ind_b]
  c_best = grid_c[ind_c]
  return(c(b_best,c_best))
}

#'Max likelihood grid search for parameter f in Model 6
#'@param thet parameter values for other parameters 
#'@param X array of death times 
#'@param mina min value for f
#'@param maxa max value for f 
#'@returnsthe optimal values of f according to the loglikelihood
grid_search_1p <- function(thet,X,mina,maxa){
  grid_a = mina + (0:10)/10*(maxa)
  la = length(grid_a)
  thet_res = thet
  res = matrix(0,la)
  for(i in 1:la){
    thet_res[1] = grid_a[i]
    res[i] = LL(X, thet_res)
  }
  ind_a = which.max(res)-floor(which.max(res)/la)*la
  ind_a = ind_a*(ind_a>0) + la*(ind_a==0)
  a_best = grid_a[ind_a]
  return(c(a_best))
}


##########Hazard rate and density functions#########

##Hazard rates 

#Smurf transition rates 
haz_smurf_gm = function(x,y=NULL,a,c,d)(a+c*exp(d*x)) #3 Gompertz-Makeham 
haz_smurf_g = function(x,y=NULL,c,d)(c*exp(d*x)) #2 Gompertz
haz_smurf_w = function(x,y=NULL,k,l)((k/l)*(x/l)^(k-1)) #3 Weibull
haz_smurf_p = function(x,y=NULL,a,b)(a*x^2+b) #2 Polynomial 

#Death rates once Smurf, in comments the number of parameters and corresponding model name
haz_c_i = function(x,y,k)(k*(x>=0)) #1 constant death rate model 1
haz_2c_i = function(x,y,k1,k2)(k1*(x<=24)+k2*(x>24)) #2 piecewise constant death rate model 2
haz_e_i = function(x,y,g1,g2,d)(g1+g2*exp(-d*x)) #3 decreasing exponential death rate model 3

haz_c_d = function(x,y,g,m,k)(exp(g*(y-m)) *k*(x>=0)) #2 constant death rate, dependence model 4
haz_2c_d = function(x,y,g,m,k1,k2)(exp(g*(y-m)) *(k1*(x<=24) + k2*(x>24))) #3 piecewise constant death rate,dependence model 5
haz_e_d = function(x,y,g,m,g1,g2,d)((exp(g*(y-m))) *(g1+g2*exp(-d*x))) #4 decreasing exponential death rate,dependence model 6

haz_c_2d = function(x,y,g,m,lim,k)((exp(g*(y-m))*(y>lim) + (y<=lim)) *k*(x>=0)) #3 constant death rate, piecewise dependence model 7
haz_2c_2d = function(x,y,g,m,lim,k1,k2)((exp(g*(y-m))*(y>lim) + (y<=lim)) *(k1*(x<=24) + k2*(x>24))) #4 piecewise constant death rate, piecewise dependence model 8
haz_e_2d = function(x,y,g,m,lim,g1,g2,d)((exp(g*(y-m))*(y>lim) + (y<=lim)) *(g1+g2*exp(-d*x))) #5 decreasing exponential death rate, piecewise dependence model 9

haz_dep_strat_c =#4 constant death rate, piecewise dependence model, piecewise base hazard 11
  function(x,y,g,m,lim,k1,k2)(exp(g*(y-m))*(y>lim)*(k1*(x>=0))+ (y<=lim)*( k2[1]*(x>=0) )) 
haz_dep_strat_2c =#6 piecewise constant death rate, piecewise dependence model, piecewise base hazard 11
  function(x,y,g,m,lim,k11,k21,k12,k22)((exp(g*(y-m))*(y>lim)*(k12*(x<=24) + k22*(x>24))+ (y<=lim)*(k11*(x<=24) + k21*(x>24)) )) 
haz_dep_strat_e = #8 decreasing exponential death rate, piecewise dependence model, piecewise base hazard 12
  function(x,y,g,m,lim,g11,g21,d1,g12,g22,d2)((exp(g*(y-m))*(y>lim)*(g12+g22*exp(-d2*x))+ (y<=lim)*(g11+g21*exp(-d1*x)) )) 



#Densities
#For time spent non-Smurf

df_ks <- function(t,a,c,d){
  return((a+ c*exp(d*t))* exp(-a*t - c/d*exp(d*t) + c/d))
} #Gompertz-Makeham
df_ks_g <- function(t,c,d){
  return((c*exp(d*t))* exp(- c/d*exp(d*t) + c/d))
}#Gompertz
df_ks_p <- function(t,a,b){
  return((b+a*t^2)*exp(- b*t - a*t^3/3))
}#Polynomial
df_ks_w <- function(t,l,k){
  return(dweibull(t,k,l))
}#Weibull


#For time spent Smurf

df_kd_i_c <- function(t,k){
  return(k*exp(-k*t))
} #Model 1
df_kd_i <- function(t,k1,k2){
  return((k1*(t<=24) + k2 *(t>24))*
           exp(-((k1*(24) + k2*(t-24))*(t>=24)
                 + k1*(t)*(t <24))))
} #Model 2
df_kd_i_e <- function(t,g1,g2,d){
  return((g1+g2*exp(-d*t))*exp(-(g1*t-g2/d*exp(-d*t) +g2/d)))
} #Model 3
df_kd <- function(t,u,g,m,k1,k2){
  return(exp(g*(u-m))* (k1*(t-u<=24) + k2*(t-u>24))*
           exp(-((k1*(24) + k2*(t-u-24))*(t-u>=24)
                 + k1*(t-u)*(t-u<24))*exp(g*(u-m)) ))
}#Model 5
df_kd_dep_e <- function(t,u,g,m,g1,g2,d){
  return(exp(g*(u-m))* (g1+g2*exp(-d*(t-u)))*
           exp(-(g1*(t-u)-g2/d*exp(-d*(t-u)) +g2/d)*exp(g*(u-m)) ))
} #Model 6
df_kd_2dep <- function(t,u,g,m,lim,k1,k2){
  return((exp(g*(u-m))*(u>=lim) + 1*(u<lim) )* (k1*((t-u)<=24) + k2 *((t-u)>24))*
           exp(-((k1*(24) + k2*(t-u-24))*((t-u)>=24)
                 + k1*(t-u)*((t-u)<24))*(exp(g*(u-m))*(u>=lim )+ 1*(u<lim ))))
} #Model 8
df_kd_2dep_cont <- function(t,u,g,m,lim,g1,g2,d){
  return((t-u>=0)*((exp(g*(u-m))*(u>=lim) + 1*(u<lim) )* (g1+g2*exp(-d*(t-u)))*
                     exp(-(g1*(t-u)-g2/d*exp(-d*(t-u)) + g2/d)*(exp(g*(u-m))*(u>=lim) + 1*(u<lim) ))))
} #Model 9
df_kd_2dep_strat1 <- function(t,u,g,m,lim,g11,g21,d1){
  return((t-u>=0)*( ( g11+g21*exp(-d1*(t-u))   )*
                      exp(-( (g11*(t-u)-g21/d1*exp(-d1*(t-u))   +g21/d1)))))
} #Model 12 strata 1
df_kd_2dep_strat <- function(t,u,g,m,lim,g11,g21,d1,g12,g22,d2){
  return((t-u>=0)*((exp(g*(u-m))*(u>=lim) + 1*(u<lim) )* ((u>=lim)*( g12+g22*exp(-d2*(t-u))) +(u<lim)* (g11+g21*exp(-d1*(t-u)) )  )*
                     exp(-((u>=lim)* (g12*(t-u)-g22/d2*exp(-d2*(t-u)) +g22/d2) +(u<lim)* (g11*(t-u)-g21/d1*exp(-d1*(t-u))   +g21/d1))*(exp(g*(u-m))*(u>=lim) + 1*(u<lim) ))))
} #Model 12


#Marginal density of time spent Smurf in models with dependence, defined as a convolution
time_smurf_dens <- function(u,dens){
  int = integrate(function(x)(dens(u,x)),0,1000,rel.tol=10^(-4))$value
  return(int)
} #auxiliary function to compute the integral in the convolution

marg_2strat <- function(u,t,au,cu,du,g1_s1,g2_s1,d_deps1,g1_s2,g2_s2,d_deps2,g_sup,m2,lim){
  return(df_ks(t,au,cu,du)*df_kd_2dep_strat(u+t,t,g_sup,m2,lim,g1_s1,g2_s1,d_deps1,g1_s2,g2_s2,d_deps2))} #Model 12
marg_2dep <- function(u,t,au,cu,du,g1_d2,g2_d2,d_dep2,g_sup,m2,lim){
  return(df_ks(t,au,cu,du)*df_kd_2dep_cont(u+t,t,g_sup,m2,lim,g1_d2,g2_d2,d_dep2))
}#Model 9
marg_dep <- function(u,t,au,cu,du,g1_d,g2_d,d_dep,g,m)(df_ks(t,au,cu,du)*df_kd_dep_e(u+t,t,g,m,g1_d,g2_d,d_dep)) #Model 6


#Densities of lifetimes
df_td <- function(t,g,m,k1,k2,a,c,d){
  return(integrate(function(u) (df_ks(u,a,c,d)*df_kd(t,u,g,m,k1,k2)), 0,t,rel.tol = 0.000001)$value)
} #Model 2
df_td_i <- function(t,g1,g2,d,f,g,h){
  return(integrate(function(u) (df_ks(u,f,g,h)*df_kd_i_e(t-u,g1,g2,d)), 0,t,rel.tol = 0.000001)$value)
} #Model 3
df_td_dep_cont <- function(t,g,m,g1,g2,d,a,b,c){
  return(integrate(function(u) (df_ks(u,a,b,c)*df_kd_dep_e(t,u,g,m,g1,g2,d)), 0,t,rel.tol = 0.000001)$value)
} #Model 6
df_td_2dep <- function(t,g,m,lim,k1,k2,a,c,d){
  return(integrate(function(u) (df_ks(u,a,c,d)*df_kd_2dep(t,u,g,m,lim,k1,k2)), 0,t,rel.tol = 0.000001)$value)
} #Model 8
df_td_2dep_cont <- function(t,g,m,lim,g1,g2,d,a,c,b){
  return(integrate(function(u) (df_ks(u,a,c,b)*df_kd_2dep_cont(t,u,g,m,lim,g1,g2,d)), 0,t,rel.tol = 0.000001)$value)
} #Model 9
df_td_2dep_strat <- function(t,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d){
  return(integrate(function(u) (df_ks(u,a,c,d)*df_kd_2dep_strat(t,u,g,m,lim,g11,g21,d1,g12,g22,d2)), 0,t,rel.tol = 0.000001)$value)
} #Model 12

#' Density of death time for Gompertz-Makeham one phase model ad+cd*exp(dd.t)
#' @param t time at which to estimate 
#' @param ad param of gm
#' @param cd param of gm
#' @param dd param of gm
#' @returns cdf of gm param at time t
dis_makeham <- function(t,ad,cd,dd){
  return((ad+cd*exp(dd*(t) )) * exp(-(ad*(t) + cd/dd*exp(dd*t)-cd/dd)))
}


#'Auxiliary function, returns the joint density of time smurf v and time 
#'non smurf t-v
#'@param v time spent smurf 
#'@param t time of death 
#'@param thet an array of parameters (f,g,h,kd1e,kd2e,de,gam) for model 6
#'@returns joint density of time smurf v and time non smurf t-v 
density_EM <- function(v,t,thet){
  a = thet[1]
  b = thet[2]
  c = thet[3]
  g1 = thet[4]
  g2 = thet[5]
  d = thet[6]
  g = thet[7]
  haz1 = Vectorize(function(x)((g1+g2*exp(-d*x))*exp(g*(t-x))))
  int_haz1 = Vectorize(function(x)(g1*x*exp(g*(t-x)) -g2/d*exp(g*(t-x))*exp(-d*x)+g2/d*exp(g*(t-x))))
  dens = function(x)(haz1(x)*exp(-int_haz1(x)))
  dens2= function(x)((a + b*exp(c*x))*exp(-a*x-(b/c)*exp(c*x) +b/c))
  return(dens(v)*dens2(t-v))
}


#'Density of death at time t for model 6 with parameters in thet
#'@param t time of death 
#'@param thet an array of parameters (f,g,h,kd1e,kd2e,de,gam) for model 6
#'@returns density of death at t
density_int <- function(t,thet){
  dens <- integrate(function(v)(density_EM(v,t,thet)),0.0001,t,subdivisions = 100L, rel.tol = 1e-13, 
                    abs.tol = 1e-13)$value
  return(dens)
}


#Cumulative distribution functions of lifetime 
F_2dep_strat <-function(t,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d){
  f_int = Vectorize(function(u)(df_td_2dep_strat(u,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d)))
  return(integrate(f_int,0,t)$value)
} #Model 12
F_dep_cont <-function(t,g,m,g1,g2,d1,a,c,d){
  f_int = Vectorize(function(u)(df_td_dep_cont(u,g,m,g1,g2,d1,a,c,d)))
  return(integrate(f_int,0,t,subdivisions = 3000)$value)
} #Model 6
F_2dep_cont <-function(t,g,m,lim,g11,g21,d1,a,c,d){
  f_int = Vectorize(function(u)(df_td_2dep_cont(u,g,m,lim,g11,g21,d1,a,c,d)))
  return(integrate(f_int,0,t,subdivisions = 3000)$value)
} #Model 9

#Apparent death hazard rates 
haz_death <- function(t,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d){
  return(df_td_2dep_strat(t,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d)/(1-F_2dep_strat(t,g,m,lim,g11,g21,d1,g12,g22,d2,a,c,d)))
} #Model 12
haz_death6 <- function(t,g,m,g1,g2,d1,a,c,d){
  return(df_td_dep_cont(t,g,m,g1,g2,d1,a,c,d)/(1-F_dep_cont(t,g,m,g1,g2,d1,a,c,d)))
} #Model 6
haz_death9 <- function(t,g,m,lim,g11,g21,d1,a,c,d){
  return(df_td_2dep_cont(t,g,m,lim,g11,g21,d1,a,c,d)/(1-F_2dep_cont(t,g,m,lim,g11,g21,d1,a,c,d)))
} #Model 9


##########Mice data treatment and estimation############

#' Transforms a data frame into two arrays with time spent non-smurf and smurf
#' for a population of mice
#' @param mice a dataframe
#' @returns a list with elements S and NS corresponding to times spent smurf and 
#' non-smurf 
mice_treatment <- function(mice){
  Non_Sm = 1:max(mice['mice'])
  Sm = 1:max(mice['mice'])
  #B = mice['Smurf_state' == 'birth']
  j=1 #counting the number of mice for which data is complete
  for (i in 1:max(mice['mice'])){ #iterating over every mouse index
    tab = mice |> subset(mice==i)
    life = tab['REMAINTIME']
    smurf = tab |> subset(Smurf_state == 'S')
    k = 1
    while( k < dim(smurf)[1] && smurf$X[k+1] - smurf$X[k] >1) {
      k = k+1 #moment where mouse goes from NS to S
    }
    if (dim(smurf)[1]>0 ){ #check that it has been scored S 
      t_s = -smurf$REMAINTIME[k]
      t_ns = -min(life)-t_s
      Non_Sm[j] = t_ns #add data to table
      Sm[j] = t_s
      j = j+1}
  }
  Non_Sm = Non_Sm[1:(j-1)] #only keep the indexes that have been filled
  Sm = Sm[1:(j-1)]
  return(list('S' = Sm, 'NS' = Non_Sm ))
}

#' Transforms a data frame into two arrays with time spent non-smurf and smurf
#' for a population of mice, uniformly assigning time spent smurf in the last
#' observation interval if the mouse was not scores smurf
#' @param mice a dataframe
#' @returns a list with elements S and NS corresponding to times spent smurf and 
#' non-smurf 
mice_treatment_add_smurf <- function(mice){
  Non_Sm = 1:max(mice['mice'])
  Sm = 1:max(mice['mice'])
  mouse = unique(mice['mice']) #mouse indexes
  j=1
  for (m in 1:dim(mouse)[1]){#iterate over mouse indexes 
    i = mouse$mice[m]
    tab = mice |> subset(mice==i)
    life = tab$REMAINTIME
    smurf = tab |> subset(Smurf_state == 'S')
    k = 1
    while( k < dim(smurf)[1] && smurf$X[k+1] - smurf$X[k] >1) {
      k = k+1 #moment where mouse goes from NS to S
    }
    if (dim(smurf)[1]>0 ){#check that it has been scored S 
      t_s = -smurf$REMAINTIME[k]
      t_ns = -min(life)-t_s
      Non_Sm[j] = t_ns#add data to table
      Sm[j] = t_s
      j = j+1}
    else {#if not scored S, time spent smurf is randomly picked in the last 
      #observation interval
      k = length(life)
      t_s = runif(1,-life[k] ,-life[k-1])
      t_ns = -min(life)-t_s
      Non_Sm[j] = t_ns #add data to table
      Sm[j] = t_s
      j = j+1}
  }
  Non_Sm = Non_Sm[1:(j-1)] #only keep the indexes that have been filled
  Sm = Sm[1:(j-1)]
  return(list('S' = Sm, 'NS' = Non_Sm ))
}


#max likelihood of ks = f + g*exp(h*t), other method
#'Find max likelihood estimations of parameters f and g and value of log likelihood 
#'of gompertz makeham smurf hazard rate f + g.exp(h.t) for a given h 
#' @param h value of h
#' @param X1 an array of time spent non-smurf
#' @returns an array with log likelihood and optimals parameters f, g of 
#' gompertz makeham smurf hazard rate for this value of k
param_Gomp_Mak2 <- function(X1,h){
  N = length(X1)
  X1ord = sort(append(0,X1,0))#reorder times
  NS = (N-1):0
  der_vrais3 <- function(h){
    int = sum((exp(h*X1ord[2:(N+1)]) - exp(h*X1ord[1:N]))*NS)/h
    xm = min(X1ord[2:N])
    find_0 <- function(j){#function to determine parameter g (=variable j)
      f = (N-j*int)/sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)
      sum1 = sum(1/(f+j*exp(h*X1))) - sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)
      #return(sum(log(f+j*exp(d*X1ord))) - sum(f*X1ord+j/d*exp(d*X1ord) - j/d)) 
      #return log-likelihood to check that it is max, commented 
      return(sum1)#derivative at j 
    }
    clim = N/(int - sum((X1ord[2:(N+1)] - X1ord[1:N])*NS)*exp(h*xm)) #to determine where to look for c2
    int2 = sum((X1ord[2:(N+1)]*exp(h*X1ord[2:(N+1)]) - X1ord[1:N]*exp(h*X1ord[1:N]))*NS)/h
    inttot = int2 - int/h
    A_test = (1:99)*(clim - clim/10000)/(99.5) +clim/10000
    apply_test = sapply(A_test, find_0)
    i_max =  max(which(apply_test < 0))
    g = try(uniroot(find_0,c(A_test[i_max],A_test[i_max+1]),tol=0.0000001)$root,silent=TRUE) #c2 where the derivative is 0
    #g = try(uniroot(find_0,c(ming,maxg),tol=0.0000001)$root,silent=TRUE)
    if(inherits(g,'try-error')){
      return(NA)
    } 
    #g = uniroot(find_0,c(clim/10,clim/1.00000001),tol=0.0000001)$root #g where the derivative is 0
    #g = uniroot(find_0,c(ming,maxg),tol=0.0000001)$root
    f = (N-g*int)/sum((X1ord[2:(N+1)] - X1ord[1:N])*NS) #f computed directly with analytic result
    sum1 = sum(X1*exp(h*X1)/(f+g*exp(h*X1))) - inttot
    return(c(f,g))
  }
  der= der_vrais3(h)
  f = der[1]
  g = der[2]
  return(c(log_lik_gm(c(f,g,h),X1),f,g ))
}#min and max of h (parameter inside the exponential)

##########Miscellaneous data treatment##########
#' Number of smurf individuals at time t for experimental data
#' @param t time at which to evaluate number of non smurf individuals
#' @param X1 array of time spent non smurf
#' @param X2 array of time spent smurf, with same indexes as X1 
#' @returns the number of smurf individuals at time t
Smurf_1_2 <- function(t,X1,X2){#counting process of smurf individuals for experimental data 
  #X1 if smurf jumping times and X2 death times
  s = sum(X1 <= t & X2+X1 > t) #already smurf, not dead yet 
  return(s)
}

#'Survival function of population with model 6 and parameter thet at time t
#'@param t time at which to estimate survival function
#'@param thet an array of parameters (f,g,h,kd1,kd2,d,gam) for model 6
#'@returns value of survival function at t 
survival_tot <- function(t,thet){
  fun_aux <- Vectorize(function(x)(density_int(x,thet)))
  aux<- integrate(fun_aux,0,t)$value
  return(1-aux)
}

#' Theoretical number of Smurf with model 6
#' @param t time at which to evaluate
#' @param thet an array of the 7 parameters for model 6 in the order
#' f g h kd1 kd2 d gam
#' @returns the theoretical number of smurf at time t
Smurf_th_1dep <- function(t,thet) { #theoretical smurf population with piece wise base hazard and conditional dependance
  a = thet[1]
  b = thet[2]
  c = thet[3]
  g1 = thet[4]
  g2 = thet[5]
  d = thet[6]
  gam = thet[7]
  integr <- function(u){
    (a+b*exp(c*u))*exp(-(g1*(t-u)-g2/d*exp(-d*(t-u))+g2/d)*exp(gam*(u))
                       - a*u-b/c*exp(c*u) +b/c)
  }
  return(integrate(integr,0,t)$value)
}

#'Transforms survival function into death times
#'@param survival function with a value for each day
#'@returns an array of death times (censored by 24h)
transf_surv <- function(DGRP){
  Times_DGRP = (0:length(DGRP))*24
  ND = -DGRP[2:length(DGRP)]+DGRP[1:(length(DGRP)-1)]
  X_DGRP = c()
  for (i in 1:(length(DGRP)-1)){
    X_DGRP = c(X_DGRP, rep(i*24, ND[i]))
  }
  return(X_DGRP)
}