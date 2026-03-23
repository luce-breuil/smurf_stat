#This file contains the code to simulate fly populations for all the models 


#' Returns a portion of code and the parameters in the correct form to simulate
#'a population for a given parametrization of smurf hazard rates
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param param an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#'@returns a named list with the object 'code', a string to input in the code by 
#' sim_pop_n and object 'par', a named list of parameters useable to simulate 
#' the population
get_code_ks = function(p,param){
  if (p=='p'){code = 'a*(I.age(t))*(I.age(t)) + b'
  par = list('a' = param[1], 'b' = param[2])
  }#polynomial
  else if (p=='g'){ code = 'c*exp(I.age(t)*d)'
  par = list('c' = param[1], 'd' = param[2])
  }#gompertz
  else if (p=='gm'){code ='a2 +c2* exp(I.age(t)*d2)'
  par = list('a2' = param[1], 'c2' = param[2], 'd2' = param[3])
  }#gompertz-makeham
  else if (p=='w'){ code='kw/lw*exp((kw-1)*log(I.age(t)/lw))'
  par = list('kw'= param[1], 'lw' = param[2])
  }#weibull
  else if (p=='np'){ code = 'gam_haz_simul_s(I.age(t))'
  par = list('gam_haz_simul_s' = param)
  }#non-parametric
  return(list('code'=code,'par'=par))
}

#' Returns a portion of code and the parameters in the correct form to simulate
#'a population for a given parametrization of death hazard rates
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param param an array with the parameters for the chosen death hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#'@returns a named list with the object 'code', a string to input in the code by 
#' sim_pop_n and object 'par', a named list of parameters useable to simulate 
#' the population
get_code_kd = function(i,param){
  if (i==1){ code = 'if (I.smurf) 
                                    {result = kDe;}
                                           else 
                                               {result = 0;}'
  par = list('kDe' = param[1])
  }#constant
  else if (i==2) {
    code = ' if (I.smurf && t-I.time_smurf <= 24 )
                                {result =kd1i;}
                                else if (I.smurf)
                                {result =kd2i; }
                              else 
                                               {result = 0;}'
    par = list('kd1i' = param[1], 'kd2i' = param[2])
  }#2 constants
  else if (i==3) {
    code = ' if (I.smurf)
                                {result =g1+g2*exp(-d_kd*(t-I.time_smurf)); }
                              else 
                                               {result = 0;}'
    par = list('g1' = param[1], 'g2' = param[2], 'd_kd' = param[3])
  }#exponential decreasing
  else if (i == 4){
    code = 'if (I.smurf )
                              {result = kdp*exp(g*(I.time_smurf-moy));}
              
                              else
                                               {result = 0;}'
    par = list('kdp'= param[1], 'g' = param[2], 'moy' = param[3])
  }#constant Cox dependence
  else if (i==5) {
    code= 'if (I.smurf && t-I.time_smurf  <= 24 )
                            {result = kd1*exp(g*(I.time_smurf-moy));}
                              else if (I.smurf)
                              {result = kd2*exp(g*(I.time_smurf-moy));}
                              else
                                               {result = 0;}'
    par = list('kd1' = param[1], 'kd2' = param[2], 'g' = param[3],'moy'= param[4])
  }#2 constants Cox dependence
  
  else if(i==6) {
    code = 'if (I.smurf )
                              {result = (g1_d+g2_d*exp(-d_dep*(t-I.time_smurf)))*exp(g*(I.time_smurf-moy));}
                              else
                                               {result = 0;}'
    par = list('g1_d' = param[1], 'g2_d' = param[2], 'd_dep' = param[3], 'g' = param[4], 'moy' = param[5])
  }#exponential decreasing Cox dependence
  else if (i==7) {
    code ='if (I.smurf && I.time_smurf >= lim)
                                {result = k_2pc*exp(g*(I.time_smurf-moy));}
                            else if (I.smurf)
                            {result =  k_2pc;}
                              else
                                               {result = 0;}'
    par = list('k_2pc' = param[1],'g' = param[2], 'moy' = param[3],'lim' = param[4] )
  }#constant piecewise Cox dependence
  else if (i==8) {
    code = 'if (I.smurf && t-I.time_smurf  <= 24 &&I.time_smurf >= lim)
                                {result = exp(g_sup*(I.time_smurf-moy))*kd1;}
                                else if (I.smurf && t-I.time_smurf  <= 24 &&I.time_smurf < lim)
                                {result =kd1;}
                            else if (I.smurf && I.time_smurf < lim)
                                {result = kd2 ;}
                            else if (I.smurf && I.time_smurf >= lim)
                            {result = kd2*exp(g_sup*(I.time_smurf-moy)) ;}
                              else
                                               {result = 0;}'
    
    par = list('kd1' = param[1], 'kd2' = param[2], 'g_sup' = param[3],'moy'= param[4],'lim' = param[5])
  }#2 constants piecewise Cox dependence
  else if (i==9) {
    code =  'if (I.smurf && I.time_smurf >= lim)
                                {result = (g1_d2+g2_d2*exp(-d_dep2*(t-I.time_smurf)))*exp(g*(I.time_smurf-moy));}
                            else if (I.smurf)
                            {result =  (g1_d2+g2_d2*exp(-d_dep2*(t-I.time_smurf)));}
                              else
                                               {result = 0;}'
    par = list('g1_d2' = param[1], 'g2_d2' = param[2],'d_dep2' = param[3],'g' = param[4], 'moy' = param[5],'lim'=param[6])
  }#exponential decreasing piece wise Cox dependence
  
  else if(i==10) {
    code = 'if (I.smurf && I.time_smurf <lim)
                              {result = (kd_s1);}
                                else if (I.smurf && I.time_smurf >= lim)
                                {result = (kd_s2)*exp(g*(I.time_smurf - moy));}

                              else
                                               {result = 0;}'
    par = list('kd_s1' = param[1], 'kd_s2' = param[2], 'g' = param[3], 'moy' = param[4], 'lim' = param[5])
  }#constant Cox dependence, different base hazards
  else if (i==11) {
    code = 'if (I.smurf && t-I.time_smurf  <= 24 &&I.time_smurf >= lim)
                                {result = exp(g*(I.time_smurf-moy))*kd1t;}
                                else if (I.smurf && t-I.time_smurf  <= 24 &&I.time_smurf < lim)
                                {result =kd1ti;}
                            else if (I.smurf && I.time_smurf < lim)
                                {result = kd2ti ;}
                            else if (I.smurf && I.time_smurf >= lim)
                            {result = kd2t*exp(g*(I.time_smurf-moy)) ;}
                              else
                                               {result = 0;}'
    par = list('kd1ti' = param[1], 'kd2ti' = param[2],'kd1t' = param[3], 'kd2t' = param[4],
               'g' = param[5],'moy'= param[6],'lim' = param[7])
  }#2 constants Cox dependence, different base hazards
  else if(i==12) {
    code = 'if (I.smurf && I.time_smurf <lim)
                              {result = (g1_s1+g2_s1*exp(-d_deps1*(t-I.time_smurf)));}
                                else if (I.smurf && I.time_smurf >= lim)
                                {result = (g1_s2+g2_s2*exp(-d_deps2*(t-I.time_smurf)))*exp(g*(I.time_smurf - moy));}

                              else
                                               {result = 0;}'
    par = list('g1_s1' = param[1], 'g2_s1' = param[2], 'd_deps1' = param[3],
               'g1_s2' = param[4], 'g2_s2' = param[5], 'd_deps2' = param[6],
               'g' = param[7], 'moy' = param[8], 'lim' = param[9])
  }#exponential decreasing Cox dependence, different base hazards
  else if(i==13) {
    code = 'if (I.smurf && I.time_smurf <lim)
                              {result = (g1_s1+g2_s1*exp(-d_deps1*(t-I.time_smurf)));}
                                else if (I.smurf && I.time_smurf >= lim)
                                {result = k2*exp(g*(I.time_smurf - moy));}

                              else
                                               {result = 0;}'
    par = list('g1_s1' = param[1], 'g2_s1' = param[2], 'd_deps1' = param[3],
               'k2' = param[4], 'g' = param[5], 'moy' = param[6], 'lim' = param[7])
  }# Cox dependence, different base hazards (constant and exponential)
  else if(i=='np') {
    code = 'if (I.smurf )
                              {result = gam_haz_simul(t-I.time_smurf);}
              
                              else
                                               {result = 0;}'
    par = list('gam_haz_simul' = param)
  }#non-parametric
  
  
  return(list('code'=code,'par'=par))
}

#' Returns a simulated population for a given parametrization of smurf and death hazard rates
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#' @param init_size population size, default 1159
#' @param N population size for simulation, default 1159
#' @param n time of simulation is n*10, default is 150
#' @returns an object of type population
sim_pop_n <- function(p,i,paramS, paramD,init_size = 1159, n = 150){
  pop_init <- IBMPopSim::population(data.frame( birth = rep(0, init_size), 
                                                death = NA, smurf = rep(FALSE,init_size), 
                                                time_smurf = rep(0, init_size)),id = TRUE)
  
  smurf_param = get_code_ks(p,paramS) #smurf rate parameters and code
  code_smurf = smurf_param$code
  paramS1 = smurf_param$par
  smurf = IBMPopSim::mk_event_individual(type = 'swap', 
                                         intensity_code = paste('if (I.smurf) {result = 0;} 
                                      else {result = ',code_smurf,';}'),
                                         kernel_code = 'I.smurf = true; I.time_smurf = t;'
  ) #defining smurf event
  death_param = get_code_kd(i,paramD) #death rate parameters and code
  code_death = death_param$code
  paramD1 = death_param$par
  params = c(paramS1,paramD1)
  death =IBMPopSim::mk_event_individual(type = 'death', 
                                        #kernel_code = 'I.death = t;',
                                        intensity_code = code_death
  )#defining death event
  birth_death <-IBMPopSim::mk_model(characteristics = IBMPopSim::get_characteristics(pop_init),
                                    events = list(death,smurf),
                                    parameters = params)
  sim_out <- IBMPopSim::popsim(model = birth_death,
                               initial_population = pop_init,
                               events_bounds = c('death' = 1, 'swap' = 1),
                               parameters = params,
                               time = (0:((n-1)))*10) #simulation
  return(sim_out$population)
}



#' Plots the pointwise 95% CI of non-smurf population for a given parametrization
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#' @param th experimental times spent non-smurf
#' @param init_size population size, default 1159
#' @param N population size for simulation, default 1159
#' @param n_rep number of simulations
#' @param array, array of times for simulation, default is (0:149)*10
#' @returns nothing, but plots the experimental non-smurf survival function and 
#' empirical 95 CI
NS_CI <- function(p,i, paramS, paramD,th,init_size =1159,nrep=100, array = NA){
  #initializing population 
  if(any(is.na(array))){array = (0:149)*10}
  n = length(array)
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
  
  M = matrix(0,nrow= (n-1),ncol=nrep)
  for (j in 1:nrep){#simulating nrep populations
    sim_out <- popsim(model = birth_death,
                      initial_population = pop_init,
                      events_bounds = c('death' = 1, 'swap' = 1),
                      parameters = params,
                      time = array) #simulation
    pop_out1 <- sim_out$population
    nsmurf_num1 <- lapply(1:(n-1), function (i) (
      (sum(pop_out1[[i]]$smurf == FALSE))))
    for (k in 1:(n-1)){
      M[k,j] =nsmurf_num1[[k]]}
  }
  CI = sapply(1:(n-1) , function(g)(quantile(M[g,],c(0.025,0.975) ) ))
  return(CI)
}



#' Plots the pointwise 95% CI of smurf survival population for a given parametrization
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#' @param th experimental time spent smurf
#' @param init_size population size, default 1159
#' @param N population size for simulation, default 1159
#' @param n_rep number of simulations
#' @param array, array of times for simulation, default is (0:149)*10
#' @returns nothing, but plots the experimental survival function and the 95 CI
S_CI <- function(p,i, paramS, paramD,th,init_size =1159, nrep = 100, array = NA){
  #initializing population 
  if(any(is.na(array))){array = (0:149)*10}
  n = length(array)
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
  
  M = matrix(0,nrow= (init_size),ncol=nrep)
  for (j in 1:nrep){#simulating nrep populations
    sim_out <- popsim(model = birth_death,
                      initial_population = pop_init,
                      events_bounds = c('death' = 1, 'swap' = 1),
                      parameters = params,
                      time = array) #simulation
    pop_out1 <- sim_out$population
    times_smurf = (pop_out1[[n-1]]$death -  pop_out1[[n-1]]$time_smurf)
    t = sort(times_smurf)
    t = t*(!is.na(t))
    for (k in 1:init_size){
      M[k,j] =t[k]}
  }
  
  CI = sapply(1:(init_size) , function(g)(quantile(M[g,],c(0.025,0.975) ,na.rm=TRUE) ))
  return(CI)
  
  
}



#' Plots the pointwise 95% CI of smurf population for a given parametrization
#' @param p to indicate which smurf hazard rate parametrization to use, either np, g, gm, p or w 
#' (for non parametric, gompertz, gompertz makeham, polynomial or weibull)
#' @param i o indicate which death hazard rate parametrization to use, either 
#' np or a number from 1 to 12 depending on the model (see report)
#' @param paramS an array with the parameters for the chosen smurf hazard rate
#'  parametrization OR an object of type linfun of IBMPopSim if p=np
#' @param paramD an array with the parameters for the chosen death rate
#'  parametrization OR an object of type linfun of IBMPopSim if i=np
#' @param th experimental smurf population on data evaluated at each hour
#' @param init_size population size, default 1159
#' @param N population size for simulation, default 1159
#' @param n time of simulation is n*10, default is 150
#' @param n_rep number of simulations
#' @returns nothing, but plots the experimental pop and the 95 CI
smurf_CI <- function(p,i, paramS, paramD,th,title=NA, init_size =1159, n = 150,nrep=100){
  M = matrix(0,nrow= (n-1),ncol=nrep)
  for (j in 1:nrep){ #simulating nrep populations
    pop_out1 <- sim_pop_n(p,i,paramS,paramD,init_size, n)
    smurf_num1 <- lapply(1:(n-1), function (g) (
      (sum(pop_out1[[g]]$smurf == TRUE) - sum(!is.na(pop_out1[[g]]$death)))))
    for (k in 1:(n-1)){
      M[k,j] = smurf_num1[[k]]}
  }
  CI = sapply(1:(n-1) , function(g)(quantile(M[g,],c(0.025,0.975) ) ))
  return(CI)
}
