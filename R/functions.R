###############################################################################
#' Set-up age groups and find ageing rate
#' Takes the global variables numMonthly and final_age, along
#' With the appropriate greater Perth data file. Assuming a non-uniform,
#' but constant population distribution
#' @param year ABS greater Perth population year
#' @param area 0 or 1. Default is 0 = "Greater Perth", 1 = "Southern WA" 
#' @return list nAges, pop_groups, age_in, age_out
#'
age_structure <- function(year, area = 0)
{
  age_vect_years <- c(seq(0, numMonthly/12, 1/12), seq(10, final_age-5, 5)) #age lower limit in years
  age_vect_months <- age_vect_years*12 #age lower limit in months
  size_cohorts_months <- c(diff(age_vect_months), final_age*12 - age_vect_months[length(age_vect_months)])
  nAges <- length(age_vect_years) #number of age groups
  
  if(!area)
  {  
    pop_data <- get(paste("data.greaterPerth.", year, sep = ""))
  }else{
    pop_data <- get(paste("data.southernWA.", year, sep = ""))
  }  
  
  pop_groups <- c(rep(as.numeric(pop_data$population[1])/numMonthly, numMonthly),
                  as.numeric(pop_data$population[2:(nAges - numMonthly +1)]))
  #uniformly distributed pop from first 5-year group into monthly groups
  birth_pop <- pop_groups[1]
  
  age_rate <- (1/pop_groups)*birth_pop #to maintain the same non-uniform age distribution
  
  age_in <- c(birth_pop, age_rate[-nAges])
  age_out <- c(age_rate[-nAges], birth_pop/pop_groups[nAges])
  
  data.pop.model <- pop_data[-1, ]
  
  if(!area)
  {  
    data.pop.model <- data.pop.model %>% tibble::add_row(lga = "Greater Perth",
                                                         lower.age.limit = age_vect_years[which(age_vect_years < 5)],
                                                         year = year,
                                                         population = pop_groups[1:numMonthly],
                                                         .before = 1)
  }else{
    data.pop.model <- data.pop.model %>% tibble::add_row(lga = "Southern WA",
                                                         lower.age.limit = age_vect_years[which(age_vect_years < 5)],
                                                         year = year,
                                                         population = pop_groups[1:numMonthly],
                                                         .before = 1)
  }  
  data.pop.model <- data.pop.model[which(data.pop.model$lower.age.limit < final_age),]
  save(data.pop.model, file = "data.pop.model.rda" )
  
  return(list("nAges" = nAges, "pop_groups" = pop_groups,
              "age_years" = age_vect_years, "age_in" = age_in, "age_out" = age_out))
}

###############################################################################
#' Get contact matrix
#' Using R package conmat, population for greater Perth and assumed
#' age groupings for model
#'
#' @param year year of ABS population data data
#' @param sym return a symmetric matrix
#' @param area 0 or 1. Default is 0 = "Greater Perth", 1 = "Southern WA"
#' @param save 1 to save matrix, 0 to return
#' @param fileName if save, give file name
#' @return matrix
#'
get_contact_matrix <- function(year, sym = 0, area = 0, save = 1, fileName = NA)
{
  nAges <- numMonthly + (final_age - numMonthly/12)/5 #5 - year age groups
  
  if(!area)
  {  
    pop_data <- get(paste("data.greaterPerth.", year, sep = ""))
  }else{
    pop_data <- get(paste("data.southernWA.", year, sep = ""))
  }  
  pop <- pop_data[which(pop_data$lower.age.limit < final_age),]
  
  setting_models <- conmat::fit_setting_contacts(
    contact_data_list = conmat::get_polymod_setting_data(countries = "United Kingdom"),
    population = conmat::get_polymod_population(countries = "United Kingdom")
  )
  
  raw_mat <- conmat::predict_setting_contacts(
    population = pop,
    contact_model = setting_models,
    age_breaks = seq(0, final_age, by = 5)
  )
  raw_mat <- raw_mat$all
  
  mixing <- matrix(0, ncol = nAges, nrow = nAges)
  rownames(mixing) <- c(rep(rownames(raw_mat)[1], numMonthly),
                        rownames(raw_mat)[2:length(rownames(raw_mat))])
  colnames(mixing) <- rownames(mixing)
  
  #make symmetric
  if(sym == 1)
  {
    temp <- matrix(0, ncol = ncol(raw_mat), nrow = nrow(raw_mat))
    
    for (i in 1:nrow(raw_mat)) {
      for (j in 1:ncol(raw_mat)){
        temp[i,j] <- (raw_mat[i,j] + raw_mat[j,i]) / 2
      }
    }
    raw_mat <- temp
  }
  
  #lower right section
  mixing[(numMonthly+1):nAges, (numMonthly+1):nAges] <- raw_mat[2:nrow(raw_mat),2:nrow(raw_mat)] #age groups 5-9y to 75-79y
  
  #upper left section
  mixing[1:numMonthly,1:numMonthly] <- raw_mat[1,1]/numMonthly #monthly age groups from 0-4y
  
  #lower left section
  for (i in 1:numMonthly){
    mixing[(numMonthly+1):nAges,i] <- raw_mat[2:nrow(raw_mat),1]
  }
  
  #upper right section
  for (i in (numMonthly+1):nAges){
    mixing[1:numMonthly,i] <- raw_mat[(i-(numMonthly-1)),1]/numMonthly
  }
  
  mixing <- t(mixing) #change to rows "To" and cols "From"
  
  if(save)
  {
    write.csv(mixing, paste(fileName, ".csv", sep = ''))
    return(1)
  }
  
  return(mixing)
}
############################################################################
#' From model output, find the time range that matches to the observed
#' data set. Do this by assuming a "peak" month, default is 7 (July)
#'
#' @param counts list of modelled detected infections in aggregated age groups,
#'               in same form as observational data
#' @param month month 1 - 12, Default is 7 (July)
#' @param obs table of monthly hospitalisation data for 4 age groups, at this
#'            stage assuming complete calendar years
#' @param biennial if 0 then observed series is annual, if 1 then observed series 
#'                 is biennial with higher peak in the last year, if 2 then
#'                 observed series is biennial with lower peak in the last year          
#' @return list(index_start, index_end, phase)
#'
findIndexRange <- function(counts, month = 7, obs, biennial = 0)
{
  nYears <- nrow(obs)/12

  aggDetInc <- apply(counts, 1, sum)
  peaks <- pracma::findpeaks(aggDetInc)

  
  if(biennial == 0) #annual series
  {
    p_index <- peaks[(nrow(peaks)-1),2] # find index of second-last peak
    index_start <- p_index - (nYears - 1)*12 - (month - 1)
    index_end <- p_index + 12 - month
    return(list("start" = index_start, "end" = index_end))
  }  
  
  if(biennial == 1) #biennial with higher peak in last year
  {  
    peak_1 <- peaks[(nrow(peaks)-1),1] #value of second-last peak
    peak_2 <-peaks[(nrow(peaks)-2),1] #value of third-last peak
  
    p_index <- ifelse(peak_2>peak_1, peaks[(nrow(peaks)-2),2], peaks[(nrow(peaks)-1),2])
    index_start <- p_index - (nYears - 1)*12 - (month - 1)
    index_end <- p_index + 12 - month
    return(list("start" = index_start, "end" = index_end))
  }
  
  #if biennial == 2 lower peak last
  peak_1 <- peaks[(nrow(peaks)-1),1] #value of second-last peak
  peak_2 <-peaks[(nrow(peaks)-2),1] #value of third-last peak
  
  p_index <- ifelse(peak_2>peak_1, peaks[(nrow(peaks)-1),2], peaks[(nrow(peaks)-2),2])
  index_start <- p_index - (nYears - 1)*12 - (month - 1)
  index_end <- p_index + 12 - month
  
  return(list("start" = index_start, "end" = index_end))
}
############################################################################
#'Find phase shift relating to a set "peak" month, default is 7 (July)
#'
#' @param counts table of modelled detected infections in aggregated age groups,
#'               when phi = 0
#' @param month month 1 - 12, Default is 7 (July)
#' @return numeric phase
#'
findPhaseShift <- function(counts, month = 7)
{
  dur <- nrow(counts)

  aggDetInc <- apply(counts, 1, sum)
  peaks <- pracma::findpeaks(aggDetInc)

  #Find last year of counts
  start_year <- (trunc(dur/12) - 1) *12 + 1
  end_year <- start_year + 11 
  int <- intersect(peaks[,2], start_year:end_year)
  cur_peak <- which((start_year:end_year) == int)

  if(cur_peak > month){
    phase <- (2 * pi / 12) * (cur_peak - month)
  }else{
    phase <- - (2 * pi / 12)*(month - cur_peak)
  }
   
  return(phase)
}
####################################################################
#' Set up likelihood function for risk model in a closure environment
#' for use with lazymcmc
#'
#' @param parTab parameter table setup for use in lazymcmc fitting functions
#' @param data observational incidence data
#' @param fixedPar fixed parameter values
#' @param PRIOR_FUNC optional prior function
#' @return function
create_likelihood_risk <- function(parTab, data, fixedPar, PRIOR_FUNC,...){

  par_names <- parTab$names

  likelihood_func  <- function(pars){

    names(pars) <- par_names

    times     <- seq(0, fixedPar$max_t, 0.25)
    old    <- ifelse(ncol(data) == 10, 1, 0)

    #Initialise y - state variables
    y0_seir <- matrix(0, fixedPar$nAges, 16)
    y0_seir[,1] <- (1 - fixedPar$alpha) * 0.99 * fixedPar$pop #S0 term
    y0_seir[,3] <- (1 - fixedPar$alpha) * 0.01 * fixedPar$pop #I0 term
    y0_seir[,9] <- fixedPar$alpha * 0.99 * fixedPar$pop #S0_bar preterm
    y0_seir[,11] <- fixedPar$alpha * 0.01 * fixedPar$pop #I0_bar preterm
    y0 <- cbind(y0_seir, matrix(0, fixedPar$nAges, 14))

    #set parameter values
    pars_ode <- list(
      b0 = pars["b0"],
      b1 = pars["b1"],
      phi = pars["phi"],
      delta = fixedPar$delta,
      gamma0 = fixedPar$gamma0,
      gamma1 = fixedPar$gamma1,
      nu = fixedPar$nu,
      omega_vect = fixedPar$omega_vect,
      omegaE = fixedPar$omegaE,
      alpha = fixedPar$alpha,
      A = pars["A"],
      B = pars["B"],
      C = pars["C"],
      D = pars["D"],
      E = pars["E"],
      age_in = fixedPar$age_in,
      age_out = fixedPar$age_out,
      age_months = fixedPar$age_months,
      sigma_vect = fixedPar$sigma_vect,
      sigmaE = fixedPar$sigmaE,
      mixing = fixedPar$mixing,
      nAges = fixedPar$nAges)

    deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm, pars_ode)

    modelOut <- aggregate_output(deSolveOut = deSolve_out,
                                 times = times,
                                 nAges = fixedPar$nAges,
                                 model = "risk",
                                 old = old)

    # "Match-up" model output with observational data
    mrange <- findIndexRange(counts = modelOut$combined, month = 7, obs = data, 
                             biennial = fixedPar$biennial)

    DetInc_1_3T = modelOut$term[mrange$start:mrange$end,"<3months"]
    DetInc_4_6T = modelOut$term[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12T = modelOut$term[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24T = modelOut$term[mrange$start:mrange$end,"12-<24months"]
    if(old) DetInc_25_60T = modelOut$term[mrange$start:mrange$end,"24-<60months"]

    DetInc_1_3PT = modelOut$pre[mrange$start:mrange$end,"<3months"]
    DetInc_4_6PT = modelOut$pre[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12PT = modelOut$pre[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24PT = modelOut$pre[mrange$start:mrange$end,"12-<24months"]
    if(old) DetInc_25_60PT = modelOut$pre[mrange$start:mrange$end,"24-<60months"]

    # calculate the log-likelihood
    loglikeli_term <- sum(data$`<3monthsT` * log(DetInc_1_3T) - DetInc_1_3T) +
      sum(data$`3-<6monthsT` * log(DetInc_4_6T) - DetInc_4_6T) +
      sum(data$`6-<12monthsT` * log(DetInc_7_12T) - DetInc_7_12T) +
      sum(data$`12-<24monthsT` * log(DetInc_13_24T) - DetInc_13_24T)
    if(old) loglikeli_term <- loglikeli_term + sum(data$`24-<60monthsT` * log(DetInc_25_60T) - DetInc_25_60T)

    loglikeli_preterm <- sum(data$`<3monthsPT` * log(DetInc_1_3PT) - DetInc_1_3PT) +
      sum(data$`3-<6monthsPT` * log(DetInc_4_6PT) - DetInc_4_6PT) +
      sum(data$`6-<12monthsPT` * log(DetInc_7_12PT) - DetInc_7_12PT) +
      sum(data$`12-<24monthsPT` * log(DetInc_13_24PT) - DetInc_13_24PT)
    if(old) loglikeli_preterm <- loglikeli_preterm + sum(data$`24-60monthsPT` * log(DetInc_25_60PT) - DetInc_25_60PT)

    loglikeli <- loglikeli_term + loglikeli_preterm

    if(!is.null(PRIOR_FUNC)) loglikeli <- loglikeli + PRIOR_FUNC(pars)

    loglikeli
  }
  return(likelihood_func)
}

####################################################################
#' Set up likelihood function for base model in a closure environment
#' for use with lazymcmc
#'
#' @param parTab parameter table setup for use in lazymcmc fitting functions
#' @param data observational incidence data
#' @param fixedPar fixed parameter values
#' @param PRIOR_FUNC optional prior function
#' @return function
create_likelihood_base <- function(parTab, data, fixedPar, PRIOR_FUNC,...){

  par_names <- parTab$names

  likelihood_func  <- function(pars){
  
    names(pars) <- par_names

    times     <- seq(0, fixedPar$max_t, 0.25)
    old <- 0
    if(ncol(data) == 5) old <- 1
    if(ncol(data) == 7) old <- 2

    #Initialise y - state variables
    y0_seir <- matrix(0, fixedPar$nAges, 8)
    y0_seir[,1] <- 0.99 * fixedPar$pop #S0 term
    y0_seir[,3] <- 0.01 * fixedPar$pop #I0 term
    y0 <- cbind(y0_seir, matrix(0, fixedPar$nAges, 6))

    #set parameter values
    pars_ode <- list(
      b0 = pars["b0"],
      b1 = pars["b1"],
      phi = pars["phi"],
      delta = fixedPar$delta,
      gamma0 = fixedPar$gamma0,
      gamma1 = fixedPar$gamma1,
      nu = fixedPar$nu,
      omega_vect = fixedPar$omega_vect,
      omegaE = fixedPar$omegaE,
      A = pars["A"],
      B = pars["B"],
      C = pars["C"],
      D = pars["D"],
      age_in = fixedPar$age_in,
      age_out = fixedPar$age_out,
      age_months = fixedPar$age_months,
      sigma_vect = fixedPar$sigma_vect,
      sigmaE = fixedPar$sigmaE,
      mixing = fixedPar$mixing,
      nAges = fixedPar$nAges)

    deSolve_out <- deSolve::ode(y0, times, deSolve_base, pars_ode)

    modelOut <- aggregate_output(deSolveOut = deSolve_out,
                                 times = times,
                                 nAges = fixedPar$nAges,
                                 model = "base",
                                 old = old)

    # "Match-up" model output with observational data
    mrange <- findIndexRange(counts = modelOut$combined, month = 7, obs = data, 
                             biennial = fixedPar$biennial)

    DetInc_1_3 = modelOut$combined[mrange$start:mrange$end,"<3months"]
    DetInc_4_6 = modelOut$combined[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12 = modelOut$combined[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24 = modelOut$combined[mrange$start:mrange$end,"12-<24months"]
    if(old == 1) DetInc_25_60 = modelOut$combined[mrange$start:mrange$end,"24-<60months"]
    if(old == 2){
      DetInc_25_36 <- modelOut$combined[mrange$start:mrange$end,"24-<36months"]
      DetInc_37_48 <- modelOut$combined[mrange$start:mrange$end,"36-<48months"]
      DetInc_49_60 <- modelOut$combined[mrange$start:mrange$end,"48-<60months"]
    }

    # calculate the log-likelihood
    loglikeli <- sum(data$`<3months` * log(DetInc_1_3) - DetInc_1_3) +
      sum(data$`3-<6months` * log(DetInc_4_6) - DetInc_4_6) +
      sum(data$`6-<12months` * log(DetInc_7_12) - DetInc_7_12) +
      sum(data$`12-<24months` * log(DetInc_13_24) - DetInc_13_24)
    if(old == 1) loglikeli <- loglikeli + sum(data$`24-<60months` * log(DetInc_25_60) - DetInc_25_60)
    if(old == 2) loglikeli <- loglikeli + sum(data$`24-<36months` * log(DetInc_25_36) - DetInc_25_36) +
                                          sum(data$`36-<48months` * log(DetInc_37_48) - DetInc_37_48) +
                                          sum(data$`48-<60months` * log(DetInc_49_60) - DetInc_49_60)

    if(!is.null(PRIOR_FUNC)) loglikeli <- loglikeli + PRIOR_FUNC(pars)

    loglikeli
  }
  return(likelihood_func)
}
####################################################################
#' Set up likelihood function for base model in a closure environment
#' for use with lazymcmc - this likelihood includes the assumption that most 
#' children by the age of 2 year old have been infected
#'
#' @param parTab parameter table setup for use in lazymcmc fitting functions
#' @param data observational incidence data
#' @param fixedPar fixed parameter values
#' @param PRIOR_FUNC optional prior function
#' @return function
create_likelihood_base_new <- function(parTab, data, fixedPar, PRIOR_FUNC,...){
  
  par_names <- parTab$names
  
  likelihood_func  <- function(pars){
    
    names(pars) <- par_names
    
    times     <- seq(0, fixedPar$max_t, 0.25)
    
    #Initialise y - state variables
    y0_seir <- matrix(0, fixedPar$nAges, 8)
    y0_seir[,1] <- 0.99 * fixedPar$pop #S0 term
    y0_seir[,3] <- 0.01 * fixedPar$pop #I0 term
    y0 <- cbind(y0_seir, matrix(0, fixedPar$nAges, 6))
    
    #set parameter values
    pars_ode <- list(
      b0 = pars["b0"],
      b1 = pars["b1"],
      phi = pars["phi"],
      delta = fixedPar$delta,
      gamma0 = fixedPar$gamma0,
      gamma1 = fixedPar$gamma1,
      nu = fixedPar$nu,
      omega_vect = fixedPar$omega_vect,
      omegaE = fixedPar$omegaE,
      A = pars["A"],
      B = pars["B"],
      C = pars["C"],
      D = pars["D"],
      age_in = fixedPar$age_in,
      age_out = fixedPar$age_out,
      age_months = fixedPar$age_months,
      sigma_vect = fixedPar$sigma_vect,
      sigmaE = fixedPar$sigmaE,
      mixing = fixedPar$mixing,
      nAges = fixedPar$nAges)
    
    deSolve_out <- deSolve::ode(y0, times, deSolve_base, pars_ode)
    
    modelOut <- aggregate_output(deSolveOut = deSolve_out,
                                 times = times,
                                 nAges = fixedPar$nAges,
                                 model = "base",
                                 old = 1)#1 = include older age group 24-<60 months
    
    # "Match-up" model output with observational data
    mrange <- findIndexRange(counts = modelOut$combined, month = 7, obs = data, 
                             biennial = fixedPar$biennial)
    
    DetInc_1_3 = modelOut$combined[mrange$start:mrange$end,"<3months"]
    DetInc_4_6 = modelOut$combined[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12 = modelOut$combined[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24 = modelOut$combined[mrange$start:mrange$end,"12-<24months"]
    DetInc_25_60 = modelOut$combined[mrange$start:mrange$end,"24-<60months"]
 
    #Find number of 2 year olds that have never been infected
    ageGroup_2yo <- which(fixedPar$age_months ==24)
    S0_2yo <- modelOut$deSolve[mrange$start:mrange$end,ageGroup_2yo,1]
    S0_5percent <- rep(0.05 * fixedPar$pop[ageGroup_2yo], 
                       length(mrange$start:mrange$end))
    
    # calculate the log-likelihood
    loglikeli <- sum(data$`<3months` * log(DetInc_1_3) - DetInc_1_3) +
      sum(data$`3-<6months` * log(DetInc_4_6) - DetInc_4_6) +
      sum(data$`6-<12months` * log(DetInc_7_12) - DetInc_7_12) +
      sum(data$`12-<24months` * log(DetInc_13_24) - DetInc_13_24)+ 
      sum(data$`24-<60months` * log(DetInc_25_60) - DetInc_25_60) +
      sum(S0_5percent * log(S0_2yo) - S0_2yo)
     
    if(!is.null(PRIOR_FUNC)) loglikeli <- loglikeli + PRIOR_FUNC(pars)
    
    loglikeli
  }
  return(likelihood_func)
}


####################################################################
#' Find average age of first infection/hospitalisation, second+
#' infection/hospitalisation for all, term and preterm
#'
#' @param modOut  model output returned from aggregate_output function
#' @param startT start t
#' @param endT end t
#' @param ageYears lower bound of age groups in years
#' @param risk if 1, then risk model, 0 for base
#' @return function
findAverageAge <- function(modOut, ageYears, startT, endT, risk = 1)
{
#NEED TO CHECK - ONLY CODED FOR RISK MODEL
  modOut_mat <- modOut$deSolve[startT:endT,,]

  if(risk)
  {  nStates <- 16
     modOut_names <-  c("dS0",
                     "dE0",
                     "dI0",
                     "dR0",
                     "dS1",
                     "dE1",
                     "dI1",
                     "dR1",
                     "dS0_bar",
                     "dE0_bar",
                     "dI0_bar",
                     "dR0_bar",
                     "dS1_bar",
                     "dE1_bar",
                     "dI1_bar",
                     "dR1_bar",
                     "dinc0", #term first infections
                     "ddetinc0",
                     "dinc1",
                     "ddetinc1",
                     "dincT",
                     "ddetincT",
                     "dinc0_bar", #preterm first infections
                     "ddetinc0_bar",
                     "dinc1_bar",
                     "ddetinc1_bar",
                     "dincPT",
                     "ddetincPT",
                     "dinc",
                     "ddetinc")
  }else{
    nStates <- 8
    modOut_names <-  c("dS0",
                     "dE0",
                     "dI0",
                     "dR0",
                     "dS1",
                     "dE1",
                     "dI1",
                     "dR1",
                     "dinc0",
                     "ddetinc0",
                     "dinc1",
                     "ddetinc1",
                     "dinc",
                     "ddetinc")
  }

  dimnames(modOut_mat)[[1]] <- 1:length(startT:endT) #t in months
  dimnames(modOut_mat)[[2]] <- age$age_years*12 #age in months
  dimnames(modOut_mat)[[3]] <- modOut_names

  midPointAge <- ageYears + (0.5*(1/12)) #in years

  #Only look at children under 5
  modOut_mat <- modOut_mat[, which(midPointAge < 5),]

  if(risk)
  {
    table_av <- matrix(NA, nrow = 6, ncol = 3)
    colnames(table_av) <- c("Term", "Preterm", "Comb")
    rownames(table_av) <- c("1st Inf", "1st Hosp",
                          "2+ Inf", "2+ Hosp",
                          "Inf", "Hosp")
    results_array <- rep(NA, 14)
    names(results_array) <- modOut_names[(nStates+1):(nStates+14)]

    for(i in 1:14)
    {
      freqAge <- t(apply(modOut_mat[,,modOut_names[nStates+i]], 1, function(x) x*midPointAge[1:60]))
      sumfreqAge <- apply(freqAge, 1, sum)
      numInfec <- apply(modOut_mat[,,modOut_names[nStates+i]], 1, sum)
      avMonth <- sumfreqAge/numInfec
      avMonth[which(is.na(avMonth))] <- 0
      results_array[i] <- sum(avMonth)/length(startT:endT)
    }

    table_av[, "Term"] <- results_array[1:6]
    table_av[, "Preterm"] <- results_array[7:12]
    table_av[c("Inf", "Hosp"), "Comb"] <- results_array[13:14]

    comb <- c("dinc0", "dinc0_bar",
            "ddetinc0", "ddetinc0_bar",
            "dinc1", "dinc1_bar",
            "ddetinc1", "ddetinc1_bar")

    results_comb <- rep(NA, 4)
    for(i in 1:4)
   {
    comb_mat<- modOut_mat[,,comb[i*2-1]] + modOut_mat[,,comb[i*2]]
    freqAge <- t(apply(comb_mat, 1, function(x) x*midPointAge[1:60]))
    sumfreqAge <- apply(freqAge, 1, sum)
    numInfec <- apply(comb_mat, 1, sum)
    avMonth <- sumfreqAge/numInfec
    avMonth[which(is.na(avMonth))] <- 0
    results_comb[i] <- sum(avMonth)/length(startT:endT)
  }

  table_av[1:4, "Comb"] <- results_comb

  return(table_av)
  }
}

####################################################################
#' Calculate log likelihood
#'
#' @param modelOut model output from aggregate_output function
#' @param data observational incidence data
#' @param risk 1 if risk model, 0 for base
#' @param biennial if 0 then observed series is annual, if 1 then observed series 
#'                 is biennial with higher peak in the last year, if 2 then
#'                 observed series is biennial with lower peak in the last year
#' @return numeric
calc_loglikelihood <- function(modelOut, data, risk = 1, biennial = 0)
{
  # "Match-up" model output with observational data
  mrange <- findIndexRange(counts = modelOut$combined, month = 7, obs = data,
                           biennial = biennial)
  
  if(risk)
  {
    DetInc_1_3T = modelOut$term[mrange$start:mrange$end,"<3months"]
    DetInc_4_6T = modelOut$term[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12T = modelOut$term[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24T = modelOut$term[mrange$start:mrange$end,"12-<24months"]
    DetInc_25_60T = modelOut$term[mrange$start:mrange$end,"24-<60months"]

    DetInc_1_3PT = modelOut$pre[mrange$start:mrange$end,"<3months"]
    DetInc_4_6PT = modelOut$pre[mrange$start:mrange$end,"3-<6months"]
    DetInc_7_12PT = modelOut$pre[mrange$start:mrange$end,"6-<12months"]
    DetInc_13_24PT = modelOut$pre[mrange$start:mrange$end,"12-<24months"]
    DetInc_25_60PT = modelOut$pre[mrange$start:mrange$end,"24-<60months"]
 
    data <- as.data.frame(data)
    # calculate the log-likelihood
    loglikeli_term <- sum(data$`<3monthsT` * log(DetInc_1_3T) - DetInc_1_3T) +
      sum(data$`3-<6monthsT` * log(DetInc_4_6T) - DetInc_4_6T) +
      sum(data$`6-<12monthsT` * log(DetInc_7_12T) - DetInc_7_12T) +
      sum(data$`12-<24monthsT` * log(DetInc_13_24T) - DetInc_13_24T) +
      sum(data$`24-<60monthsT` * log(DetInc_25_60T) - DetInc_25_60T)

    loglikeli_preterm <- sum(data$`<3monthsPT` * log(DetInc_1_3PT) - DetInc_1_3PT) +
      sum(data$`3-<6monthsPT` * log(DetInc_4_6PT) - DetInc_4_6PT) +
      sum(data$`6-<12monthsPT` * log(DetInc_7_12PT) - DetInc_7_12PT) +
      sum(data$`12-<24monthsPT` * log(DetInc_13_24PT) - DetInc_13_24PT) +
      sum(data$`24-<60monthsPT` * log(DetInc_25_60PT) - DetInc_25_60PT)
   
    loglikeli <- loglikeli_term + loglikeli_preterm

    return(loglikeli)
  }

  DetInc_1_3 = modelOut$combined[mrange$start:mrange$end,"<3months"]
  DetInc_4_6 = modelOut$combined[mrange$start:mrange$end,"3-<6months"]
  DetInc_7_12 = modelOut$combined[mrange$start:mrange$end,"6-<12months"]
  DetInc_13_24 = modelOut$combined[mrange$start:mrange$end,"12-<24months"]
  DetInc_25_60 = modelOut$combined[mrange$start:mrange$end,"24-<60months"]
 
  data <- as.data.frame(data)
  # calculate the log-likelihood
  loglikeli <- sum(data$`<3months` * log(DetInc_1_3) - DetInc_1_3) +
    sum(data$`3-<6months` * log(DetInc_4_6) - DetInc_4_6) +
    sum(data$`6-<12months` * log(DetInc_7_12) - DetInc_7_12) +
    sum(data$`12-<24months` * log(DetInc_13_24) - DetInc_13_24) + 
    sum(data$`24-<60months` * log(DetInc_25_60) - DetInc_25_60)

  return(loglikeli)
}

####################################################################
#' Find average age of first infection/hospitalisation, second+
#' infection/hospitalisation for all, term and preterm
#' NOTE: Only coded for risk model
#'
#' @param modOut  model output returned from aggregate_output function or
#'                observed hospitalisations in table form
#' @param startT start month, default NA for observed hospitalisations
#' @param endT end month, default NA for observed hospitalisations
#' @param ageYears lower bound of age groups in years
#' @param ageRange age range of each age group, default is 1 month
#' @param risk if 1, then risk model, 0 for base #not coded for 0 yet
#' @param obs if 1, then finding average age in observed hospitalisations
#'            (combined terms and preterms)
#' @return list
# findAverageAgePerMonth <- function(modOut, ageYears, ageRange = 1/12,
#                                    startT = NA, endT = NA, risk = 1, obs = 0)
# {
#   if(obs)
#   {
#     midPointAge <- ageYears + (0.5*ageRange) #in years
# 
#     freqAge <- t(apply(modOut, 1, function(x) x*midPointAge))
#     sumfreqAge <- apply(freqAge, 1, sum)
#     numInfec <- apply(modOut, 1, sum)
#     avMonth <- sumfreqAge/numInfec
# 
#     #find average age over peak months June - August and over year
#     totalMonths <- nrow(modOut)
#     avPeak <- rep(NA, totalMonths)
#     avYear <- rep(NA, totalMonths)
#     for(i in 1:(totalMonths/12))
#     {
#      peakMonths <- ((i* 12 - 12) + 6) : ((i* 12 - 12) + 8)
#      monthsYear <- (i* 12 - 11) : (i* 12)
# 
#      numInfecPeak <- sum(numInfec[peakMonths])
#      avPeak[monthsYear] <- sum(sumfreqAge[peakMonths])/numInfecPeak
# 
#      numInfecYear <- sum(numInfec[monthsYear])
#      avYear[monthsYear] <- sum(sumfreqAge[monthsYear])/numInfecYear
# 
#     }
# 
#     avMonthNoNA <- avMonth
#     avMonthNoNA[which(is.na(avMonthNoNA))] <- 0
# 
#     return(list("AvMonth" = avMonth, "avMonthNoNA" = avMonthNoNA,
#                 "AvPeak" = avPeak, "AvYear" = avYear))
#   }
# 
#   if(risk)nStates <- 16
# 
#   modOut_mat <- modOut$deSolve[startT:endT,,]
# 
#   modOut_names <-  c("dS0",
#                      "dE0",
#                      "dI0",
#                      "dR0",
#                      "dS1",
#                      "dE1",
#                      "dI1",
#                      "dR1",
#                      "dS0_bar",
#                      "dE0_bar",
#                      "dI0_bar",
#                      "dR0_bar",
#                      "dS1_bar",
#                      "dE1_bar",
#                      "dI1_bar",
#                      "dR1_bar",
#                      "dinc0", #term first infections
#                      "ddetinc0",
#                      "dinc1",
#                      "ddetinc1",
#                      "dincT",
#                      "ddetincT",
#                      "dinc0_bar", #preterm first infections
#                      "ddetinc0_bar",
#                      "dinc1_bar",
#                      "ddetinc1_bar",
#                      "dincPT",
#                      "ddetincPT",
#                      "dinc",
#                      "ddetinc")
# 
#   dimnames(modOut_mat)[[1]] <- 1:length(startT:endT) #t in months
#   dimnames(modOut_mat)[[2]] <- ageYears*12 #age in months
#   dimnames(modOut_mat)[[3]] <- modOut_names
# 
#   midPointAge <- ageYears + (0.5*ageRange) #in years
# 
#   #Only look at children under 5
#   modOut_mat <- modOut_mat[, which(midPointAge < 5),]
# 
#   table_av <- matrix(NA, nrow = 6, ncol = length(startT:endT))
#   colnames(table_av) <- 1:length(startT:endT)
#   rownames(table_av) <- c("1st Inf",
#                           "1st Hosp",
#                           "2+ Inf",
#                           "2+ Hosp",
#                           "Inf",
#                           "Hosp")
# 
#   comb <- c("dinc0", "dinc0_bar",
#             "ddetinc0", "ddetinc0_bar",
#             "dinc1", "dinc1_bar",
#             "ddetinc1", "ddetinc1_bar")
# 
#   for(i in 1:4)
#   {
#     comb_mat<- modOut_mat[,,comb[i*2-1]] + modOut_mat[,,comb[i*2]]
#     freqAge <- t(apply(comb_mat, 1, function(x) x*midPointAge[1:60]))
#     sumfreqAge <- apply(freqAge, 1, sum)
#     numInfec <- apply(comb_mat, 1, sum)
#     avMonth <- sumfreqAge/numInfec
#     avMonth[which(is.na(avMonth))] <- 0
#     table_av[i,] <- avMonth
#   }
# 
# 
#   freqAge <- t(apply(modOut_mat[,,"dinc"], 1, function(x) x*midPointAge[1:60]))
#   sumfreqAge <- apply(freqAge, 1, sum)
#   numInfec <- apply(modOut_mat[,,"dinc"], 1, sum)
#   avMonth <- sumfreqAge/numInfec
#   avMonth[which(is.na(avMonth))] <- 0
#   table_av[5,] <- avMonth
# 
#   freqAge <- t(apply(modOut_mat[,,"ddetinc"], 1, function(x) x*midPointAge[1:60]))
#   sumfreqAge <- apply(freqAge, 1, sum)
#   numInfec <- apply(modOut_mat[,,"ddetinc"], 1, sum)
#   avMonth <- sumfreqAge/numInfec
#   avMonth[which(is.na(avMonth))] <- 0
#   table_av[6,] <- avMonth
# 
#   hosp_Under5 <- apply(modOut_mat[,,"ddetinc"], 1, sum)
# 
#   return(list("table_av" = table_av, "hosp_Under5" = hosp_Under5))
# }