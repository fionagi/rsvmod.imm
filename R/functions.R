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

