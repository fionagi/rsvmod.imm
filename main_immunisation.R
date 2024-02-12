rm(list = ls())

library(deSolve)
library(conmat)
library(dplyr)
library(pracma)
library(ggplot2)

setwd("C:/Users/FGiannini/OneDrive/TKI 2021/RSV/R code/rsvmod - exploratory/")
load("data/data.greaterPerth.2011.rda")
load("data/data.greaterPerth.2021.rda")
source("R/functions.r")
source("R/parameters.r")
source("R/models.r")
source("R/immunisation.r")
source("R/diagnostic.r")

age <- age_structure(year = "2011")
mixing <- as.matrix(read.csv("data/conmat_2011UKSym_daily.csv", header = TRUE))*365/12

#create vectors
#reduced infectiousness - currently assuming no reduction
omega_vect <- as.vector(rep(1, age$nAges))
omega_vect[!(age$age_years < 10)] <- omega
#reduced susceptibility due to age - natural maternal immunity
sigma_vect <- as.vector(c(sigma, rep(1, (age$nAges-6))))

# output time steps
run_in <- 12*12 #run-in the model for 12 years
model_time <- 13*12 #then run model for 13 years
n_times <- run_in + model_time
times     <- seq(0, n_times-1, 0.25)

#100% mAb coverage
kappaM <- 1

#Year-round mAb at birth (after model run-in)
kappaM_yearRound <- c(rep(0, run_in), rep(kappaM, model_time))

#Seasonal mAb at birth (after model run-in)
seasonal <- rep(0, 12) # Jan - Dec
#seasonal[4:9] <- 1 #mAb at birth for those born in April - Sept.
seasonal[3:8] <- 1 #mAb at birth for those born in April - Sept.
                   #NEED TO CHANGE EQUATIONS TIME STEP BEFORE
kappaM_seasonal <-  c(rep(0, run_in), kappaM * rep(seasonal, model_time/12))

#Turn off coverage as model check
kappaM_zero <- rep(0, n_times)

#Seasonal with catch-up
kappaM_seasWithCatch <- c(rep(0, run_in), kappaM * rep(seasonal, model_time/12))
catchUp <- rep(0, 12) # Jan - Dec
catchUp[3] <- 1 #mAb catchup month April
kappaM_catchUpApril <- c(rep(0, run_in), kappaM * rep(catchUp, model_time/12))
age_catchUp <- 5 #Up to and including 5 month olds in catchup month [0,6)


kappaM_noCatchUp <- rep(0, n_times)

#Second RSV season dose
kappaM_dose2 <- rep(0, n_times)
dose2Age <- 8

numImm <- 20 #number of compartments per age group for base immunisation model
             # SEIR x 2 exposures x 2 states (protected or not) + (k - 1) x 2 exposures
             # where k is number of protected Susceptible compartments per exposure
numOutB <- 18

######################################################################################
# #Write function to fit phase
# y0_seir <- matrix(0, age$nAges, 8)
# y0_seir[,1] <- 0.99 * age$pop_groups #S0 term
# y0_seir[,3] <- 0.01 * age$pop_groups #I0 term
# y0 <- cbind(y0_seir, matrix(0, age$nAges, 6))
#
# #set parameter values
# pars_base <- list(
#   b0 = b0_fit,
#   b1 = b1_fit,
#   phi = 0,
#   delta = delta,
#   gamma0 = gamma0,
#   gamma1 = gamma1,
#   nu = nu,
#   omega_vect = omega_vect,
#   omegaE = omegaE,
#   A = A_fit,
#   B = B_fit,
#   C = C_fit,
#   D = D_fit,
#   age_in = age$age_in,
#   age_out = age$age_out,
#   age_months = age$age_years*12,
#   sigma_vect = sigma_vect,
#   sigmaE = sigmaE,
#   mixing = mixing,
#   nAges = age$nAges)
#
# deSolve_out <- deSolve::ode(y0, times, deSolve_base, pars_base)
#
# modelOut_base <- aggregate_output(deSolveOut = deSolve_out,
#                                   times = times,
#                                   nAges = age$nAges,
#                                   model = "base",
#                                   old = 1)
# phase <- findPhaseShift(counts = modelOut_base$combined, month = 7)
#
######################################################################################
#Run ORIGINAL BASE MODEL (no immunisation)
#Initialise y - state variables - pre-terms
y0_seir <- matrix(0, age$nAges, 8)
y0_seir[,1] <- 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- 0.01 * age$pop_groups #I0 term
y0 <- cbind(y0_seir, matrix(0, age$nAges, 6))

#set parameter values
pars_base <- list(
  b0 = b0_fit,
  b1 = b1_fit,
  phi = phi_fit,
  delta = delta,
  gamma0 = gamma0,
  gamma1 = gamma1,
  nu = nu,
  omega_vect = omega_vect,
  omegaE = omegaE,
  A = A_fit,
  B = B_fit,
  C = C_fit,
  D = D_fit,
  age_in = age$age_in,
  age_out = age$age_out,
  age_months = age$age_years*12,
  sigma_vect = sigma_vect,
  sigmaE = sigmaE,
  mixing = mixing,
  nAges = age$nAges)

deSolve_out <- deSolve::ode(y0, times, deSolve_base, pars_base)

modelOut_base <- aggregate_output(deSolveOut = deSolve_out,
                                  times = times,
                                  nAges = age$nAges,
                                  model = "base",
                                  old = 1)

comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_base$combined[(run_in+1):n_times,])

#####################################################################################
#Run BASE IMMUNISATION MODEL with coverage turned off
#Initialise y - state variables - pre-terms
y0_seir <- matrix(0, age$nAges, numImm)
y0_seir[,1] <- 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- 0.01 * age$pop_groups #I0 term
y0 <- cbind(y0_seir, matrix(0, age$nAges, numOutB))

#set parameter values
pars_base_imm <- list(
  b0 = b0_fit,
  b1 = b1_fit,
  phi = phi_fit,
  delta = delta,
  gamma0 = gamma0,
  gamma1 = gamma1,
  nu = nu,
  nuM = nuM,
  omega_vect = omega_vect,
  omegaE = omegaE,
  kappaM_Birth = kappaM_zero,
  kappaM_catchUp = kappaM_noCatchUp,
  catchAge = age_catchUp,
  kappaM_dose2 = kappaM_dose2,
  dose2Age = dose2Age,
  A = A_fit,
  B = B_fit,
  C = C_fit,
  D = D_fit,
  rho = 0, #mAb efficacy
  age_in = age$age_in,
  age_out = age$age_out,
  age_months = age$age_years*12,
  sigma_vect = sigma_vect,
  sigmaE = sigmaE,
  mixing = mixing,
  nAges = age$nAges
)

deSolve_out <- deSolve::ode(y0, times, deSolve_base_imm, pars_base_imm)

modelOut_ImmNoCover <- aggregate_output_immi (deSolveOut = deSolve_out,
                                              times = times,
                                              nAges = age$nAges,
                                              nState = numImm,
                                              nRet = numImm + numOutB,
                                              old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_ImmNoCover$all_combined[(run_in+1):n_times,])

#####################################################################################
#Run BASE IMMUNISATION MODEL year-round mAb
#3 mAb compartments + average durability of 150 days gives efficacy of ~0.78

pars_base_imm$kappaM_Birth <- kappaM_yearRound
pars_base_imm$rho <- 1 #0 hospitalisation while in mAb compartment

deSolve_out <- deSolve::ode(y0, times, deSolve_base_imm, pars_base_imm)

modelOut_yearRound <- aggregate_output_immi (deSolveOut = deSolve_out,
                                              times = times,
                                              nAges = age$nAges,
                                              nState = numImm,
                                              nRet = numImm + numOutB,
                                              old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_yearRound$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb (colour .) and year round mAb at birth (black --)")

####################################################################################
#Run BASE IMMUNISATION MODEL seasonal mAb

pars_base_imm$kappaM_Birth <- kappaM_seasonal

deSolve_out <- deSolve::ode(y0, times, deSolve_base_imm, pars_base_imm)

modelOut_seasonal <- aggregate_output_immi (deSolveOut = deSolve_out,
                                             times = times,
                                             nAges = age$nAges,
                                             nState = numImm,
                                             nRet = numImm + numOutB,
                                             old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonal$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb (colour .) and seasonal at birth mAb (black --)")

#Compare year-round with seasonal
comp_data(obs1 = modelOut_yearRound$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonal$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: year-round mAb (colour .) and seasonal at birth mAb (black --)")

comp_monthly(obs1 = modelOut_base$combined[277:300, 1],
             obs2 = modelOut_seasonal$all_combined[277:300, 1],
             ylabel = "Hospitalisations children < 3months",
             text = "No immunisation red, Seasonal blue")

comp_monthly(obs1 = modelOut_yearRound$all_combined[277:300, 1],
             obs2 = modelOut_seasonalCUApril$all_combined[277:300, 1],
             obs3 = modelOut_seasonal$all_combined[277:300, 1],
             ylabel = "Hospitalisations children < 3months",
             text = "Year round black, Seasonal CU blue, Seasonal red")

####################################################################################
#Run BASE IMMUNISATION MODEL seasonal with catch-up mAb in
#April, all infants under 6 months
pars_base_imm$kappaM_Birth <- kappaM_seasonal #April - September at birth

pars_base_imm$kappaM_catchUp <- kappaM_catchUpApril
pars_base_imm$catchAge <- age_catchUp

deSolve_out <- deSolve::ode(y0, times, deSolve_base_imm, pars_base_imm)

modelOut_seasonalCUApril <- aggregate_output_immi (deSolveOut = deSolve_out,
                                              times = times,
                                              nAges = age$nAges,
                                              nState = numImm,
                                              nRet = numImm + numOutB,
                                              old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb (colour .) and seasonal with catchup April mAb (black --)")

comp_monthly(obs1 = modelOut_base$combined[277:300, 1],
             obs2 = modelOut_seasonalCUApril$all_combined[277:300, 1],
             ylabel = "Hospitalisations children < 3months",
             text = "No immunisation red, Seasonal with CU blue")

#Compare year-round with seasonal w catchup
comp_data(obs1 = modelOut_yearRound$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: year-round mAb (colour .) and seasonal with catchup April mAb (black --)")

comp_monthly(obs1 = modelOut_yearRound$all_combined[277:300, 1],
             obs2 = modelOut_seasonalCUApril$all_combined[277:300, 1],
             ylabel = "Hospitalisations children < 3months",
             text = "Year-round red, Seasonal with CU blue")

comp_monthly(obs1 = modelOut_base$combined[277:300, 1],
             obs2 = modelOut_yearRound$all_combined[277:300, 1],
             obs3 = modelOut_seasonalCUApril$all_combined[277:300, 1],
             ylabel = "Hospitalisations children < 3months",
             text = "No immunisation black, Year-round blue, Seasonal with CU red")

comp_monthly(obs1 = modelOut_base$combined[277:300, 4],
             obs2 = modelOut_yearRound$all_combined[277:300, 4],
             obs3 = modelOut_seasonalCUApril$all_combined[277:300, 4],
             ylabel = "Hospitalisations children > 3months, < 6months",
             text = "No immunisation black, Year-round blue, Seasonal with CU red")


#Compare seasonal with seasonal w catchup
comp_data(obs1 = modelOut_seasonal$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: seasonal mAb (colour .) and seasonal with catchup April mAb (black --)")

#Compare seasonal w catchup April with seasonal w catchup
comp_data(obs1 = modelOut_seasonalCU$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: seasonal with catchup May mAb (colour .) and seasonal with catchup April mAb (black --)")

#################################################################################
#Run BASE IMMUNISATION MODEL with second dose - US equivalent scenario
#Seasonal w catchup in April, catchup age <7 months for first season,
# <19 months for second season
pars_base_imm$catchAge <- 6
pars_base_imm$kappaM_dose2 <- kappaM_catchUpApril
pars_base_imm$dose2Age <- 18

deSolve_out <- deSolve::ode(y0, times, deSolve_base_imm, pars_base_imm)

modelOut_USequ <- aggregate_output_immi (deSolveOut = deSolve_out,
                                                   times = times,
                                                   nAges = age$nAges,
                                                   nState = numImm,
                                                   nRet = numImm + numOutB,
                                                   old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_USequ$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb (colour .) and US equiv mAb (black --)")


#Compare seasonal w catchup April with US equiv
comp_data(obs1 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_USequ$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: seasonal with catchup mAb (colour .) and US equiv mAb (black --)")

#Compare year-round with US equiv
comp_data(obs1 = modelOut_yearRound$all_combined[(run_in+1):n_times,],
          obs2 = modelOut_USequ$all_combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: year-round mAb (colour .) and US equiv mAb (black --)")

###############################################################################
#Run original RISK model (no immunisation)
#Initialise y - state variables - pre-terms
y0_seir <- matrix(0, age$nAges, 16)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0 <- cbind(y0_seir, matrix(0, age$nAges, 14))


#set parameter values
pars_pT <- list(b0 = b0_fit,
                b1 = b1_fit,
                phi = phi_fit,
                delta = delta,
                gamma0 = gamma0,
                gamma1 = gamma1,
                nu = nu,
                omega_vect = omega_vect,
                omegaE = omegaE,
                alpha = alpha,
                A = A_fit,
                B = B_fit,
                C = C_fit,
                D = D_fit,
                E = E_fit,
                age_in = age$age_in,
                age_out = age$age_out,
                age_months = age$age_years*12,
                sigma_vect = sigma_vect,
                sigmaE = sigmaE,
                mixing = mixing,
                nAges = age$nAges)

deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm, pars_pT)

modelOut_risk <- aggregate_output(deSolveOut = deSolve_out,
                                times = times,
                                nAges = age$nAges,
                                model = "risk",
                                old = 1)

#Compare output with original base model
comp_data(obs1 = modelOut_base$combined[(run_in+1):n_times,],
          obs2 = modelOut_risk$combined[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb base (colour .) and risk (black --)")

#################################################################################
#Run risk (preterm) immunisation model with coverage turned off to check
y0_seir <- matrix(0, age$nAges, numImm * 2)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0 <- cbind(y0_seir, matrix(0, age$nAges, 6))


pars_pT_imm <- list(b0 = b0_fit,
                b1 = b1_fit,
               phi = phi_fit,
             delta = delta,
            gamma0 = gamma0,
            gamma1 = gamma1,
                nu = nu,
               nuM = nuM, #durability of mAb
        omega_vect = omega_vect,
            omegaE = omegaE,
      kappaM_Birth = kappaM_zero, #coverage for terms
    kappaM_BirthPT = kappaM_zero, #coverage for preterms
    kappaM_catchUp = kappaM_zero, #coverage + month catch-up for terms
  kappaM_catchUpPT = kappaM_zero, #coverage + month catch-up preterms
          catchAge = 0, #catch-up covers birth to this age group for all
      kappaM_dose2 = kappaM_zero, #coverage + month for 2nd dose terms
    kappaM_dose2PT = kappaM_zero, #coverage + month for 2nd dose preterms
          dose2Age = 0, #2nd dose covers catchAge +1 to this age group for all
             alpha = alpha,
                 A = A_fit,
                 B = B_fit,
                 C = C_fit,
                 D = D_fit,
                 E = E_fit,
               rho = 0, #efficacy for terms
            rho_PT = 0, #efficacy for preterms
            age_in = age$age_in,
           age_out = age$age_out,
        age_months = age$age_years*12,
        sigma_vect = sigma_vect,
            sigmaE = sigmaE,
            mixing = mixing,
             nAges = age$nAges)

deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm_imm, pars_pT_immi)

modelOut_preTerm <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                    times = times,
                                                    nAges = age$nAges,
                                                    numImm = numImm,
                                                    old = 1)

#Compare output with original risk model - combined
comp_data(obs1 = modelOut_risk$combined[(run_in+1):n_times,],
          obs2 = modelOut_preTerm$hosp_comb[(run_in+1):n_times,],
          text = "RSV hospitalisations: original risk(colour .) and no mAb imm model (black --)")

#Terms only
comp_data(obs1 = modelOut_risk$term[(run_in+1):n_times,],
          obs2 = modelOut_preTerm$hos_term[(run_in+1):n_times,],
          text = "RSV hospitalisations: original term only (colour .) and no mAb imm model term only (black --)")

#Preterms only
comp_data(obs1 = modelOut_risk$pre[(run_in+1):n_times,],
          obs2 = modelOut_preTerm$hos_preterm[(run_in+1):n_times,],
          text = "RSV hospitalisations: original preterm only (colour .) and no mAb imm model preterm only (black --)")

##################################################################################
#Run RISK IMMUNISATION MODEL
#US equiv
#Seasonal w catchup in April, catchup age <7 months for first season ALL INFANTS,
# <19 months for second season ONLY PRETERMS
pars_pT_immi$kappaM_Birth <- kappaM_seasonal
pars_pT_immi$kappaM_BirthPT <- kappaM_seasonal
pars_pT_immi$kappaM_catchUp <- kappaM_catchUpApril
pars_pT_immi$kappaM_catchUpPT <- kappaM_catchUpApril
pars_pT_immi$catchAge <- 6
pars_pT_immi$kappaM_dose2 <- kappaM_zero
pars_pT_immi$kappaM_dose2PT <- kappaM_catchUpApril
pars_pT_immi$dose2Age <- 18
pars_pT_immi$rho <- 1
pars_pT_immi$rho_PT <- 1


deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm_imm1, pars_pT_immi)

modelOut_USequPT <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                             times = times,
                                             nAges = age$nAges,
                                             numImm = numImm,
                                             old = 1)

#Compare output with original risk model - combined
comp_data(obs1 = modelOut_risk$combined[(run_in+1):n_times,],
          obs2 = modelOut_USequPT$hosp_comb[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb (colour .) and preterm US equiv imm model (black --)")

#Terms only
comp_data(obs1 = modelOut_risk$term[(run_in+1):n_times,],
          obs2 = modelOut_USequPT$hos_term[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb term only (colour .) and no mAb imm model term only (black --)")

#Preterms only
comp_data(obs1 = modelOut_risk$pre[(run_in+1):n_times,],
          obs2 = modelOut_USequPT$hos_preterm[(run_in+1):n_times,],
          text = "RSV hospitalisations: no mAb preterm only (colour .) and no mAb imm model preterm only (black --)")




# #Compare year-round with seasonal w catchup
# comp_data(obs1 = modelOut_yearRound$all_combined[(run_in+1):n_times,],
#           obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
#           text = "RSV hospitalisations: year-round mAb (colour .) and seasonal with catchup April mAb (black --)")
#
# #Compare seasonal with seasonal w catchup
# comp_data(obs1 = modelOut_seasonal$all_combined[(run_in+1):n_times,],
#           obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
#           text = "RSV hospitalisations: seasonal mAb (colour .) and seasonal with catchup April mAb (black --)")
#
# #Compare seasonal w catchup April with seasonal w catchup
# comp_data(obs1 = modelOut_seasonalCU$all_combined[(run_in+1):n_times,],
#           obs2 = modelOut_seasonalCUApril$all_combined[(run_in+1):n_times,],
#           text = "RSV hospitalisations: seasonal with catchup May mAb (colour .) and seasonal with catchup April mAb (black --)")
#
#
