rm(list = ls())

library(deSolve)
library(conmat)
library(dplyr)
library(pracma)
library(ggplot2)

setwd("C:/Users/FGiannini/OneDrive/TKI 2021/RSV/R code/rsvmod_imm/")
load("data/data.southernWA.2016.rda")
load("data/data.contactDaily.ABSMicro.rda")

source("R/data.r")
source("R/functions.r")
source("R/parameters.r")
source("R/models.r")
source("R/immunisation.r")
source("R/diagnostic.r")

##############################################################################
#Set up model
age <- age_structure(year = "2016", area = 1)
mixing <- data.contactDaily.ABSMicro * 365/12

#create vectors
#reduced infectiousness - currently assuming no reduction
omega_vect <- as.vector(rep(1, age$nAges))
omega_vect[!(age$age_years < 10)] <- omega
#reduced susceptibility due to natural maternal immunity
sigma_vect <- as.vector(c(sigma, rep(1, (age$nAges-12))))

# output time steps
run_in <- 12*12 #run-in the model for 12 years
model_time <- 5*12 #then run model for 5 years
n_times <- run_in + model_time
times     <- seq(0, n_times-1, 0.25)

##############################################################################
#MODEL FIT

#Plot age-to-risk function using fixed and fitted values
plot_ageToRisk(A = A_fit, B = B_fit, C = C_fit, D = D_fit,
               E = E_fit, maxAge = 24)


#Run original RISK MODEL (no immunisation, but stratification from preterms)
#Run risk model
#Initialise y - state variables - pre-terms
y0_seir <- matrix(0, age$nAges, 16)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0 <- cbind(y0_seir, matrix(0, age$nAges, 14))

#set parameter values
pars_pT <- list(
  b0 = b0_fit,
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

#Solve ODEs
deSolve_out <- deSolve::ode(y0, times, deSolve_preTerm, pars_pT)

#Aggregate output into age groups for plots
modelOut_pT <- aggregate_output(deSolveOut = deSolve_out,
                                times = times,
                                nAges = age$nAges,
                                model = "risk",
                                old = 1)

#Plot to compare model output to observed. Observed data can not be provided
#due to confidentiality. Here we have put in model output twice
dummyObs <- cbind(modelOut_pT$term[(run_in+1):n_times,1],
                  modelOut_pT$pre[(run_in+1):n_times,1] )
for(i in 2:5) dummyObs <- cbind(dummyObs,
                                modelOut_pT$term[(run_in+1):n_times,i],
                                modelOut_pT$pre[(run_in+1):n_times,i])

all_plots <- plot_riskFit(term = modelOut_pT$term[(run_in+1):n_times,],
                          pTerm = modelOut_pT$pre[(run_in+1):n_times,], 
                          obs = dummyObs,
                          years = 2015:2020)

#First plot gives preterm and term predictions on same plot
all_plots[[1]]
#Second, just term
all_plots[[2]]
#Third, just preterm
all_plots[[3]]

##############################################################################
#Set immunisation parameters

#Seasonal mAb at birth (after model run-in)
#WA program: At birth for all babies born between 1 May and 30 September 2024
seasonal <- rep(0, 12) # Jan - Dec
seasonal[4:8] <- 1 #mAb at birth for those born in May - Sept, adjusting for t = 0
atBirth_April <- seasonal
atBirth_April[3] <- 1 #mAb at birth for those born in April - Sept
kappaM_seasonal <-  c(rep(0, run_in), kappaM * rep(seasonal, model_time/12))
kappaM_seasonalApril <- c(rep(0, run_in), kappaM * rep(atBirth_April, model_time/12))

#Turn off coverage as model check
kappaM_zero <- rep(0, n_times)

#Catch-up mAb
#WA program: Babies born from 1 October 2023 to 30 April 2024 are able 
#to receive nirsevimab from 1 April 2024 to 30 September 2024
catchUp <- rep(0, 12) # Jan - Dec
catchUp[3] <- 1 #assume optimal catchup occurring in April
kappaM_catchUp <- c(rep(0, run_in), kappaM * rep(catchUp, model_time/12))

age_catchUp <- 6 #Up to 7 month olds in first catchup month April [6,7)

#Second RSV season dose
age_catchUp_2seasons <- 18 #Up to and including 18 month olds [0,19) 

numImm_Risk <- 40 #number of compartments per age group for immunisation model
# 4 SEIR x 2 exposures x 2 states (protected or not) x 2 stratifications (preterm/term) 
#+ (k - 1) x 2 stratifications x 2 exposures, where k is 
#number of protected Susceptible compartments per exposure
numOutB_Risk <- 6 #model output variables

#################################################################################
#Run immunisation model with no immunisation
y0_seir <- matrix(0, age$nAges, numImm_Risk)
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
                    rho = rho, #efficacy for terms
                    rho_PT = rho_PT, #efficacy for preterms
                    age_in = age$age_in,
                    age_out = age$age_out,
                    age_months = age$age_years*12,
                    sigma_vect = sigma_vect,
                    sigmaE = sigmaE,
                    mixing = mixing,
                    nAges = age$nAges)

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_NoVaxx <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                 times = times,
                                                 nAges = age$nAges,
                                                 nState = numImm_Risk,
                                                 nRet = numImm_Risk + numOutB_Risk,
                                                 old = 1)

##################################################################################
#Run RISK IMMUNISATION MODEL
#Seasonal-at-birth component of WA program May - September
#Run with 100% coverage and 85% coverage
#3 mAb compartments + average durability of 150 days gives efficacy of ~0.78

#100%
pars_pT_imm$kappaM_Birth <- kappaM_seasonal
pars_pT_imm$kappaM_BirthPT <- kappaM_seasonal

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirth100 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                             times = times,
                                                             nAges = age$nAges,
                                                             nState = numImm_Risk,
                                                             nRet = numImm_Risk + numOutB_Risk,
                                                             old = 1)

#WA coverage
pars_pT_imm$kappaM_Birth <- 0.85 * kappaM_seasonal
pars_pT_imm$kappaM_BirthPT <- 0.85 * kappaM_seasonal

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirth85 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                            times = times,
                                                            nAges = age$nAges,
                                                            nState = numImm_Risk,
                                                            nRet = numImm_Risk + numOutB_Risk,
                                                            old = 1)

##################################################################################
#Run RISK IMMUNISATION MODEL
#Seasonal-at-birth component of WA program May - September
#Include April at birth CU component
#100% coverage and 85% coverage
#3 mAb compartments + average durability of 150 days gives efficacy of ~0.78

#100%
pars_pT_imm$kappaM_Birth <- kappaM_seasonalApril
pars_pT_imm$kappaM_BirthPT <- kappaM_seasonalApril

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirthApril100 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                  times = times,
                                                                  nAges = age$nAges,
                                                                  nState = numImm_Risk,
                                                                  nRet = numImm_Risk + numOutB_Risk,
                                                                  old = 1)

#WA coverage
pars_pT_imm$kappaM_Birth <- 0.85 * kappaM_seasonalApril
pars_pT_imm$kappaM_BirthPT <- 0.85 * kappaM_seasonalApril

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirthApril85 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                 times = times,
                                                                 nAges = age$nAges,
                                                                 nState = numImm_Risk,
                                                                 nRet = numImm_Risk + numOutB_Risk,
                                                                 old = 1)

################################################################################
#Run RISK IMMUNISATION MODEL seasonal with catch-up mAb WA program
#100% coverage and
#At birth if born between May - September, catch-up 85% coverage
#catch-up in April for infants under 7 months in April 65% coverage (April at
#birth is 85% coverage)

#100%
pars_pT_imm$kappaM_Birth <- kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirthCU <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                            times = times,
                                                            nAges = age$nAges,
                                                            nState = numImm_Risk,
                                                            nRet = numImm_Risk + numOutB_Risk,
                                                            old = 1)

#WA coverage
pars_pT_imm$kappaM_Birth <- 0.85 * kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- 0.85 * kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirth85CU65 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                times = times,
                                                                nAges = age$nAges,
                                                                nState = numImm_Risk,
                                                                nRet = numImm_Risk + numOutB_Risk,
                                                                old = 1)

################################################################################
#Run RISK IMMUNISATION MODEL seasonal with catch-up mAb WA program, 2nd dose for PT
#At birth if born between May - September, catch-up in AT 85% coverage 
#catch-up in April for infants under 7 months in April AT 65% coverage (April at
#birth is 80% coverage)
#PT 2nd dose 25%
#PLUS the above with 100% coverage

#100%
pars_pT_imm$kappaM_Birth <- kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp

pars_pT_imm$kappaM_dose2PT <- kappaM_catchUp
pars_pT_imm$dose2Age <- age_catchUp_2seasons

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirthCU_2nd <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                times = times,
                                                                nAges = age$nAges,
                                                                nState = numImm_Risk,
                                                                nRet = numImm_Risk + numOutB_Risk,
                                                                old = 1)

#WA coverage
pars_pT_imm$kappaM_Birth <- 0.85 * kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- 0.85 * kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp

pars_pT_imm$kappaM_dose2PT <- 0.25 * kappaM_catchUp
pars_pT_imm$dose2Age <- age_catchUp_2seasons

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirth85CU65_2nd25 <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                      times = times,
                                                                      nAges = age$nAges,
                                                                      nState = numImm_Risk,
                                                                      nRet = numImm_Risk + numOutB_Risk,
                                                                      old = 1)

################################################################################
#Run RISK IMMUNISATION MODEL seasonal with catch-up mAb WA program, 2nd dose for ALL
#At birth if born between May - September, catch-up in AT 85% coverage 
#catch-up in April for infants under 7 months in April AT 65% coverage (April at
#birth is 80% coverage)
#2nd dose 25%
#PLUS the above with 100% coverage
pars_pT_imm$kappaM_Birth <- 0.85 * kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- 0.85 * kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- 0.65 * kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp

pars_pT_imm$kappaM_dose2 <-  0.25 * kappaM_catchUp
pars_pT_imm$kappaM_dose2PT <- 0.25 * kappaM_catchUp
pars_pT_imm$dose2Age <- age_catchUp_2seasons

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirth85CU65_2nd25ALL <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                         times = times,
                                                                         nAges = age$nAges,
                                                                         nState = numImm_Risk,
                                                                         nRet = numImm_Risk + numOutB_Risk,
                                                                         old = 1)

pars_pT_imm$kappaM_Birth <- kappaM_seasonalApril #April - September at birth
pars_pT_imm$kappaM_BirthPT <- kappaM_seasonalApril #April - September at birth

pars_pT_imm$kappaM_catchUp <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$kappaM_catchUpPT <- kappaM_catchUp #only for non-newborn infants
pars_pT_imm$catchAge <- age_catchUp_2seasons

#Easiest to turn off second dose and just increase age eligibility for dose above
pars_pT_imm$kappaM_dose2 <-  kappaM_zero
pars_pT_imm$kappaM_dose2PT <- kappaM_zero
pars_pT_imm$dose2Age <- 0

deSolve_out <- deSolve::ode(y0, times, deSolve_risk_imm, pars_pT_imm)

modelOut_Risk_seasonalAtBirthCU_2ndALL <- aggregate_output_immiPT (deSolveOut = deSolve_out,
                                                                   times = times,
                                                                   nAges = age$nAges,
                                                                   nState = numImm_Risk,
                                                                   nRet = numImm_Risk + numOutB_Risk,
                                                                   old = 1)
###########################################################################################
#CREATE COMPARATIVE PLOTS 

#Figure 3
pMonths <- 12
maxAge <- 24

noVaxx_PT <- modelOut_Risk_NoVaxx$deSolve[(run_in + 1):(run_in + 12), 1:60, 43]
WAVaxx_PT <- modelOut_Risk_seasonalAtBirth85CU65_2nd25$deSolve[(run_in + 1):(run_in + 12), 1:60, 43]
WAVaxx_100 <- modelOut_Risk_seasonalAtBirthCU_2nd$deSolve[(run_in + 1):(run_in + 12), 1:60, 43]

hosp <- noVaxx_PT[1:pMonths,1:maxAge]
yMax <- max(apply(hosp, 1, sum)) #total averted per month over all ages
melt_hosp <- reshape2::melt(hosp)
colnames(melt_hosp) <- c("Month", "Age", "Hospitalisations")

hosp_total_noVAxx <- apply(hosp, 1, sum)
hosp_total_noVAxx <- as.data.frame(cbind(1:12, 1, hosp_total_noVAxx ))
colnames(hosp_total_noVAxx) <- c("Month", "Age", "Hospitalisations")

xlabel <- "Month"
ylabel <- "Hospitalisations <2yo"

p1 <- ggplot2::ggplot(data=melt_hosp, ggplot2::aes(x = Month, y = Hospitalisations, fill = as.numeric(Age) )) +
  ggplot2::geom_bar(position="stack", stat="identity")+
  ggplot2::scale_x_continuous(name = "Month", breaks = 1:12, 
                              labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, plyr::round_any(yMax, 5, ceiling)),
                              breaks = seq(0, plyr::round_any(yMax, 5, ceiling),5))+
  ggplot2::scale_fill_gradientn(name = "Age (months)", colours = c(viridisLite::viridis(13), rep(viridisLite::viridis(13)[13], maxAge - 13)))+
  theme(legend.title=element_text(size=1))+
  ggplot2::geom_line(data = hosp_total_noVAxx, aes(x = Month, y = Hospitalisations), lwd = 1)+
  ggplot2::theme_bw()
p1

hosp_WA <- WAVaxx_PT[1:pMonths,1:maxAge]
melt_hosp_WA <- reshape2::melt(hosp_WA)
colnames(melt_hosp_WA) <- c("Month", "Age", "Hospitalisations")

hosp_total_WA <- apply(hosp_WA, 1, sum)
hosp_total_WA <- as.data.frame(cbind(1:12, 1, hosp_total_WA ))
colnames(hosp_total_WA) <- c("Month", "Age", "Hospitalisations")


p2 <- ggplot2::ggplot(data=melt_hosp_WA, ggplot2::aes(x = Month, y = Hospitalisations, fill = as.numeric(Age) )) +
  ggplot2::geom_bar(position="stack", stat="identity")+
  ggplot2::scale_x_continuous(name = "Month", breaks = 1:12, 
                              labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, plyr::round_any(yMax, 5, ceiling)),
                              breaks = seq(0, plyr::round_any(yMax, 5, ceiling),5))+
  ggplot2::scale_fill_gradientn(name = "Age (months)", colours = c(viridisLite::viridis(13), rep(viridisLite::viridis(13)[13], maxAge - 13)))+
  theme(legend.title=element_text(size=1))+
  ggplot2::geom_line(data = hosp_total_WA, aes(x = Month, y = Hospitalisations), lwd = 1, lty = "dashed")+
  ggplot2::geom_line(data = hosp_total_noVAxx, aes(x = Month, y = Hospitalisations), lwd = 1)+
  ggplot2::theme_bw()
p2

#Figure 4
maxAgeProp <- 24 #Proportion of hospitalisations under 2yo

noVaxx_PT <- modelOut_Risk_NoVaxx$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirth <- modelOut_Risk_seasonalAtBirth100$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirth85 <- modelOut_Risk_seasonalAtBirth85$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirthApril <- modelOut_Risk_seasonalAtBirthApril100$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirthApril85 <- modelOut_Risk_seasonalAtBirthApril85$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirthCU <- modelOut_Risk_seasonalAtBirthCU$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirth85CU65 <- modelOut_Risk_seasonalAtBirth85CU65$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirthCU_2nd <- modelOut_Risk_seasonalAtBirthCU_2nd$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirth85CU65_2nd25 <- modelOut_Risk_seasonalAtBirth85CU65_2nd25$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirth85CU65_2nd25ALL <- modelOut_Risk_seasonalAtBirth85CU65_2nd25ALL$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]
seasonalAtBirthCU_2ndALL <- modelOut_Risk_seasonalAtBirthCU_2ndALL$deSolve[(run_in + 1):(run_in + 12), 1:maxAgeProp, 43]

totalHospPerYear <- sum(noVaxx_PT)

avertSeasonalAtBirth <- noVaxx_PT - seasonalAtBirth
avertSeasonalAtBirth[which(avertSeasonalAtBirth < 0)] <- 0
avertSeasonalAtBirth <- apply(avertSeasonalAtBirth, 2, sum)/totalHospPerYear

avertSeasonalAtBirth85 <- noVaxx_PT - seasonalAtBirth85
avertSeasonalAtBirth85[which(avertSeasonalAtBirth85 < 0)] <- 0
avertSeasonalAtBirth85 <- apply(avertSeasonalAtBirth85, 2, sum)/totalHospPerYear

avertSeasonalAtBirthApril <- noVaxx_PT - seasonalAtBirthApril
avertSeasonalAtBirthApril[which(avertSeasonalAtBirthApril < 0)] <- 0
avertSeasonalAtBirthApril <- apply(avertSeasonalAtBirthApril, 2, sum)/totalHospPerYear

avertSeasonalAtBirthApril85 <- noVaxx_PT - seasonalAtBirthApril85
avertSeasonalAtBirthApril85[which(avertSeasonalAtBirthApril85 < 0)] <- 0
avertSeasonalAtBirthApril85 <- apply(avertSeasonalAtBirthApril85, 2, sum)/totalHospPerYear

avertSeasonalAtBirthCU <- noVaxx_PT - seasonalAtBirthCU
avertSeasonalAtBirthCU[which(avertSeasonalAtBirthCU < 0)] <- 0
avertSeasonalAtBirthCU <- apply(avertSeasonalAtBirthCU, 2, sum)/totalHospPerYear

avertSeasonalAtBirth85CU65 <- noVaxx_PT - seasonalAtBirth85CU65
avertSeasonalAtBirth85CU65[which(avertSeasonalAtBirth85CU65 < 0)] <- 0
avertSeasonalAtBirth85CU65 <- apply(avertSeasonalAtBirth85CU65, 2, sum)/totalHospPerYear

avertSeasonalAtBirthCU_2nd <- noVaxx_PT - seasonalAtBirthCU_2nd
avertSeasonalAtBirthCU_2nd[which(avertSeasonalAtBirthCU_2nd < 0)] <- 0
avertSeasonalAtBirthCU_2nd <- apply(avertSeasonalAtBirthCU_2nd, 2, sum)/totalHospPerYear

avertSeasonalAtBirth85CU65_2nd25 <- noVaxx_PT - seasonalAtBirth85CU65_2nd25
avertSeasonalAtBirth85CU65_2nd25[which(avertSeasonalAtBirth85CU65_2nd25 < 0)] <- 0
avertSeasonalAtBirth85CU65_2nd25 <- apply(avertSeasonalAtBirth85CU65_2nd25, 2, sum)/totalHospPerYear

avertSeasonalAtBirthCU_2ndALL <- noVaxx_PT - seasonalAtBirthCU_2ndALL
avertSeasonalAtBirthCU_2ndALL[which(avertSeasonalAtBirthCU_2ndALL < 0)] <- 0
avertSeasonalAtBirthCU_2ndALL <- apply(avertSeasonalAtBirthCU_2ndALL, 2, sum)/totalHospPerYear

avertSeasonalAtBirth85CU65_2nd25ALL <- noVaxx_PT - seasonalAtBirth85CU65_2nd25ALL
avertSeasonalAtBirth85CU65_2nd25ALL[which(avertSeasonalAtBirth85CU65_2nd25ALL < 0)] <- 0
avertSeasonalAtBirth85CU65_2nd25ALL <- apply(avertSeasonalAtBirth85CU65_2nd25ALL, 2, sum)/totalHospPerYear

maxAge <- 24

plotData <- as.data.frame(rbind(c("SeasonalAtBirth", "WA coverage", sum(avertSeasonalAtBirth85[1:maxAge])),
                                c("SeasonalAtBirth", "100% coverage", sum(avertSeasonalAtBirth[1:maxAge])),
                                c("SeasonalAtBirthApril", "WA coverage", sum(avertSeasonalAtBirthApril85[1:maxAge])),
                                c("SeasonalAtBirthApril", "100% coverage", sum(avertSeasonalAtBirthApril[1:maxAge])),
                                c("SeasonalAtBirthCU", "WA coverage", sum(avertSeasonalAtBirth85CU65[1:maxAge])),
                                c("SeasonalAtBirthCU", "100% coverage", sum(avertSeasonalAtBirthCU[1:maxAge])),
                                c("SeasonalAtBirthCU_2nd", "WA coverage", sum(avertSeasonalAtBirth85CU65_2nd25[1:maxAge])),
                                c("SeasonalAtBirthCU_2nd", "100% coverage", sum(avertSeasonalAtBirthCU_2nd[1:maxAge])),
                                c("SeasonalAtBirthCU_2ndALL", "WA coverage", sum(avertSeasonalAtBirth85CU65_2nd25ALL[1:maxAge])),
                                c("SeasonalAtBirthCU_2ndALL", "100% coverage", sum(avertSeasonalAtBirthCU_2ndALL[1:maxAge]))))

colnames(plotData) <- c("Scenario", "Coverage", "Proportions")
plotData$Proportions <- as.numeric(plotData$Proportions)
plotData$Coverage <- factor(plotData$Coverage, levels = c("WA coverage", "100% coverage"))

ggplot2::ggplot(data=plotData, ggplot2::aes(x = Scenario, y = Proportions, fill = Coverage)) +
  ggplot2::geom_bar(position="dodge", stat="identity", width = 0.9)+
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, 0.55),
                              breaks = seq(0, 0.55, 0.05))+
  ggplot2::xlab("")+
  ggplot2::scale_fill_manual(values= c(viridisLite::viridis(22)[14], viridisLite::viridis(22)[6]))+
  ggplot2::theme_bw()+
  ggplot2::coord_flip()





