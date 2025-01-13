#'BASE ODE model with multiple exposures and immunisation strategies
#'MONOCLONAL Year-round at birth, seasonal at birth, seasonal at birth with
#'catchup month before peak season
#'3 monoclonal M_S0 compartments
#'This function is used with R package deSolve to solve the ODEs
#'
#' @param t a vector of times at which model output is wanted
#' @param y a matrix of the initial values of the state variables.
#'          Each column corresponds to a state variable.
#' @param parms the required model parameters and their values. For this model,
#'                        parms <-  list(b0 = b0,
#'                                       b1 = b1,
#'                                       phi = phi,
#'                                       delta = delta,
#'                                       gamma0 = gamma0,
#'                                       gamma1 = gamma1,
#'                                       nu = nu,
#'                                       nuM = nuM, #durability of mAb / mat vaxx
#'                                       omega_vect = omega_vect,
#'                                       omegaE = omegaE,
#'                                       kappaM_Birth = kappaM_Birth, #coverage at birth time-dependent vector
#'                                       kappaM_catchUp = kappaM_catchUp, #coverage of catch-up month time-dependent vector
#'                                       catchAge = catchAge, #catch-up vaxx covers birth to this age group
#'                                       A = A,
#'                                       B = B,
#'                                       C = C,
#'                                       D = D,
#'                                       rho = rho, #efficacy
#'                                       age_in = age_in,
#'                                       age_out = age_out,
#'                                       age_months = age_months,
#'                                       sigma_vect = sigma_vect,
#'                                       sigmaE = sigmaE,
#'                                       mixing = mixing,
#'                                       nAges = nAges)
#' @return list
#' @export
deSolve_base_imm <- function(t, y, parms) {

  y <- matrix(y, nrow = 75, ncol = 38)


  with(as.list(c(y, parms)), {

    #INITIALISE
    #unprotected
    S0 <- y[, 1]
    E0 <- y[, 2]
    I0 <- y[, 3]
    R0 <- y[, 4]
    S1 <- y[, 5]
    E1 <- y[, 6]
    I1 <- y[, 7]
    R1 <- y[, 8]
    #monoclonal vaxx
    M_S0_1 <- y[, 9]
    M_S0_2 <- y[, 10]
    M_S0_3 <- y[, 11]
    M_E0 <- y[, 12]
    M_I0 <- y[, 13]
    M_R0 <- y[, 14]
    M_S1_1 <- y[, 15]
    M_S1_2 <- y[, 16]
    M_S1_3 <- y[, 17]
    M_E1 <- y[, 18]
    M_I1 <- y[, 19]
    M_R1 <- y[, 20]

    #total population
    N <- S0 + E0 + I0 + R0 + S1 + E1 + I1 + R1 +
      M_S0_1 + M_S0_2 + M_S0_3 + M_E0 + M_I0 + M_R0 +
      M_S1_1 + M_S1_2 + M_S1_3 + M_E1 + M_I1 + M_R1
    #age rate
    tau_shift <- c(age_in[1], age_in[-1])
    tau <- c(age_out[-nAges], age_out[nAges])
    
    #Error catching
    if(dose2Age < catchAge) dose2Age <- catchAge

    #Force of infection
    temp <- omega_vect * (I0 + M_I0 + omegaE*(I1 + M_I1)) / N
    s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
    lambda <- b0 * (1 + b1 * cos(2* pi *t / 12 + phi)) * rowSums(s)

    #NEWLY INFECTED
    #unprotected
    infect0 <- lambda * sigma_vect * S0
    infect1 <- lambda * sigma_vect * sigmaE * S1
    #monoclonal or maternal vaxx
    infectM0_1 <- lambda * sigma_vect * M_S0_1
    infectM0_2 <- lambda * sigma_vect * M_S0_2
    infectM0_3 <- lambda * sigma_vect * M_S0_3
    infectM0 <- infectM0_1 + infectM0_2 + infectM0_3
    infectM1_1 <- lambda *sigma_vect * sigmaE * M_S1_1
    infectM1_2 <- lambda *sigma_vect * sigmaE * M_S1_2
    infectM1_3 <- lambda *sigma_vect * sigmaE * M_S1_3
    infectM1 <- infectM1_1 + infectM1_2 + infectM1_3

    #AGEING
    #Newborns split between unprotected and monoclonal depending on coverage kappaM
    #unprotected
    S0_shift <- c(1 - kappaM_Birth[trunc(t)+1],
                  c(rep(1 - kappaM_catchUp[trunc(t)+1], catchAge), rep(1 - kappaM_dose2[trunc(t)+1], dose2Age - catchAge), rep(1, nAges - dose2Age - 1)) * S0[-nAges])
    E0_shift <- c(0, E0[-nAges])
    I0_shift <- c(0, I0[-nAges])
    R0_shift <- c(0, R0[-nAges])
    S1_shift <- c(0,
                  c(rep(1 - kappaM_catchUp[trunc(t)+1], catchAge), rep(1 - kappaM_dose2[trunc(t)+1], dose2Age - catchAge), rep(1, nAges - dose2Age - 1)) * S1[-nAges])
    E1_shift <- c(0, E1[-nAges])
    I1_shift <- c(0, I1[-nAges])
    R1_shift <- c(0, R1[-nAges])
    #monoclonal
    M_S0_1_shift <- c(kappaM_Birth[trunc(t)+1],
                      c(rep(kappaM_catchUp[trunc(t)+1], catchAge), rep(kappaM_dose2[trunc(t)+1], dose2Age - catchAge), rep(0, nAges - dose2Age - 1)) * S0[-nAges]
                      + M_S0_1[-nAges])
    M_S0_2_shift <- c(0, M_S0_2[-nAges])
    M_S0_3_shift <- c(0, M_S0_3[-nAges])
    M_E0_shift <- c(0, M_E0[-nAges])
    M_I0_shift <- c(0, M_I0[-nAges])
    M_R0_shift <- c(0, M_R0[-nAges])
    M_S1_1_shift <- c(0,
                    c(rep(kappaM_catchUp[trunc(t)+1], catchAge), rep(kappaM_dose2[trunc(t)+1], dose2Age - catchAge), rep(0, nAges - dose2Age - 1)) * S1[-nAges]
                    + M_S1_1[-nAges])
    M_S1_2_shift <- c(0, M_S1_2[-nAges])
    M_S1_3_shift <- c(0, M_S1_3[-nAges])
    M_E1_shift <- c(0, M_E1[-nAges])
    M_I1_shift <- c(0, M_I1[-nAges])
    M_R1_shift <- c(0, M_R1[-nAges])

    #ODEs
    #unprotected
    dS0 <- tau_shift * S0_shift - infect0 + nuM * M_S0_3 - tau * S0
    dE0 <- tau_shift * E0_shift + infect0 - delta * E0 - tau * E0
    dI0 <- tau_shift * I0_shift + delta * E0 - gamma0 * I0 - tau * I0
    dR0 <- tau_shift * R0_shift + gamma0 * I0 - nu * R0 - tau * R0
    dS1 <- tau_shift * S1_shift - infect1 + nu * (R1 + R0 + M_R0 + M_R1) + nuM * M_S1_3 - tau * S1
    dE1 <- tau_shift * E1_shift + infect1 - delta * E1 - tau * E1
    dI1 <- tau_shift * I1_shift + delta*E1 - gamma1 * I1 - tau * I1
    dR1 <- tau_shift * R1_shift + gamma1 * I1 - nu * R1 - tau * R1
    #monoclonal
    dM_S0_1 <- tau_shift * M_S0_1_shift - infectM0_1 - nuM * M_S0_1 - tau * M_S0_1
    dM_S0_2 <- tau_shift * M_S0_2_shift - infectM0_2 + nuM * (M_S0_1 - M_S0_2) - tau * M_S0_2
    dM_S0_3 <- tau_shift * M_S0_3_shift - infectM0_3 + nuM * (M_S0_2 - M_S0_3) - tau * M_S0_3
    dM_E0 <- tau_shift * M_E0_shift + infectM0 - delta * M_E0 - tau * M_E0
    dM_I0 <- tau_shift * M_I0_shift + delta*M_E0 - gamma0 * M_I0 - tau * M_I0
    dM_R0 <- tau_shift * M_R0_shift + gamma0 * M_I0 - nu * M_R0  - tau * M_R0
    dM_S1_1 <- tau_shift * M_S1_1_shift - infectM1_1 - nuM * M_S1_1 - tau * M_S1_1
    dM_S1_2 <- tau_shift * M_S1_2_shift - infectM1_2 + nuM * (M_S1_1 - M_S1_2) - tau * M_S1_2
    dM_S1_3 <- tau_shift * M_S1_3_shift - infectM1_3 + nuM * (M_S1_2 - M_S1_3) - tau * M_S1_3
    dM_E1 <- tau_shift * M_E1_shift + infectM1 - delta * M_E1  - tau * M_E1
    dM_I1 <- tau_shift * M_I1_shift + delta * M_E1 - gamma1 * M_I1 - tau * M_I1
    dM_R1 <- tau_shift * M_R1_shift + gamma1 * M_I1 - nu * M_R1  - tau * M_R1

    #OUTPUT (dinc - incidence / no. of infections, ddetinc - detected incidence / hospitalisations)
    #Unprotected
    dinc0 <- infect0 #First exposure
    ddetinc0 <- (A * exp(-B * age_months) + C) * infect0
    dinc1 <- infect1 #Second exposure
    ddetinc1 <- (A * exp(-B * age_months) + C) * D * infect1
    dinc <- infect0 + infect1 #Combined
    ddetinc <- ddetinc0 + ddetinc1
    #Monoclonal
    dincM0 <- infectM0 #First exposure
    ddetincM0 <- (A * exp(-B * age_months) + C) * (1 - rho) * infectM0
    dincM1 <- infectM1 #Second exposure
    ddetincM1 <- (A * exp(-B * age_months) + C) * (1 - rho) * D * infectM1
    dincM <- infectM0 + infectM1 #Combined
    ddetincM <- ddetincM0 + ddetincM1
    #ALL
    dinc0_comb <- dinc0 + dincM0 #First exposure
    ddetinc0_comb <- ddetinc0 + ddetincM0
    dinc1_comb <- dinc1 + dincM1 #Second exposure
    ddetinc1_comb <- ddetinc1 + ddetincM1
    dinc_comb <- dinc0_comb + dinc1_comb #Combined
    ddetinc_comb <- ddetinc0_comb + ddetinc1_comb

    return(list(c(dS0, #unprotected
                  dE0,
                  dI0,
                  dR0,
                  dS1,
                  dE1,
                  dI1,
                  dR1,
                  dM_S0_1, #monoclonal
                  dM_S0_2,
                  dM_S0_3,
                  dM_E0,
                  dM_I0,
                  dM_R0,
                  dM_S1_1,
                  dM_S1_2,
                  dM_S1_3,
                  dM_E1,
                  dM_I1,
                  dM_R1,
                  dinc0, #unprotected
                  ddetinc0,
                  dinc1,
                  ddetinc1,
                  dinc,
                  ddetinc,
                  dincM0, #monoclonal
                  ddetincM0,
                  dincM1,
                  ddetincM1,
                  dincM,
                  ddetincM,
                  dinc0_comb, #combined
                  ddetinc0_comb,
                  dinc1_comb,
                  ddetinc1_comb,
                  dinc_comb,
                  ddetinc_comb
    )))
  })
}

#'PRETERM IMMUNISATION ODE model with multiple exposures
#'MONOCLONAL Very high-risk, Very high-risk 2 dose, Year-around (all) + extra
#'dose preterm
#'3 monoclonal M_S0 compartments
#'
#'This function is used with R package deSolve to solve the ODEs
#'
#' @param t a vector of times at which model output is wanted
#' @param y a matrix of the initial values of the state variables.
#'          Each column corresponds to a state variable.
#' @param parms the required model parameters and their values. For this model,
#'              parms <- list(b0 = b0,
#'                            b1 = b1,
#'                            phi = phi,
#'                            delta = delta,
#'                            gamma0 = gamma0,
#'                            gamma1 = gamma1,
#'                            nu = nu,
#'                            nuM = nuM, #durability of mAb
#'                            omega_vect = omega_vect,
#'                            omegaE = omegaE,
#'                            kappaM_Birth = kappaM_Birth, #coverage at birth for terms
#'                            kappaM_BirthPT = kappaM_BirthPT, #coverage at birth for preterms
#'                            kappaM_catchUp = kappaM_catchUp, #coverage and month of catch-up for terms
#'                            KappaM_catchUpPT = kappaM_catchUpPT, #coverage and month of catch-up for preterms
#'                            catchAge = catchAge, #catch-up vaxx covers birth to this age group for all infants
#'                            alpha = alpha, #proportion preterm
#'                            A = A,
#'                            B = B,
#'                            C = C,
#'                            D = D,
#'                            E = E,
#'                            rho = rho, #efficacy for terms
#'                            rho_PT = rho_PT, #efficacy for preterms
#'                            age_in = age_in,
#'                            age_out = age_out,
#'                            age_months = age_months,
#'                            sigma_vect = sigma_vect,
#'                            sigmaE = sigmaE,
#'                            mixing = mixing,
#'                            nAges = nAges)
#' @return list
#' @export
deSolve_risk_imm <- function(t, y, parms) {

  y <- matrix(y, nrow = 75, ncol = 46)


  with(as.list(c(y, parms)), {

    #################################################################
    #INITIALISE#
    #UNPROTECTED
    #term population
    S0 <- y[, 1]
    E0 <- y[, 2]
    I0 <- y[, 3]
    R0 <- y[, 4]
    S1 <- y[, 5]
    E1 <- y[, 6]
    I1 <- y[, 7]
    R1 <- y[, 8]
    #pre-term population
    S0_bar <- y[, 9]
    E0_bar <- y[, 10]
    I0_bar <- y[, 11]
    R0_bar <- y[, 12]
    S1_bar <- y[, 13]
    E1_bar <- y[, 14]
    I1_bar <- y[, 15]
    R1_bar <- y[, 16]
    #MONOCLONAL
    #term population
    M_S0_1 <- y[, 17]
    M_S0_2 <- y[, 18]
    M_S0_3 <- y[, 19]
    M_E0 <- y[, 20]
    M_I0 <- y[, 21]
    M_R0 <- y[, 22]
    M_S1_1 <- y[, 23]
    M_S1_2 <- y[, 24]
    M_S1_3 <- y[, 25]
    M_E1 <- y[, 26]
    M_I1 <- y[, 27]
    M_R1 <- y[, 28]
    #preterm population
    M_S0_1_bar <- y[, 29]
    M_S0_2_bar <- y[, 30]
    M_S0_3_bar <- y[, 31]
    M_E0_bar <- y[, 32]
    M_I0_bar <- y[, 33]
    M_R0_bar <- y[, 34]
    M_S1_1_bar <- y[, 35]
    M_S1_2_bar <- y[, 36]
    M_S1_3_bar <- y[, 37]
    M_E1_bar <- y[, 38]
    M_I1_bar <- y[, 39]
    M_R1_bar <- y[, 40]

    #total population
    N <- S0 + E0 + I0 + R0 + S1 + E1 + I1 + R1 +
      S0_bar + E0_bar + I0_bar + R0_bar + S1_bar + E1_bar + I1_bar + R1_bar +
      M_S0_1 + M_S0_2 + M_S0_3 + M_E0 + M_I0 + M_R0 +
      M_S1_1 + M_S1_2 + M_S1_3 + M_E1 + M_I1 + M_R1 +
      M_S0_1_bar + M_S0_2_bar + M_S0_3_bar + M_E0_bar + M_I0_bar + M_R0_bar +
      M_S1_1_bar + M_S1_2_bar + M_S1_3_bar + M_E1_bar + M_I1_bar + M_R1_bar

    ###########################################################################
    #Error catching
    if(dose2Age < catchAge) dose2Age <- catchAge
    
    #age rate
    tau_shift <- c((1 - alpha) * age_in[1], age_in[-1])
    tau <- c(age_out[-nAges], age_out[nAges])
    tauBar_shift <- c(alpha * age_in[1], age_in[-1])
    tauBar <- c(age_out[-nAges], age_out[nAges])

    #Force of infection
    temp <- omega_vect * (I0 + I0_bar + M_I0 + M_I0_bar  +
                            omegaE * (I1 + I1_bar + M_I1 + M_I1_bar)) /N
    s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
    lambda <- b0 * (1 + b1 * cos(2* pi *t / 12 + phi)) * rowSums(s)

    #############################################################

    #NEWLY INFECTED#
    #unprotected
    infect0 <- lambda * sigma_vect * S0
    infect1 <- lambda * sigma_vect * sigmaE * S1
    infect0_bar <- lambda * sigma_vect * S0_bar
    infect1_bar <- lambda * sigma_vect * sigmaE * S1_bar
    #monoclonal
    infectM0_1 <- lambda * sigma_vect * M_S0_1
    infectM0_2 <- lambda * sigma_vect * M_S0_2
    infectM0_3 <- lambda * sigma_vect * M_S0_3
    infectM0 <- infectM0_1 + infectM0_2 + infectM0_3
    infectM1_1 <- lambda *sigma_vect * sigmaE * M_S1_1
    infectM1_2 <- lambda *sigma_vect * sigmaE * M_S1_2
    infectM1_3 <- lambda *sigma_vect * sigmaE * M_S1_3
    infectM1 <- infectM1_1 + infectM1_2 + infectM1_3
    infectM0_1_bar <- lambda * sigma_vect * M_S0_1_bar
    infectM0_2_bar <- lambda * sigma_vect * M_S0_2_bar
    infectM0_3_bar <- lambda * sigma_vect * M_S0_3_bar
    infectM0_bar <- infectM0_1_bar + infectM0_2_bar + infectM0_3_bar
    infectM1_1_bar <- lambda *sigma_vect * sigmaE * M_S1_1_bar
    infectM1_2_bar <- lambda *sigma_vect * sigmaE * M_S1_2_bar
    infectM1_3_bar <- lambda *sigma_vect * sigmaE * M_S1_3_bar
    infectM1_bar <- infectM1_1_bar + infectM1_2_bar + infectM1_3_bar

    ############################################################################

    #AGEING#
    #Newborn preterms split between unprotected and mAb depending on coverage
    #kappaM_PT and terms depending on coverage kappaM

    #unprotected
    S0_shift <-c(1 - kappaM_Birth[trunc(t)+1],
                 c(rep(1 - kappaM_catchUp[trunc(t)+1], catchAge), 
                   rep(1 - kappaM_dose2[trunc(t)+1], 
                       dose2Age - catchAge), 
                   rep(1, nAges - dose2Age - 1)) * S0[-nAges])
    E0_shift <- c(0, E0[-nAges])
    I0_shift <- c(0, I0[-nAges])
    R0_shift <- c(0, R0[-nAges])
    S1_shift <- c(0,
                  c(rep(1 - kappaM_catchUp[trunc(t)+1], catchAge), 
                    rep(1 - kappaM_dose2[trunc(t)+1], dose2Age - catchAge), 
                    rep(1, nAges - dose2Age - 1)) * S1[-nAges])
    E1_shift <- c(0, E1[-nAges])
    I1_shift <- c(0, I1[-nAges])
    R1_shift <- c(0, R1[-nAges])
    S0bar_shift <- c(1 - kappaM_BirthPT[trunc(t)+1],
                     c(rep(1 - kappaM_catchUpPT[trunc(t)+1], catchAge), 
                       rep(1 - kappaM_dose2PT[trunc(t)+1], dose2Age - catchAge), 
                       rep(1, nAges - dose2Age - 1)) * S0_bar[-nAges])
    E0bar_shift <- c(0, E0_bar[-nAges])
    I0bar_shift <- c(0, I0_bar[-nAges])
    R0bar_shift <- c(0, R0_bar[-nAges])
    S1bar_shift <- c(0,
                     c(rep(1 - kappaM_catchUpPT[trunc(t)+1], catchAge), 
                       rep(1 - kappaM_dose2PT[trunc(t)+1], dose2Age - catchAge), 
                       rep(1, nAges - dose2Age - 1)) * S1_bar[-nAges])
    E1bar_shift <- c(0, E1_bar[-nAges])
    I1bar_shift <- c(0, I1_bar[-nAges])
    R1bar_shift <- c(0, R1_bar[-nAges])
    #monoclonal
    M_S0_1_shift <- c(kappaM_Birth[trunc(t)+1],
                      c(rep(kappaM_catchUp[trunc(t)+1], catchAge), 
                        rep(kappaM_dose2[trunc(t)+1], dose2Age - catchAge), 
                        rep(0, nAges - dose2Age - 1)) * S0[-nAges]
                      + M_S0_1[-nAges])
    M_S0_2_shift <- c(0, M_S0_2[-nAges])
    M_S0_3_shift <- c(0, M_S0_3[-nAges])
    M_E0_shift <- c(0, M_E0[-nAges])
    M_I0_shift <- c(0, M_I0[-nAges])
    M_R0_shift <- c(0, M_R0[-nAges])
    M_S1_1_shift <- c(0,
                    c(rep(kappaM_catchUp[trunc(t)+1], catchAge), 
                      rep(kappaM_dose2[trunc(t)+1], dose2Age - catchAge), 
                      rep(0, nAges - dose2Age - 1)) * S1[-nAges]
                    + M_S1_1[-nAges])
    M_S1_2_shift <- c(0, M_S1_2[-nAges])
    M_S1_3_shift <- c(0, M_S1_3[-nAges])
    M_E1_shift <- c(0, M_E1[-nAges])
    M_I1_shift <- c(0, M_I1[-nAges])
    M_R1_shift <- c(0, M_R1[-nAges])
    M_S0bar_1_shift <- c(kappaM_BirthPT[trunc(t)+1],
                         c(rep(kappaM_catchUpPT[trunc(t)+1], catchAge), 
                           rep(kappaM_dose2PT[trunc(t)+1], dose2Age - catchAge), 
                           rep(0, nAges - dose2Age - 1)) * S0_bar[-nAges]
                         + M_S0_1_bar[-nAges])
    M_S0bar_2_shift <- c(0, M_S0_2_bar[-nAges])
    M_S0bar_3_shift <- c(0, M_S0_3_bar[-nAges])
    M_E0bar_shift <- c(0, M_E0_bar[-nAges])
    M_I0bar_shift <- c(0, M_I0_bar[-nAges])
    M_R0bar_shift <- c(0, M_R0_bar[-nAges])
    M_S1bar_1_shift <- c(0,
                       c(rep(kappaM_catchUpPT[trunc(t)+1], catchAge), 
                         rep(kappaM_dose2PT[trunc(t)+1], dose2Age - catchAge), 
                         rep(0, nAges - dose2Age - 1)) * S1_bar[-nAges]
                       + M_S1_1_bar[-nAges])
    M_S1bar_2_shift <- c(0, M_S1_2_bar[-nAges])
    M_S1bar_3_shift <- c(0, M_S1_3_bar[-nAges])
    M_E1bar_shift <- c(0, M_E1_bar[-nAges])
    M_I1bar_shift <- c(0, M_I1_bar[-nAges])
    M_R1bar_shift <- c(0, M_R1_bar[-nAges])

    ############################################################################

    #ODEs
    #unprotected
    dS0 <- tau_shift * S0_shift - infect0 + nuM * M_S0_3 - tau * S0
    dE0 <- tau_shift * E0_shift + infect0 - delta * E0 - tau * E0
    dI0 <- tau_shift * I0_shift + delta * E0 - gamma0 * I0 - tau * I0
    dR0 <- tau_shift * R0_shift + gamma0 * I0 - nu * R0 - tau * R0
    dS1 <- tau_shift * S1_shift - infect1 + nu * (R1 + R0 + M_R0 + M_R1) + nuM * M_S1_3 - tau * S1
    dE1 <- tau_shift * E1_shift + infect1 - delta * E1 - tau * E1
    dI1 <- tau_shift * I1_shift + delta*E1 - gamma1 * I1 - tau * I1
    dR1 <- tau_shift * R1_shift + gamma1 * I1 - nu * R1 - tau * R1

    dS0_bar <- tauBar_shift * S0bar_shift - infect0_bar + nuM * M_S0_3_bar- tauBar * S0_bar
    dE0_bar <- tauBar_shift * E0bar_shift + infect0_bar - delta * E0_bar - tauBar * E0_bar
    dI0_bar <- tauBar_shift * I0bar_shift + delta*E0_bar - gamma0 * I0_bar - tauBar * I0_bar
    dR0_bar <- tauBar_shift * R0bar_shift + gamma0 * I0_bar - nu * R0_bar - tauBar * R0_bar
    dS1_bar <- tauBar_shift * S1bar_shift - infect1_bar + nu * (R1_bar + R0_bar + M_R0_bar + M_R1_bar) + nuM * M_S1_3_bar - tauBar * S1_bar
    dE1_bar <- tauBar_shift * E1bar_shift + infect1_bar - delta * E1_bar - tauBar * E1_bar
    dI1_bar <- tauBar_shift * I1bar_shift + delta*E1_bar - gamma1 * I1_bar - tauBar * I1_bar
    dR1_bar <- tauBar_shift * R1bar_shift + gamma1 * I1_bar - nu * R1_bar - tauBar * R1_bar

    #monoclonal
    dM_S0_1 <- tau_shift * M_S0_1_shift - infectM0_1 - nuM * M_S0_1 - tau * M_S0_1
    dM_S0_2 <- tau_shift * M_S0_2_shift - infectM0_2 + nuM * (M_S0_1 - M_S0_2) - tau * M_S0_2
    dM_S0_3 <- tau_shift * M_S0_3_shift - infectM0_3 + nuM * (M_S0_2 - M_S0_3) - tau * M_S0_3
    dM_E0 <- tau_shift * M_E0_shift + infectM0 - delta * M_E0 - tau * M_E0
    dM_I0 <- tau_shift * M_I0_shift + delta * M_E0 - gamma0 * M_I0 - tau * M_I0
    dM_R0 <- tau_shift * M_R0_shift + gamma0 * M_I0 - nu * M_R0  - tau * M_R0
    dM_S1_1 <- tau_shift * M_S1_1_shift - infectM1_1 - nuM * M_S1_1 - tau * M_S1_1
    dM_S1_2 <- tau_shift * M_S1_2_shift - infectM1_2 + nuM * (M_S1_1 - M_S1_2) - tau * M_S1_2
    dM_S1_3 <- tau_shift * M_S1_3_shift - infectM1_3 + nuM * (M_S1_2 - M_S1_3) - tau * M_S1_3
    dM_E1 <- tau_shift * M_E1_shift + infectM1 - delta * M_E1 - tau * M_E1
    dM_I1 <- tau_shift * M_I1_shift + delta * M_E1 - gamma1 * M_I1 - tau * M_I1
    dM_R1 <- tau_shift * M_R1_shift + gamma1 * M_I1 - nu * M_R1 - tau * M_R1
    
    dM_S0_1_bar <- tauBar_shift * M_S0bar_1_shift - infectM0_1_bar - nuM * M_S0_1_bar - tauBar * M_S0_1_bar
    dM_S0_2_bar <- tauBar_shift * M_S0bar_2_shift - infectM0_2_bar + nuM * (M_S0_1_bar - M_S0_2_bar) - tauBar * M_S0_2_bar
    dM_S0_3_bar <- tauBar_shift * M_S0bar_3_shift - infectM0_3_bar + nuM * (M_S0_2_bar - M_S0_3_bar) - tauBar * M_S0_3_bar
    dM_E0_bar <- tauBar_shift * M_E0bar_shift + infectM0_bar - delta * M_E0_bar - tauBar * M_E0_bar
    dM_I0_bar <- tauBar_shift * M_I0bar_shift + delta * M_E0_bar - gamma0 * M_I0_bar - tauBar * M_I0_bar
    dM_R0_bar <- tauBar_shift * M_R0bar_shift + gamma0 * M_I0_bar - nu * M_R0_bar - tauBar * M_R0_bar
    dM_S1_1_bar <- tauBar_shift * M_S1bar_1_shift - infectM1_1_bar - nuM * M_S1_1_bar - tauBar * M_S1_1_bar
    dM_S1_2_bar <- tauBar_shift * M_S1bar_2_shift - infectM1_2_bar + nuM * (M_S1_1_bar - M_S1_2_bar) - tauBar * M_S1_2_bar
    dM_S1_3_bar <- tauBar_shift * M_S1bar_3_shift - infectM1_3_bar + nuM * (M_S1_2_bar - M_S1_3_bar) - tauBar * M_S1_3_bar
    dM_E1_bar <- tauBar_shift * M_E1bar_shift + infectM1_bar - delta * M_E1_bar - tauBar * M_E1_bar
    dM_I1_bar <- tauBar_shift * M_I1bar_shift + delta*M_E1_bar - gamma1 * M_I1_bar - tauBar * M_I1_bar
    dM_R1_bar <- tauBar_shift * M_R1bar_shift + gamma1 * M_I1_bar - nu * M_R1_bar - tauBar * M_R1_bar
    
    #############################################################################

    #OUTPUT (dinc - incidence, ddetinc - detected incidence / hospitalisations)
    #unprotected
    #Term unprotected
    dinc0 <- infect0
    ddetinc0 <- (A * exp(-B * age_months) + C) * infect0
    dinc1 <- infect1
    ddetinc1 <- (A * exp(-B * age_months) + C) * D * infect1
    dincT <- dinc0 + dinc1
    ddetincT <- ddetinc0 + ddetinc1 #Term hospitalisations unprotected
    #Preterm unprotected
    dinc0_bar <- infect0_bar
    ageRiskF <- (A * exp(-B * age_months) + C) * E
    ageRiskF[which(ageRiskF >1)] <- 1
    ddetinc0_bar <- ageRiskF * infect0_bar
    dinc1_bar <- infect1_bar
    ageRiskS <- (A * exp(-B * age_months) + C) * E * D
    ageRiskS[which(ageRiskS >1)] <- 1
    ddetinc1_bar <- ageRiskS * infect1_bar
    dincPT <- dinc0_bar + dinc1_bar
    ddetincPT <- ddetinc0_bar + ddetinc1_bar #Preterm hospitalisations unprotected
    #Term + Preterm: unprotected
    dinc <- dincT + dincPT
    ddetinc <- ddetincT + ddetincPT #Combined hospitalisations unprotected
    #Monoclonal
    #Term mAb
    dincM0 <- infectM0
    ddetincM0 <- (A * exp(-B * age_months) + C) * (1 - rho) * infectM0
    dincM1 <- infectM1
    ddetincM1 <- (A * exp(-B * age_months) + C) * (1 - rho) * D * infectM1
    dincMT <- dincM0 + dincM1
    ddetincMT <- ddetincM0 + ddetincM1 #Term hospitalisations mAb
    #Preterm mAb
    dincM0_bar <- infectM0_bar
    ageRiskF <- (A * exp(-B * age_months) + C) * (1 - rho_PT) * E
    ageRiskF[which(ageRiskF >1)] <- 1
    ddetincM0_bar <- ageRiskF * infectM0_bar
    dincM1_bar <- infectM1_bar
    ageRiskS <- (A * exp(-B * age_months) + C) * (1 - rho_PT) * E * D
    ageRiskS[which(ageRiskS >1)] <- 1
    ddetincM1_bar <- ageRiskS * infectM1_bar
    dincMPT <- dincM0_bar + dincM1_bar
    ddetincMPT <- ddetincM0_bar + ddetincM1_bar #Preterm hospitalisations mAb
    #Term + Preterm: mAb
    dincM <- dincMT + dincMPT
    ddetincM <- ddetincMT + ddetincMPT #Combined hospitalisations mAb

    #Total hospitalisations for term
    hos_term <- ddetincT + ddetincMT
    #Total hospitalisations for preterm
    hos_preterm <- ddetincPT + ddetincMPT
    #Total hospitalisations combined
    hos_comb <- ddetinc + ddetincM
    #Number of doses term
    doses_term <- c(tau_shift[1] * M_S0_1_shift[1], rep(0, nAges - 1)) #only doses at birth
    #Number of doses preterm
    doses_preterm <- c(tauBar_shift[1] * M_S0bar_1_shift[1], rep(0, nAges - 1)) #only doses at birth
    #Number of doses combined
    doses_comb <- doses_term + doses_preterm

    ############################################################################

    return(list(c(dS0, #unprotected
                  dE0,
                  dI0,
                  dR0,
                  dS1,
                  dE1,
                  dI1,
                  dR1,
                  dS0_bar,
                  dE0_bar,
                  dI0_bar,
                  dR0_bar,
                  dS1_bar,
                  dE1_bar,
                  dI1_bar,
                  dR1_bar,
                  dM_S0_1, #mAb
                  dM_S0_2,
                  dM_S0_3,
                  dM_E0,
                  dM_I0,
                  dM_R0,
                  dM_S1_1,
                  dM_S1_2,
                  dM_S1_3,
                  dM_E1,
                  dM_I1,
                  dM_R1,
                  dM_S0_1_bar,
                  dM_S0_2_bar,
                  dM_S0_3_bar,
                  dM_E0_bar,
                  dM_I0_bar,
                  dM_R0_bar,
                  dM_S1_1_bar,
                  dM_S1_2_bar,
                  dM_S1_3_bar,
                  dM_E1_bar,
                  dM_I1_bar,
                  dM_R1_bar,
                  hos_term, #output
                  hos_preterm,
                  hos_comb,
                  doses_term,
                  doses_preterm,
                  doses_comb
    )))
  })
}
