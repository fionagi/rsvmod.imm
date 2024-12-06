% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immunisation.R
\name{deSolve_risk_imm}
\alias{deSolve_risk_imm}
\title{PRETERM IMMUNISATION ODE model with multiple exposures
MONOCLONAL Very high-risk, Very high-risk 2 dose, Year-around (all) + extra
dose preterm
3 monoclonal M_S0 compartments}
\usage{
deSolve_risk_imm(t, y, parms)
}
\arguments{
\item{t}{a vector of times at which model output is wanted}

\item{y}{a matrix of the initial values of the state variables.
Each column corresponds to a state variable.}

\item{parms}{the required model parameters and their values. For this model,
parms <- list(b0 = b0,
b1 = b1,
phi = phi,
delta = delta,
gamma0 = gamma0,
gamma1 = gamma1,
nu = nu,
nuM = nuM, #durability of mAb
omega_vect = omega_vect,
omegaE = omegaE,
kappaM_Birth = kappaM_Birth, #coverage at birth for terms
kappaM_BirthPT = kappaM_BirthPT, #coverage at birth for preterms
kappaM_catchUp = kappaM_catchUp, #coverage and month of catch-up for terms
KappaM_catchUpPT = kappaM_catchUpPT, #coverage and month of catch-up for preterms
catchAge = catchAge, #catch-up vaxx covers birth to this age group for all infants
alpha = alpha, #proportion preterm
A = A,
B = B,
C = C,
D = D,
E = E,
rho = rho, #efficacy for terms
rho_PT = rho_PT, #efficacy for preterms
age_in = age_in,
age_out = age_out,
age_months = age_months,
sigma_vect = sigma_vect,
sigmaE = sigmaE,
mixing = mixing,
nAges = nAges)}
}
\value{
list
}
\description{
This function is used with R package deSolve to solve the ODEs
}
