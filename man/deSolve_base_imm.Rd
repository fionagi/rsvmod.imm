% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immunisation.R
\name{deSolve_base_imm}
\alias{deSolve_base_imm}
\title{BASE ODE model with multiple exposures and immunisation strategies
MONOCLONAL Year-round at birth, seasonal at birth, seasonal at birth with
catchup month before peak season
3 monoclonal M_S0 compartments
This function is used with R package deSolve to solve the ODEs}
\usage{
deSolve_base_imm(t, y, parms)
}
\arguments{
\item{t}{a vector of times at which model output is wanted}

\item{y}{a matrix of the initial values of the state variables.
Each column corresponds to a state variable.}

\item{parms}{the required model parameters and their values. For this model,
parms <-  list(b0 = b0,
b1 = b1,
phi = phi,
delta = delta,
gamma0 = gamma0,
gamma1 = gamma1,
nu = nu,
nuM = nuM, #durability of mAb / mat vaxx
omega_vect = omega_vect,
omegaE = omegaE,
kappaM_Birth = kappaM_Birth, #coverage at birth time-dependent vector
kappaM_catchUp = kappaM_catchUp, #coverage of catch-up month time-dependent vector
catchAge = catchAge, #catch-up vaxx covers birth to this age group
A = A,
B = B,
C = C,
D = D,
rho = rho, #efficacy
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
BASE ODE model with multiple exposures and immunisation strategies
MONOCLONAL Year-round at birth, seasonal at birth, seasonal at birth with
catchup month before peak season
3 monoclonal M_S0 compartments
This function is used with R package deSolve to solve the ODEs
}
