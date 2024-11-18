#Set parameter values

#Age group structure
final_age <- 80
numMonthly <- 5*12 #number of age groups (from age 0), in 1 month groups
                   #Assume the first 5 years are in monthly groups

#Fixed values
omega <- 1 #reduced infectiousness in older age groups
omegaE <- 0.7 #reduced infectiousness from prior exposure
nu <- 0.132 #immunity period 230 days 1/(230/(365.25/12))
delta <- 7.610 #latent period 4 days 1/(4/(365.25/12))
gamma0 <- 3.044 #infectious period after first exposure 10 days 1/(10/(365.25/12))
gamma1 <- 4.348 #infectious period after subsequent exposures 7 days 1/(7/(365.25/12))

sigma <- as.vector(c(0.6674784, 0.7314309, 0.7830837, 0.8248023,
            0.8584973, 0.8857119, 0.9076925, 0.9254456,
            0.9397843, 0.9513653, 0.9607190, 0.9682738))  #maternal protection for first 12 months
                                                                  #assuming only 37% of mothers pass on protection,
                                                                  #scaled susceptibility

sigmaE <- 0.77 #reduced susceptibility from prior exposure
alpha <- 0.092  #proportion of births born pre-term (<37 weeks)

#Fitted values - for 95% infected by age 2y.o.
b0_fit <- 0.049356195
b1_fit <- 0.143672161
phi_fit <- -pi

#Risk of Hosp. = (Ae^(-Bi)+C)D*E
A_fit <- 0.025995568 #max. risk of hospitalisation (at age 0)
B_fit <- 0.185155023 #decay constant
C_fit <- 0.007098669 #min. risk across all ages
D_fit <- 0.2 #scaling factor for previous hospitalisations
E_fit <- 2.548210491 #scaling of risk for pre-terms

#Immunisation parameters
dur <- 150 #average durability of mAb is 150 days
nuM <- 0.60875 #waning rate of mAb in each of 3 compartments 1/((150/3)/(365.25/12))
kappaM <- 1 #coverage of mAb
rho <- 1 # reduction of hospitalisation risk of those protected by mAb,
         # taking into account Erlang-3 distribution waning
rho_PT <-1 # reduction of hospitalisation risk of those protected by mAb,
           # taking into account Erlang-3 distribution waning
