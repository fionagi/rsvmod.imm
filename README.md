# rsvmod.imm
A mathematical compartmental transmission model for RSV with nirsevimab immunisation of infants

Based on foundational work described here https://doi.org/10.1186/s12879-024-09400-2, with code here https://github.com/fionagi/rsvmod, 
we fit a dynamic compartmental model of RSV transmission to contemporary data (2015-2019) from temperate Western Australia (WA) and model the administration of nirsevimab at birth and/or in childhood, for both at-risk infants (i.e., preterm infants) and the wider infant population. The code here was used to simulate the 2024 WA RSV immunisation program, predicting impact on RSV-hospitalisations under the estimated coverage levels in 2024. We also modelled the potential impact of higher nirsevimab coverage for different age and risk groups, and identified the components of the program that have the most impact on reducing RSV-hospitalisations. The wrapper file, **main_imm.r** provides examples of how to run the model and produce some of the key plots.

A basic summary of the R code files follows. The functions are described in further detail in the comments above the relevant code.

## data.r ##
This file contains descriptions of the data files used. Note that the timeseries of monthly WA RSV-hospitalisations used in model fitting is not provided due to confidentiality issues.

## parameters.r ##
Key parameter values are set here, including the estimated fitted values of the seasonal forcing function and exponential function that relates infections (by age) to hospitalisations (the age-to-risk function).

## models.r ##
Contains the functions that define the model ODEs without immunisation. These functions are used by the R package _deSolve_ to solve the ODEs. The deSolve_base model does not stratify the population by gestational age whereas deSolve_preTerm stratifies the population into preterm (born <37 weeks gestational age) and full-term.

## immunisation.r ##
Contains the functions that define the model ODEs with immunisation. These functions are used by the R package _deSolve_ to solve the ODEs. The deSolve_base_imm model does not stratify the population by gestational age whereas deSolve_preTerm_imm stratifies the population into preterm (born <37 weeks gestational age) and full-term. 


## functions.r ##
This file contains functions to set up the age structure used in the model,and functions used to set up the likelihood functions for model fitting using an mcmc process as implemented in the R package _lazymcmc_.

## diagnostics.r ##
This file contains functions to produce a range of different plots.

## Linked publications ##

If you wish to use any part of this code, please reference

Giannini F, Hogan AB, Cameron E, Le H, Minney-Smith C, Richmond P, Blyth CC, Glass K, Moore HC; STAMP RSV investigator team. Estimating the impact of Western Australia's first respiratory syncytial virus immunisation program for all infants: A mathematical modelling study. Vaccine. 2025 May 22;56:127155. doi: 10.1016/j.vaccine.2025.127155. Epub 2025 May 7. PMID: 40339485.

