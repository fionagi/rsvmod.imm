# data.R - This file includes documentation of the data used.

#' ABS population of Greater Perth by 5-year age group - 2016
#'
#' Data is formatted to be easily used with the R package conmat
#' to find the applicable contact matrix
#' 
#' Children in Kimberley & Pilbara regions identified via the following postcodes were excluded:
#' - Kimberley: 6725, 6726, 6728, 6731, 6733, 6740, 6743, 6765, 6770
#' - Pilbara: 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6718, 6720, 6721,
#            6722, 6723, 6751, 6753, 6754, 6758, 6760, 6761, 6762
#'
#' @format A data table containing 18 observations of 4 variables.
#' \describe{
#'   \item{lga}{Local Government Area classification}
#'   \item{lower.age.limit}{}
#'   \item{year}{census data year}
#'   \item{population}{population in age group}
#' }
#' @source \url{https://abs.gov.au/census/find-census-data/quickstats/2021/5GPER}
"data.southernWA.2016"

#' Daily contact matrix
#'
#' Daily number of contacts between age groups within the population.
#' This contact matrix was generated using the R package conmat based on the
#' POLYMOD study results for the UK, normalised to 2016 metropolitan Perth
#' population demographics but with finer age-stratification of mixing between 
#' household contacts based on the demographics of the 2016 ABS 5% microdata 
#' sample for the Greater Perth region and between school/childcare contacts of 
#' children under 5 based on the Childhood Education and Care Survey 2017.
#' The <5 year age group is divided into monthly age groups in the model, 
#' with contacts being uniformly distributed into these resulting age groups.
#'
#' @format A 75 x 75 data matrix
#' @source \url{https://rdrr.io/github/njtierney/conmat/}
"data.contactDaily.ABSMicro"
