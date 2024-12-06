% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.R
\name{makePlotAvert}
\alias{makePlotAvert}
\title{Create plots of hospitalisations averted}
\usage{
makePlotAvert(
  noVacc_data,
  vacc_data,
  maxAge = 6,
  pMonths = 24,
  yMax = NA,
  prop_yMax = NA
)
}
\arguments{
\item{noVacc_data}{hospitalisations for pre-vaccine scenario}

\item{vacc_data}{hospitalisations for vaccine scenario}

\item{maxAge}{maximum age (in months) to plot}

\item{pMonths}{no. of months to plot}
}
\value{
ggplot
}
\description{
Create plots of hospitalisations averted
}
