% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic.R
\name{plot_riskFit}
\alias{plot_riskFit}
\title{Plot fit to data for risk model (with preterms)}
\usage{
plot_riskFit(
  term,
  pTerm,
  obs,
  years = 2010:2020,
  lineCol = "black",
  lineType = "dashed",
  lineSize = 1,
  plotText = NA,
  fileName = NA
)
}
\arguments{
\item{term}{table of MODEL predicted monthly hospitalisations
for term infants}

\item{pTerm}{table of MODEL predicted monthly hospitalisations
for preterm infants}

\item{obs}{table of OBSERVED term and preterm population monthly
hospitalisation data}

\item{years}{year labels for x-axis}

\item{lineCol}{default line colour is black}

\item{lineType}{default line type is dashed}

\item{lineSize}{default line size is 1}

\item{plotText}{if plotText is NA (as default), plotText is "RSV
hospitalisations observed (dots) and model (dashed)"}

\item{fileName}{if fileName is NA (as default), resulting plots are not
saved. Otherwise fileName is base for 3 resulting jpeg files}
}
\value{
ggplot
}
\description{
Plot fit to data for risk model (with preterms)
}
