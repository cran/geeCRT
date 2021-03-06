% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleSWCRTLarge.R
\docType{data}
\name{sampleSWCRTLarge}
\alias{sampleSWCRTLarge}
\title{simulated large SW-CRT data}
\format{
A data frame with 1508 rows and 10 variables:
\describe{
  \item{period1}{indicator of being at period 1}
  \item{period2}{indicator of being at period 2}
  \item{period3}{indicator of being at period 3}
  \item{period4}{indicator of being at period 4}
  \item{period5}{indicator of being at period 5}
  \item{treatment}{indicator of being treated}
  \item{id}{cluster identification number}
  \item{period}{period order number}
  \item{y_bin}{binary outcome variable}
  \item{y_con}{continous outcome variable}
}
}
\description{
Simulated cross-sectional individual-level SW-CRT data with 12 clusters and 5 periods. The cluster-period size is uniformly distributed between 20 and 30.
  The correlated binary and continuous outcomes are used for analysis as examples.
}
