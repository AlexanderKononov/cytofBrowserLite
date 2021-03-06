% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_preparation.R
\name{exprs_outlier_squeezing}
\alias{exprs_outlier_squeezing}
\title{Detect outliers from exprs data and squeeze from 0 to 1}
\usage{
exprs_outlier_squeezing(exprs_data, use_markers, quantile = 0.01)
}
\arguments{
\item{exprs_data}{matrix with expression data}

\item{quantile}{number from 0 to 1 to set the quantile which is of cutoff threshold}

\item{use_marker}{list of markers for transformation}
}
\value{

}
\description{
Detect outliers from exprs data and squeeze from 0 to 1
}
