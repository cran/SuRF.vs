% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SubsampB.R
\name{Subsample_B}
\alias{Subsample_B}
\title{Subsample_B}
\usage{
Subsample_B(B, data, fold, Alpha, prop, weights, ncores, family)
}
\arguments{
\item{B}{the number of sub-samplings to run (e.g., B=1000)}

\item{data}{the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is in terms of 0/1)}

\item{fold}{fold used in lasso cross validation to select the tuning parameter}

\item{Alpha}{1 for Lasso,0 for ridgeression}

\item{prop}{percentage of samples left out for each subsamping}

\item{family}{The distribution family of the response variable}

\item{weights}{In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)}

\item{ncores}{the number of cores to use for parallel computation}


}
\value{
Class.Err: mis-classification error on the left out ones over B runs. A vector of length B.

Lambda: tuning parameters selected from B runs. It is a vector of length B

BETA: It is a matrix used to save the beta coefficients from all B runs
#' @export
}
\description{
This function is to run sub-sampling B times
}
