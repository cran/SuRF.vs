% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.dev.R
\name{update_dev}
\alias{update_dev}
\title{update_dev}
\usage{
update_dev(data, vslist, C, weights, ncores, family)
}
\arguments{
\item{data}{the variable 'data' within seqcutoff()}

\item{vslist}{a vector of selected variables}

\item{C}{the number of permutation times}

\item{family}{family=stats::gaussian(link="identity"));family=stats::binomial(link="logit");family=list(family="cox");etc.}

\item{weights}{In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)}

\item{ncores}{the number of cores to use for parallel computation}
}
\value{
dev: a vector of deviance after C permutations (length OF this vector is C)
}
\description{
This function is to derive the deviance distribution based on the permutation method
This function is not to be used independently but will be called by the function selectnew()
}
