\name{SURF}
\alias{SURF}
\title{SURF}
\usage{
SURF(
  Xo,
  y,
  X = NULL,
  fold = 10,
  Alpha = 1,
  prop = 0.1,
  weights = FALSE,
  B = 1000,
  C = 200,
  ncores = 1,
  display.progress = TRUE,
  family = stats::binomial(link = "logit"),
  alpha_u = 0.1,
  alpha = 0.05
)
}
\arguments{
\item{Xo}{- other type of predictor variables}

\item{y}{- response variable, a vecotr for most families. For family="cox", y will should be a matrix of the response variable in column1 and censoring status in column 2.}

\item{X}{- count data, need to be converted to proportion}

\item{fold}{- number of folds for cross-validation in Lasso}

\item{Alpha}{- Alpha parameter for elastic net}

\item{prop}{- proportion of observations left out in subsampling}

\item{weights}{- use weighted regression: for unbalanced class sizes (bimomial family only) or weighted sample for other families;In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)}

\item{B}{- number of subsamples to take}

\item{C}{- number of permutations used to estimate null distribution}

\item{display.progress}{- whether SuRF should print a message on completion of each}

\item{alpha_u}{- the upper bound of significance level for the permutation test: alpha_u has to be in the range of (0,1). The large of this value, the longer the program will run;}

\item{alpha}{- the alpha value of interest (alpha >0 and must be <=alpha_u). It can be a single value or a vector.If missing, by default it is 0.05.}

\item{ncores}{whether SuRF should compute in parallel: 1 indicates NOT; anything greater will compute in parallel}

\item{family}{The distribution family of the response variable}
}
\value{
Bmod: sub-sampling results

trdata: data frame including both X and y

ranklist: ranking table

modpath: variable selection path (along the alpha range)

selmod: model results at the selected alpha(s)

family: model family used
}
\description{
SuRF is a sparse variable selection method with uses a subsampling
approach an LASSO to rank variables before applying forward selection
using a permutation test. The function is able to give results at a
range of significance levels simultaneously.
}

\details{
SuRF consists of two steps. In the first step, LASSO variable selection
is applied to a large number of subsamples of the data set, to provide a
list of selected variables for each subsample. This list is used to rank
the variables, based on the number of subsamples in which each variable
is selected, so that variables that are selected in more subsamples are
ranked more highly. In the second step, this list is used as a basis for
forward selection, with variables higher on the list tried first. If a
highly-ranked variable is not selected, later variables are tried, and
after each variable is selected, the variables not yet selected (even
previously non-selected variables) are tried in order of the ranking from
Step 1. The decision whether to include a variable is based on a
permutation test for the deviance statistic.

Full details of the SuRF method are in the paper:

  Lihui Liu, Hong Gu, Johan Van Limbergen, Toby Kenney (2020) SuRF: A new
method for sparse variable selection, with application in microbiome
data analysis  Statistics in Medicine 40 897-919

doi: https://onlinelibrary.wiley.com/doi/10.1002/sim.8809
}
\examples{
library(survival)
library(glmnet)
library(SuRF.vs)
N=100;p=200
nzc=p/3
X=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=X[,seq(nzc)]\%*\%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator (1 or 0)
Xo=NULL
B=20
Alpha=1
fold=5
ncores=1
prop=0.1
C=3
alpha_u=0.2
alpha=seq(0.01,0.1,len=20)

#binomial model
XX=X[,1:2]
f=1+XX\%*\%c(2,1.5)
p=exp(f)/(1+exp(f))
y=rbinom(100,1,p)
weights=FALSE
family=stats::binomial(link="logit")

\donttest{
surf_binary=SURF(Xo=X,y=y,fold=5,weights=weights,B=10,C=5,family=family,alpha_u=0.1,alpha=alpha)
}

#linear regression
y=1+XX\%*\%c(0.1,0.2)
family=stats::gaussian(link="identity")
\donttest{
surf_lm=SURF(Xo=X,y=y,fold=5,weights=weights,B=10,C=5,family=family,alpha_u=0.1,alpha=alpha)
}

#cox proportional model
y=cbind(time=ty,status=1-tcens)
weights=rep(1,100)
rseed=floor(runif(20,1,100))
weights[rseed]=2
family=list(family="cox")
\donttest{
surf_cox=SURF(Xo=X,y=y,fold=5,weights=weights,B=10,C=5,family=family,alpha_u=alpha_u,alpha=alpha)
}
}
\keyword{LASSO}
\keyword{forward}
\keyword{selection}
\keyword{variable}
