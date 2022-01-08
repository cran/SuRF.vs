## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------

set.seed(12345)

library(SuRF.vs)
N=50
p=20
nzc=p/3
X=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=X[,seq(nzc)]%*%beta/3
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

alpha=seq(0.01,0.1,len=5)


 
 #binomial model
 XX=X[,1:2]
 f=1+XX%*%c(2,1.5)
 p=exp(f)/(1+exp(f))
 y=rbinom(100,1,p)
 weights=FALSE
 family=stats::binomial(link="logit")
 surf_binary=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=weights,B=50,C=10,ncores=1,display.progress=TRUE,family=family,alpha_u=0.1,alpha=alpha)
 
 
 #linear regression
 y=1+XX%*%c(0.1,0.2)
 
 family=stats::gaussian(link="identity")
 surf_lm=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=weights,B=100,C=15,ncores=1,display.progress=TRUE,family=family,alpha_u=0.1,alpha=alpha)
 
 
 #cox proportional model
 y=cbind(time=ty,status=1-tcens)
 family=list(family="cox")
 surf_cox=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=FALSE,B=50,C=5,ncores=1,display.progress=TRUE,family=family,alpha_u=alpha_u,alpha=alpha)

