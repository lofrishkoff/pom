## ************************************************************
## run analysis on simulated data
## ************************************************************
setwd('~/Dropbox/ccb_banding/manuscript_PDOM/pom-github')
rm(list=ls())

library(ape)
library(MASS)
library(rjags)
library(R2jags)

source('simulate_data.R')
source('pom.R')
## ************************************************************
set.seed(12345)

num.sp <- 32

tree.b <- rcoal(num.sp, tip.label=1:num.sp)

dd.model <- make.data.JAGS(tree=tree.b,
                           nsp=num.sp,
                           sigma.psi.env=1,
                           sigma.psi.site=1,
                           sigma.psi.season=1,
                           mu.p.0=-1,
                           sigma.p.0=1,
                           beta.trait=1,
                           lambda=1)

z.init <- dd.model$Z
my.inits <- function() {
  list(Z=z.init)
}

dd <- list(data=dd.model,
           inits=my.inits,
           params=c('p.0',
             'p.env',
             'mu.p.0',
             'mu.p.env',
             'sigma.p.0',
             'sigma.p.env',
             'psi.0',
             'psi.beta.sp',
             'beta.0',
             'sigma.psi.env',
             'sigma.psi.site',
             'sigma.psi.season',
             'beta.trait',
             'lambda'))
  
## Run POM - for real problems ni, nb, and nt should all be increased.
res.lam <- run.R2jags.model(dd, ni=1000, nb=100, nt=2, nc=3)

## Examine All Parameter Estimates
sum.lam <- res.lam$BUGSoutput$summary
cols <- c('mean', '2.5%', '97.5%', 'Rhat', 'n.eff')
round(sum.lam[,cols], 3)

## Examine beta.trait estimate -- Different from 0
sum.lam['beta.trait', cols]

## Extract species' slope parameters
betas <- sum.lam[grep('psi.beta.sp\\[', row.names(sum.lam)),]

## Extract lambda
lamdba <- sum.lam['lambda',]

## Examine simulation matrix to calculate bounded highest posterior
## probability around lambda
sm <- res.lam$BUGSoutput$sims.matrix
source('bounded_hpp.R')
bounded.hpp(sm[,'lambda'])

hist(sm[,'lambda'], col='red')

## sigma.p.env term overlaps with zero, so could remove
## species-specific detection effect along environmental gradient if
## model selection is desired.
bounded.hpp(sm[,'sigma.p.env'])

## Visually check for convergence and staionarity of parameter estimates. Many parameters show autocorrelation, signallying that mixing is poor, and confirming that ni and nt should be increased.
traceplot(res.lam)

