run.R2jags.model <- function(d,
                             ni=1100,
                             nb=100,
                             nt=10,
                             nc=3) {


  model.jags<-function(){

    ## *******************************************
    ##          DETECTION COMPONENT
    ## *******************************************

    mu.p.0.pre ~ dunif(0,1)
    mu.p.0 <- logit(mu.p.0.pre)
    
    sigma.p.0 ~ dunif(0,20)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)

    ## Effect of continuous detection covariate
    mu.p.env ~ dnorm(0, 0.01)

    sigma.p.env ~ dunif(0, 20)
    tau.p.env <- 1/(sigma.p.env* sigma.p.env)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      p.env[sp] ~ dnorm(mu.p.env,tau.p.env)
    }
    
    ## *******************************************
    ##          OCCUPANCY COMPONENT
    ## *******************************************

    ## random effect of site on occupancy
    sigma.psi.site ~ dunif(0,20)
    tau.psi.site <- 1/(sigma.psi.site*sigma.psi.site)
    for (site in 1:nsite) {
      psi.site[site] ~ dnorm(0, tau.psi.site)
    }

    ## random effect of season on occupancy
    sigma.psi.season ~ dunif(0,20)
    tau.psi.season <- 1/(sigma.psi.season*sigma.psi.season)
    for(season in 1:nseason) {
      psi.season[season] ~ dnorm(0, tau.psi.season)
    }
    
    ## species-specific intercepts are modeled as fixed effects
    for(sp in 1:nsp) {
      psi.0.sp.pre[sp] ~ dunif(0,1)
      psi.0[sp] <- logit(psi.0.sp.pre[sp])
    }
    
    ## create species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    beta.0 ~ dnorm(0,0.001)
    beta.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      mu.beta.sp[sp] <- beta.0 + beta.trait*trait[sp]
    }

    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  lambda scales
    ## the phylogenetic correlation matrix
    lambda ~ dunif(0, 1)
    
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    
    ## convert to precision (or co-precision) matrix - While lambda
    ## above affects the relative weighting on the diagonal in
    ## precision matrix (correlations), sigma.psi.env (below) 
    ## determines the overall magnitudes of the entries thereby 
    ## converting the correlation matrix to a covariance matrix.
    sigma.psi.env ~ dunif(0,20)
    tau.beta.sp.mat[1:nsp,1:nsp] <-inverse((sigma.psi.env^2)*beta.mat[,])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])

    ## *******************************************
    ##          Core Likelihood Function
    ## *******************************************

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(season in 1:nseason) {        
          
          logit(psi[sp,site,season]) <-
            psi.0[sp] + ## Species intercept
            psi.beta.sp[sp]*env[site] + ## Species’ response to environment 
            psi.site[site] + ## Random site effect
            psi.season[season]  ## Random season effect

          Z[sp,site,season] ~ dbern(psi[sp,site,season]) ## Draw occupancy states

          logit(p[sp,site, season]) <-
            p.0[sp] + ## Species detection intercept
            p.env[sp] * env[site] ## Detection responses 
          ## to environment          

          E[sp,site,season] <- Z[sp,site,season]*p[sp,site,season] ## Expected detection 
          ## probability

          for(rep in 1:nrep[sp,site,season]) {
            X[sp,site,season,rep] ~ dbern(E[sp,site,season])
          } # /rep
        } # /season
      } # /site
    } # /sp
  } # /model.jags


  d$data <- d$data[c('X', 'nsp', 'nsite', 'nseason', 'nrep', 'env',
                     'VCOV', 'ID', 'trait')]
  res <- jags(data=d$data,
              inits=d$inits,
              parameters.to.save=d$params,
              model.file=model.jags,
              n.chains=nc,
              n.thin=nt,
              n.iter=ni,
              n.burnin=nb,
              working.directory=NULL)
  
  res
}
