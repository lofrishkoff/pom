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
    mu.p.beta ~ dnorm(0, 0.01)

    sigma.p.beta ~ dunif(0, 20)
    tau.p.beta <- 1/(sigma.p.beta* sigma.p.beta)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp] ~ dnorm(mu.p.0,tau.p.0)
      p.beta[sp] ~ dnorm(mu.p.beta,tau.p.beta)
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

    ## random effect of year on occupancy
    sigma.psi.year ~ dunif(0,20)
    tau.psi.year <- 1/(sigma.psi.year*sigma.psi.year)
    for (yr in 1:nyr) {
      psi.year[yr] ~ dnorm(0, tau.psi.year)
    }
    
    ## species-specific intercepts are modeled as fixed effects
    for(sp in 1:nsp) {
      psi.0.sp.pre[sp] ~ dunif(0,1)
      psi.0[sp] <- logit(psi.0.sp.pre[sp])
    }
    
    ## create species-specific slopes with regard to the environment
    ## (phylogeny is incorporated here)
    mu.psi.beta ~ dnorm(0,0.001)
    beta.trait ~ dnorm(0,0.001)
    for(sp in 1:nsp) {
      mu.beta.sp[sp] <- mu.psi.beta + beta.trait*trait[sp]
    }

    ## incorporate phylogenetic covariance structure into
    ## species-specific slope expectations.  lambda scales
    ## the phylogenetic correlation matrix
    lambda ~ dunif(0, 1)
    
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    
    ## convert to precision (or co-precision) matrix - While lambda
    ## above affects the relative weighting on the diagonal in
    ## precision matrix (correlations), sigma.psi.beta (below) 
    ## determines the overall magnitudes of the entries thereby 
    ## converting the correlation matrix to a covariance matrix.
    sigma.psi.beta ~ dunif(0,20)
    tau.beta.sp.mat[1:nsp,1:nsp] <-inverse((sigma.psi.beta^2)*beta.mat[,])
    ## draw species-specific slopes
    psi.beta.sp[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])

    ## *******************************************
    ##          Core Likelihood Function
    ## *******************************************

    for(sp in 1:nsp) {
      for(site in 1:nsite) {
        for(yr in 1:nyr) {        
          
          logit(psi[sp,site,yr]) <-
            psi.0[sp] + ## Species intercept
              psi.beta.sp[sp]*env[site] + ## Species’ response to environment 
                psi.site[site] + ## Random site effect
                  psi.year[yr]  ## Random season effect

          Z[sp,site,yr] ~ dbern(psi[sp,site,yr]) ## Draw occupancy states

          logit(p[sp,site, yr]) <- p.0[sp] + ## Species detection intercept
            p.beta[sp] * env[site] ## Detection responses 
          ## to environment          

          E[sp,site,yr] <- Z[sp,site,yr]*p[sp,site,yr] ## Expected detection 
          ## probability

          for(rep in 1:nrep[sp,site,yr]) {
            X[sp,site,yr,rep] ~ dbern(E[sp,site,yr])
          } # /rep
        } # /yr
      } # /site
    } # /sp
  } # /model.jags


  d$data <- d$data[c('X', 'nsp', 'nsite', 'nyr', 'nrep', 'env',
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
