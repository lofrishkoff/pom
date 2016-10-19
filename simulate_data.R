## ************************************************************
## Simulate data
## ************************************************************

## Logit function
logit <- function(x) log((x)/(1-x))
## Inverse Logit function
invlogit <- function(x) exp(x)/(exp(x)+1)

## Create a balanced tree with branch lengths
balanced.tree <- function(nsp) {
  tree <- stree(n=nsp, tip.label=1:nsp, type="balanced")
  tree$edge.length <- node.depth(tree)[tree$edge[,2]]
  tree
}

make.data <- function(tree=NULL, ## input a tree
                      nsp=64, ## number of species
                      nsite=15, ## number of sites
                      nseason=5, ## number of seasons
                      nrep=5, ## number of samples within a season
                      mu.psi.0=0, ## mean intercept for occupancy
                      sigma.psi.0=0.5, ## variation in intercepts for
                                       ## occupancy 
                      mu.p.0=0, ## mean interecpt for detectability
                      sigma.p.0=0.5, ## variation in intercepts for
                                     ## detectability
                      mu.p.env=0, ## mean detection response to
                                     ## environment 
                      sigma.p.env=0, ## variance of detection
                                        ## responses to environment
                      sigma.psi.site=0.1, ## variation in site random
                                          ## effect for occupancy
                      sigma.psi.season=0.1, ## variation in season random
                                            ## effect for occupancy
                      sigma.psi.env=1, ## variation in species'
                                       ## responses to the
                                       ## environment.  This scales
                                       ## the width of the covariance
                                       ## matrix that is used to draw
                                       ## species' responses to the
                                       ## environment.
                      beta.0=0, ## mean response to environmental
                                ## gradient 
                      beta.trait=0, ## effect of trait on species'
                                    ## responses to the environment
                      lambda=1, ## weighting on phylogenetic
                                ## covariance matrix (e.g., with
                                ## lambda=1, the covariance matrix is
                                ## given full weight and with
                                ## lambda=0, the model reduces to a
                                ## standard random effect)
                      lambda.trait=0, ## the degree to wich the trait
                                      ## is phylogenetically conserved
                      psi.image=FALSE, ## show plots while creating data
                      occ.image=FALSE,
                      det.image=FALSE,
                      plot.tree=FALSE) {

  ## Create list with all inputs
  inputs <- as.list(environment())
  
  ## Generate phylogenetic relatedness for all speceis and
  ## variance-covariance matrices
  if(is.null(tree)==TRUE)
    tree <- rcoal(nsp, tip.label=1:nsp)
  tree.vcv <- vcv(tree, corr=T)
  i <- order(as.numeric(rownames(tree.vcv)))
  j <- order(as.numeric(colnames(tree.vcv)))
  tree.vcv <- tree.vcv[i,j]
  VCOV.mat       <- tree.vcv*lambda       + (1-lambda)*diag(nsp)
  VCOV.mat.trait <- tree.vcv*lambda.trait + (1-lambda.trait)*diag(nsp)

  ## Generate scaled environmental gradient
  site.env <- seq(-1,1, length.out=nsite)
  site.env <- as.numeric(scale(site.env))

  ## Generate scaled trait values for all species
  spp.trait <- mvrnorm(mu=rep(0, nsp), Sigma=VCOV.mat.trait)
  spp.trait <- as.numeric(scale(spp.trait))
  
  ## Generate species-specific occupancy intercepts
  alpha.spp <- rnorm(nsp, mean=mu.psi.0, sd=sigma.psi.0)

  ## Generate species-specific responses to the environmental
  ## gradient using phylogenetic variance-covariance matrix
  psi.beta.sp.rand <- mvrnorm(mu=rep(0, nsp), 
                              Sigma=(sigma.psi.env^2)*VCOV.mat)
  beta.spp <- beta.0 + beta.trait*spp.trait + psi.beta.sp.rand

  ## Generate random effect of site on occupancy
  site.effect <- rnorm(nsite, mean=0, sd=sigma.psi.site)

  ## Generate random effect of season on occupance
  season.effect <- rnorm(nseason, mean=0, sd=sigma.psi.season)
  
  ## Combine the above components to generate array of expected
  ## occupancies
  ##
  psi.indices <- expand.grid(sp=1:nsp, site=1:nsite, season=1:nseason)
  ##
  alpha    <- alpha.spp[psi.indices$sp]
  beta     <- beta.spp[psi.indices$sp]
  env      <- site.env[psi.indices$site]
  ran.site <- site.effect[psi.indices$site]
  ran.season   <- season.effect[psi.indices$season]
  ## occupancy probability for each species/site/season
  psi <- array(invlogit(alpha + beta*env + ran.site + ran.season),
               dim=c(nsp,nsite,nseason)) 
  
  ## Plot image of mean occupancy probability across seasons if desired
  if(psi.image==TRUE) {
    quartz()
    image(t(apply(psi, c(1,2) , mean)),
          main='Mean Psi across seasons',
          xlab='Sites', ylab='Spp',
          col=rev(gray.colors(12))) 
  }

  ## Simulate true occupancy states.
  ZZ <- array(rbinom(nsp*nsite*nseason, 1, psi), dim=c(nsp,nsite,nseason))

  ## Plot image of mean occupancy states across seasons if desired
  if(occ.image==TRUE) {
    quartz()
    image(t(apply(ZZ, c(1,2) , mean)),
          main='Mean Occupancy States accros seasons',
          xlab='Sites', ylab='Spp',
          col=rev(gray.colors(12))) 
  }

  ## Generate array of expected detectabilities
  ##  
  p.indices <- expand.grid(sp=1:nsp, site=1:nsite)
  ##
  det.beta <- rnorm(nsp, mean=mu.p.env, sd=sigma.p.env)
  det.spp  <- rnorm(nsp, mean=mu.p.0,   sd=sigma.p.0)
  det.pre  <- invlogit(det.spp[p.indices$sp] +
                       det.beta[p.indices$sp] * site.env[p.indices$site])
  det <- array(det.pre, dim=c(nsp,nsite,nseason,nrep))


  ## Generate observation matrix
  EE <- array(ZZ, dim=dim(det)) * det
  XX <- array(rbinom(nsp*nsite*nseason*nrep, 1, prob=EE), dim=dim(det))
  ## Add dimension names
  dimnames(XX) <- list(species=paste('sp.', 1:nsp, sep=''),
                       site=paste('site.', 1:nsite, sep=''),
                       season=1:nseason,
                       rep=1:nrep)

  ## Plot image of mean detection history across seasons if desired
  if(det.image==TRUE) {
    quartz()
    image(t(apply(XX, c(1,2) , mean)),
          main='Mean Detection States across replicates and seasons',
          xlab='Sites', ylab='Spp', col=rev(gray.colors(12))) 
  }
  
  ## Plot tree colored by species' occupancy slopes if desired.
  if(plot.tree==T) {
    quartz()
    colsTips <- colorize(beta.spp,
                         colors=c('red','gold','darkgreen'),
                         min= min(beta.spp), max=max(beta.spp)) 

    phycomp <- psi.beta.sp.rand 
    
    ## Use ancestral state reconstruction to quickly generate 
    ## expectation at nodes							
    suppressWarnings(anres <- ace(phycomp, phy=tree, method='REML', CI=FALSE))
    
    allnodes <- c(phycomp, anres$ace)

    colsEdges <- colorize(allnodes,
                          colors=c('red','gold','darkgreen'),
                          min=min(beta.spp), max=max(beta.spp))

    plot.phylo(tree, type='p', show.tip.label=F,
               edge.color= colsEdges[match(tree$edge[,1], names(allnodes))]) 
    tiplabels(pch=16, col=colsTips[tree$tip.label])

    print('Tree tips colored according to overall species response to environmental gradient (i.e., beta.spp), branches colored according to the phylogenetic component of (i.e., sigma.psi.env * tree.vcv)')
  }

  ## create a long-form data-set for use with lme4 or pglmm in pez
  df <- data.frame(presence=as.vector(XX),
                   expand.grid(spp=1:nsp,
                               site=1:nsite,
                               season=1:nseason,
                               rep=1:nrep))
  df$env <- site.env[df$site] 

  return(c(inputs,
           list(env=site.env,
                trait=spp.trait,
                tree=tree,
                tree.vcv=tree.vcv,
                XX=XX,
                ZZ=ZZ,
                df=df)))
}

## Package up the data and create additional structures so that it is
## ready for use in the JAGS model
make.data.JAGS <- function(...) {

  dd <- make.data(...)

  nrep <- apply(dd$XX, 1:3, function(x) sum(x>=0,na.rm=TRUE))

  Z <- apply(dd$XX, 1:3,
             function(x) (sum(x,na.rm=TRUE)>0)*1)
  Z[apply(dd$XX, 1:3, function(x) !any(!is.na(x)))] <- NA

  list(inputs=dd,
       nsp=dd$nsp,
       nsite=dd$nsite, 
       nseason=dd$nseason,
       nrep=nrep,
       trait=dd$trait,
       env=dd$env,
       X=dd$XX,
       Z=Z,
       ID=diag(1, nrow=dd$nsp, ncol=dd$nsp),
       VCOV=dd$tree.vcv)
}
