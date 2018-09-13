fourjump <- function(lev,n){
  c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-2*lev, n/5), rep(0,n/5) )
}

## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R
dosim <- function(lev, ichunk, nsim, n=200, meanfun=fourjump, mc.cores=1,
                  numSteps=4, filename=NULL, sigma = 1, sigma.add=0.5, type = "bsfs",
                  outputdir = "../output", locs=1:n) {

  onesim <- function(isim){

    ## Generate data
    set.seed(isim)
    mn = meanfun(lev=lev, n=n)
    y = mn + rnorm(n, 0, sigma)

    y.addnoise = rnorm(n, 0, sigma.add)
    results = list()

    ## Global settings
    max.numSteps = 10
    allsteps = allsteps.plain = 2:max.numSteps

    ## Plain BS inference
    if(any(type=="bsfs")){tryCatch({
      obj = bsfs(y, numSteps=max.numSteps)
      poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
      res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma, locs=locs)
      results$bsfs = res$pvs.by.step
      results$bsfs_zero = res$zeros.by.step
      results$bsfs_cps = obj$cp * obj$cp.sign
    })}

    results.list = mclapply(nsim, onesim, mc.cores=mc.cores)
    return(results.list)
  }

  library(binseginf)
  levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
  nsim = 10000
  nchunk = 10
  locs = unlist(lapply(c(40,80,120,160), function(loc) loc + c(-2, -1, 0, 1 ,2)))
  for(ichunk in 1:nchunk){
    dosim(lev=levs[1], ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=1, locs=locs,
          type="bsfs")
  }
}
