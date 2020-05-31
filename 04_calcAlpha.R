###############################################################################	
#
# Script to calculate hypervolumes from fossil pollen samples
#
# Alpha: calculates hypervolumes at sites
# 
# 20170702: Added code for taxa specific selection
###############################################################################	

start.time = Sys.time()

library(hypervolume)
library(dplyr)
library(foreach)
library(parallel)

cl <- parallel::makeCluster(4, setup_timeout = 0.5)
doParallel::registerDoParallel(cl)

randomTraits <- function(n = 1, means, sds) {
  if (length(means) != length(sds)) {stop("Means and SDs not equal")}
  ntraits = length(means)
  trait.out = rep(NA, ntraits)
  
  for (i in 1:ntraits) {
    if (sds[i] > 0 & !is.na(sds[i])) {
      trait.out[i] <- rnorm(n, means[i], sds[i])
    } else {
      trait.out[i] <- means[i]
    }
  }
  return(trait.out)
}

## First get the trait values
load("./caranaPollTraitVals.RData")

## Threshold for taxa selection (0.5%)
tth = 0.005

## Means and standard deviations
trait.mn = c(mean(dat$lsla.mn, na.rm = TRUE), 
             mean(dat$lhgt.mn, na.rm = TRUE), 
             mean(dat$lsdm.mn, na.rm = TRUE))

trait.sd = c(sd(dat$lsla.mn, na.rm = TRUE), 
             sd(dat$lhgt.mn, na.rm = TRUE), 
             sd(dat$lsdm.mn, na.rm = TRUE))

## Missing values
dat = dat %>%
  filter(!is.na(lsla.mn) & !is.na(lhgt.mn) & !is.na(lsdm.mn))

## Only trees
dat = dat %>% 
  filter(grp=="TRSH" | grp=="LIAN" | grp=="DWAR" | grp=="HERB")

sites = unique(sort(dat$ent))
nsites = length(sites)

## Resampling iterations
nbit = 500

## Output
out.hv = rep(NA, nsites)
out.rnd = array(NA, dim = c(nsites, nbit))
out.crds = data.frame(lon=rep(NA, nsites), 
                      lat=rep(NA, nsites), 
                      alt=rep(NA, nsites))
trait.mean = array(NA, dim = c(nsites, 3))
trait.rnd = array(NA, dim = c(nsites, 3, nbit))

## Now loop through by site
for (i in 1:nsites) {
  # for (i in 1:2) {
  print(paste("Doing site", i, sites[i]))
  site.dat = dat %>% 
    filter(ent == sites[i])
  
  ## Loop to remove low values
  nsd = dim(site.dat)[1]
  keepID = NULL

  for (j in 1:nsd) {
    if (site.dat$val[j] > tth) {
      keepID = c(keepID,j)
    }
  }
  # site.dat = site.dat[which(site.dat$val >= tth),]
  site.dat = site.dat[keepID,]
  nsd = dim(site.dat)[1]

  if (nsd > 5) {
    out.crds$lon[i] = site.dat$lon[1]
    out.crds$lat[i] = site.dat$lat[1]
    out.crds$alt[i] = site.dat$alt[1]
    
    trait.dat = data.frame(lsla = site.dat$lsla.mn,
                           lhgt = site.dat$lhgt.mn,
                           lsdm = site.dat$lsdm.mn)
    trait.dat = unique(trait.dat)
    
    samp.hv = hypervolume_box(trait.dat,
                              samples.per.point = 1000, verbose = FALSE)
    
    out.hv[i] = samp.hv@Volume
    trait.mean[i,] <- apply(trait.dat, 2, mean, na.rm = TRUE)
    
    ## Resampling loop - runs in parallel using foreach 
    ## Set up cores above ^
    # out.rnd[i,j,] <- foreach(k=1:nbit, .combine = 'c') %do% {
    oper <- foreach (j=1:nbit, .combine = 'c', .packages='hypervolume') %dopar% {
      ## Get random traits
      trait.dat = data.frame(lsla = randomTraits(means = site.dat$lsla.mn, 
                                                 sds = site.dat$lsla.sd),
                             lhgt = randomTraits(means = site.dat$lhgt.mn, 
                                                 sds = site.dat$lhgt.sd),
                             lsdm = randomTraits(means = site.dat$lsdm.mn, 
                                                 sds = site.dat$lsdm.sd))
      # trait.dat = (trait.dat - trait.mn) / trait.sd
      trait.dat = unique(trait.dat)
      trait.rnd.tmp <- apply(trait.dat, 2, mean, na.rm = TRUE)
      # trait.rnd[i,,j] <- apply(trait.dat, 2, mean, na.rm = TRUE)
      print(j)
      print(apply(trait.dat, 2, mean, na.rm = TRUE))
      print(trait.rnd[i,,j])
      
      # samp.hv = hypervolume_svm(trait.dat,
      #                           samples.per.point = 1000, verbose = FALSE)
      
      # samp.hv = hypervolume_gaussian(trait.dat,
      #                                kde.bandwidth = 0.5,
      #                                quantile.requested = 0.95,
      #                                quantile.requested.type = "probability",
      #                                samples.per.point = 1000, verbose = FALSE)
      # 
      samp.hv = hypervolume_box(trait.dat,
                                samples.per.point = 1000, verbose = FALSE)
      
      # out.rnd[i,j,k] = samp.hv@Volume
      list(samp.hv@Volume, trait.rnd.tmp)
      
    } # Resample loop 
    out.rnd[i,] <- oper[[1]]
    trait.rnd[i,,] <- oper[[2]]
  } # Sample check
  
  
} # Site loop
save(sites, out.hv, out.rnd, out.crds, trait.mean, trait.rnd,
     #nbit=1000
     # file="nampdTraits_alpha_all.RData")
     file="nampdTraits_alpha_all_nbit500.RData")



end.time = Sys.time()

parallel::stopCluster(cl)

