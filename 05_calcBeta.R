###############################################################################	
#
# Script to calculate hypervolumes from fossil pollen samples
#
# Beta: calculates hypervolumes at samples and the overlap between them
# 
# 20170702: Added code for taxa specific selection
###############################################################################	

start.time = Sys.time()

library(hypervolume)
library(dplyr)
library(foreach)
library(parallel)

# cl <- parallel::makeCluster(4, setup_timeout = 0.5)
# doParallel::registerDoParallel(cl)

## Function to sample from multiple PDFs
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

## Function to fix the output from foreach
comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
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

samps = unique(sort(dat$sample))
nsamps = length(samps)

## Resampling iterations
nbit = 10

## Output
ovlp.hv = rep(NA, nsamps)
ovlp.rnd = array(NA, dim = c(nsamps, nbit))

## Now loop through by site
for (i in 1:(nsamps-1)) {
  # for (i in 1:2) {
  print(paste("Doing sample", i, samps[i]))
  
  ## Which samples?
  age1 = i
  age2 = i + 1
  
  ## Sample 1
  samp1.dat = dat %>% 
    filter(sample == samps[i])
  
  ## Loop to remove low values
  nsd1 = dim(samp1.dat)[1]
  keepID = NULL
  
  for (j in 1:nsd1) {
    if (samp1.dat$val[j] > tth) {
      keepID = c(keepID,j)
    }
  }
  # site.dat = site.dat[which(site.dat$val >= tth),]
  samp1.dat = samp1.dat[keepID,]
  nsd1 = dim(samp1.dat)[1]
  
  ## Sample 2
  samp2.dat = dat %>% 
    filter(sample == samps[(i+1)])
  
  ## Loop to remove low values
  nsd2 = dim(samp2.dat)[1]
  keepID = NULL
  
  for (j in 1:nsd2) {
    if (samp2.dat$val[j] > tth) {
      keepID = c(keepID,j)
    }
  }
  # site.dat = site.dat[which(site.dat$val >= tth),]
  samp2.dat = samp2.dat[keepID,]
  nsd2 = dim(samp2.dat)[1]
  
  ## Check sample sizes
  if (nsd1 >= 5 && nsd2 >= 5) {
    ## Sample one data
    trait1.dat = data.frame(lsla = samp1.dat$lsla.mn,
                            lhgt = samp1.dat$lhgt.mn,
                            lsdm = samp1.dat$lsdm.mn)
    trait1.dat = unique(trait1.dat)
    
    samp1.hv = hypervolume_box(trait1.dat,
                               samples.per.point = 1000, verbose = FALSE)
    
    
    ## Sample two data
    trait2.dat = data.frame(lsla = samp2.dat$lsla.mn,
                            lhgt = samp2.dat$lhgt.mn,
                            lsdm = samp2.dat$lsdm.mn)
    trait2.dat = unique(trait2.dat)
    
    samp2.hv = hypervolume_box(trait2.dat,
                               samples.per.point = 1000, verbose = FALSE)
    
    ## HV set
    samp.hvset = hypervolume_set(samp1.hv, samp2.hv, check.memory = FALSE)
    ovlp.hv[i] = hypervolume_overlap_statistics(samp.hvset)[2]
    
    # for (j in 1:nbit) {
    ovlp.rnd[i,] <- foreach(j=1:nbit, .combine = 'c') %do% {
      # oper <- foreach (j=1:nbit, .combine='comb', .multicombine=TRUE, .packages='hypervolume') %dopar% {
      ## Get random traits
      trait1.dat = data.frame(lsla = randomTraits(means = samp1.dat$lsla.mn, 
                                                  sds = samp1.dat$lsla.sd),
                              lhgt = randomTraits(means = samp1.dat$lhgt.mn, 
                                                  sds = samp1.dat$lhgt.sd),
                              lsdm = randomTraits(means = samp1.dat$lsdm.mn, 
                                                  sds = samp1.dat$lsdm.sd))
      trait1.dat = unique(trait1.dat)
      
      samp1.hv = hypervolume_box(trait1.dat,
                                 samples.per.point = 1000, verbose = FALSE)
      
      trait2.dat = data.frame(lsla = randomTraits(means = samp2.dat$lsla.mn, 
                                                  sds = samp2.dat$lsla.sd),
                              lhgt = randomTraits(means = samp2.dat$lhgt.mn, 
                                                  sds = samp2.dat$lhgt.sd),
                              lsdm = randomTraits(means = samp2.dat$lsdm.mn, 
                                                  sds = samp2.dat$lsdm.sd))
      trait2.dat = unique(trait2.dat)
      
      samp2.hv = hypervolume_box(trait2.dat,
                                 samples.per.point = 1000, verbose = FALSE)
      
      samp.hvset = hypervolume_set(samp1.hv, samp2.hv, 
                                   verbose = FALSE, check.memory = FALSE)
      
      hypervolume_overlap_statistics(samp.hvset)[2]
    } # Resample loop 
    
  } # Sample check
  
} # Site loop

save(samps, ovlp.hv, ovlp.rnd, 
     #nbit=1000
     # file="nampdTraits_alpha_all.RData")
     file="carana_beta.RData")

end.time = Sys.time()

# parallel::stopCluster(cl)

