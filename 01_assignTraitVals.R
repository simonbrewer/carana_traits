###############################################################################
## Assign BIEN traits to pollen taxa

## Adding scaling here, prior to attribution to pollen

require(ggplot2)
require(GGally)
library(dplyr)

## Get the trait values and bin them for each species
dat = read.csv("./NativeCarana_6traitsBIEN.csv")

dat = dat %>%
  select(sppID = ID, scrubbed_species_binomial, leafarea, plantheight, seedmass)

# dat$leafarea = as.numeric(dat$leafarea)
# dat$plantheight[which(dat$plantheight == ".")] <- "-9999"
# dat$plantheight = as.numeric(dat$plantheight)
# dat$seedmass = as.numeric(dat$seedmass)
# 
dat$leafarea[dat$leafarea == -9999] <- NA
dat$plantheight[dat$plantheight == -9999] <- NA
dat$seedmass[dat$seedmass == -9999] <- NA

dat$leafarea[dat$leafarea == 0] <- NA
dat$plantheight[dat$plantheight <= 0] <- NA
dat$seedmass[dat$seedmass == 0] <- NA

dat$lsla = scale(log10(dat$leafarea))
dat$lhgt = scale(log10(dat$plantheight))
dat$lseedm = scale(log10(dat$seedmass))

lsla.med = tapply(dat$lsla, dat$sppID, median, na.rm=TRUE)
lhgt.med = tapply(dat$lhgt, dat$sppID, median, na.rm=TRUE)
lseedm.med = tapply(dat$lseedm, dat$sppID, median, na.rm=TRUE)

lsla.avg = tapply(dat$lsla, dat$sppID, mean, na.rm=TRUE)
lhgt.avg = tapply(dat$lhgt, dat$sppID, mean, na.rm=TRUE)
lseedm.avg = tapply(dat$lseedm, dat$sppID, mean, na.rm=TRUE)

lsla.sd = tapply(dat$lsla, dat$sppID, sd, na.rm=TRUE)
lhgt.sd = tapply(dat$lhgt, dat$sppID, sd, na.rm=TRUE)
lseedm.sd = tapply(dat$lseedm, dat$sppID, sd, na.rm=TRUE)

## Read the species/taxa matrix
dat2 = read.table("Carana_sppXpollentaxa_v4.csv", sep=',')

plev = dat2[1,-c(1,2)]
pvar = dat2[2,-c(1,2)]
pvar = as.numeric(as.matrix(pvar))
paccvar = dat2[3,-c(1,2)]
pvarcode = dat2[4,-c(1,2)]
pvarcode = as.character(as.matrix(pvarcode))
pvarname = dat2[5,-c(1,2)]
pvarname = as.character(as.matrix(pvarname))
pgroup = dat2[6,-c(1,2)]
pgroup = as.character(as.matrix(pgroup))

tryID = dat2[-seq(1:6),1]
tryname = dat2[-seq(1:6),2]

dat2 = dat2[-seq(1:6),-c(1:2)]
ttable = matrix(as.numeric(as.matrix(dat2)), ncol=length(pvar))

rownames(ttable) = tryID
colnames(ttable) = as.numeric(as.matrix(pvar))

npvar = length(pvar)
poll.trait = data.frame(lsla.med = rep(NA, npvar),
                        lsla.mn = rep(NA, npvar),
                        lsla.sd = rep(NA, npvar),
                        lsla.n = rep(NA, npvar),
                        lhgt.med = rep(NA, npvar),
                        lhgt.mn = rep(NA, npvar),
                        lhgt.sd = rep(NA, npvar),
                        lhgt.n = rep(NA, npvar),
                        lsdm.med = rep(NA, npvar),
                        lsdm.mn = rep(NA, npvar),
                        lsdm.sd = rep(NA, npvar),
                        lsdm.n = rep(NA, npvar))

## Now loop through pollen taxa and assign trait values
cat("Traits to PVar assignment",file="taxaAssign.log")
for (i in 1:npvar) {
# for (i in 1:10) {
  print(paste(i,pvar[i],pvarname[i]))
  tvec = ttable[,i]
  
  tmpTryID = tryID[which(tvec==1)]
  ntry = length(tmpTryID)
  cat(as.character(i),file="taxaAssign.log", append=TRUE)
  cat("\n",file="taxaAssign.log", append=TRUE)
  cat(paste(pvar,pvarcode,pvarname),file="taxaAssign.log", append=TRUE)
  cat("\n",file="taxaAssign.log", append=TRUE)
  cat(as.character(tryname[which(tvec==1)]),file="taxaAssign.log", append=TRUE)
  cat("\n",file="taxaAssign.log", append=TRUE)
  
  if (ntry > 0) {
    cat(as.character(tryname[which(tvec==1)]),file="taxaAssign.log", append=TRUE)
    cat("\n",file="taxaAssign.log", append=TRUE)
    # tmpsla = rep(NA, ntry)
    # tmphgt = rep(NA, ntry)
    # tmpseedm = rep(NA, ntry)
    tmpsla = NULL
    tmphgt = NULL
    tmpsdm = NULL
    
    for (j in 1:ntry) {
      #print(j)
      dat.sub = dat[which(dat$sppID == tmpTryID[j]), ]
      # tmpsla[j] = median(dat$lsla[which(dat$sppID == tmpTryID[j])], na.rm=TRUE)
      # tmphgt[j] = median(dat$lhgt[which(dat$sppID == tmpTryID[j])], na.rm=TRUE)
      # tmpseedm[j] = median(dat$lseedm[which(dat$sppID == tmpTryID[j])], na.rm=TRUE)
      tmpsla = c(tmpsla, dat.sub$lsla)
      tmphgt = c(tmphgt, dat.sub$lhgt)
      tmpsdm = c(tmpsdm, dat.sub$lseedm)
    }
    poll.trait$lsla.med[i] = median(tmpsla, na.rm=TRUE)
    poll.trait$lsla.mn[i] = mean(tmpsla, na.rm=TRUE)
    poll.trait$lsla.sd[i] = sd(tmpsla, na.rm=TRUE)
    poll.trait$lsla.n[i] = length(which(!is.na(tmpsla)))
    poll.trait$lhgt.med[i] = median(tmphgt, na.rm=TRUE)
    poll.trait$lhgt.mn[i] = mean(tmphgt, na.rm=TRUE)
    poll.trait$lhgt.sd[i] = sd(tmphgt, na.rm=TRUE)
    poll.trait$lhgt.n[i] = length(which(!is.na(tmphgt)))
    poll.trait$lsdm.med[i] = median(tmpsdm, na.rm=TRUE)
    poll.trait$lsdm.mn[i] = mean(tmpsdm, na.rm=TRUE)
    poll.trait$lsdm.sd[i] = sd(tmpsdm, na.rm=TRUE)
    poll.trait$lsdm.n[i] = length(which(!is.na(tmpsdm)))

  } else {
    cat("No matches",file="taxaAssign.log", append=TRUE)
    cat("\n",file="taxaAssign.log", append=TRUE)
    
  }
  
}

out.df = data.frame(pvar = pvar, pvarcode = pvarcode, pgroup = pgroup, 
                    pvarname = pvarname, poll.trait)
write.csv(out.df, "pvarTraitVals.csv", row.names = FALSE)
