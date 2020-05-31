## Test trait assignment

## Get the trait values and bin them for each species
dat = read.csv("./NativeCarana_6traitsBIEN.csv")

dat = dat %>%
  select(sppID = ID, scrubbed_species_binomial, leafarea, plantheight, seedmass)

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

## Try 2 and 4
myi = 4

print(paste(myi,pvar[myi],pvarname[myi]))
tvec = ttable[,myi]

tmpTryID = tryID[which(tvec==1)]
ntry = length(tmpTryID)
sd(c(5.85, 5.85, 5.85, 222, 106))

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
median(tmpsla, na.rm=TRUE)
mean(tmpsla, na.rm=TRUE)
sd(tmpsla, na.rm=TRUE)
length(which(!is.na(tmpsla)))
median(tmphgt, na.rm=TRUE)
mean(tmphgt, na.rm=TRUE)
sd(tmphgt, na.rm=TRUE)
length(which(!is.na(tmphgt)))
median(tmpsdm, na.rm=TRUE)
mean(tmpsdm, na.rm=TRUE)
sd(tmpsdm, na.rm=TRUE)
length(which(!is.na(tmpsdm)))
