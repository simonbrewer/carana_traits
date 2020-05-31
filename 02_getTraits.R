###############################################################################
## Assigns trait values down cord

library(dplyr)

## Trait values by pvar
ptraitvals = read.csv("pvarTraitVals.csv")
## Taxonomic groups assignments
## Output data frame
out.df = NULL

## Log file
if (file.exists("taxontrait.log")) {
  file.remove("taxontrait.log")
}

log_con <- file("taxontrait.log", open="a")

missing_taxa = NULL

## Get surface samples
dat = read.csv("Carana_percentages_v1.1.csv")
pvars = as.numeric(substr(colnames(dat)[-seq(1:7)], 2, 10))
nvars = length(pvars)
ptaxoncode = dat[2,-seq(1:7)]
grps = unlist(dat[4,-seq(1:7)])

ent = dat[-seq(1:4),4]
lon = dat[-seq(1:4),6]
lat = dat[-seq(1:4),5]
alt = dat[-seq(1:4),7]

poll = dat[-seq(1:4),-seq(1:7)]
poll <- mutate_all(poll, function(x) as.numeric(as.character(x)))

nsamp = nrow(poll)

## Make up temporary ID list
traitID = rep(NA, nvars)
# grps = rep(0, nvars)
for (j in 1:nvars) {
  tmpID = which(ptraitvals$pvar == pvars[j])
  # grpID = which(pgrp$Var == pvars[j])
  if (length(tmpID) != 1) { 
    # stop("HELP!!!")
    print(paste(j, "Missing taxon:", pvars[j]))
    cat(paste(j, "Missing taxon:", pvars[j]), file = log_con, sep="\n")
    missing_taxa = c(missing_taxa, pvars[j])
  } else {
    traitID[j] = tmpID
    # grps[j] = pgrp$Group[grpID]
    # print(paste(pgrp$Group[grpID], pgrp$VarName[grpID]))
  }
}

## Vectorize the pollen
nsamp = dim(poll)[1]
poll.df = data.frame(ent = rep(ent, each = nvars),
                     sample = rep(1, each = nsamp * nvars),
                     agebp = rep(0, each = nsamp * nvars),
                     lon = rep(lon, each = nvars),
                     lat = rep(lat, each = nvars),
                     alt = rep(alt, each = nvars),
                     pvar = rep(pvars, nsamp), 
                     val = c(t(as.matrix(poll))) / 100,
                     grp = rep(grps, nsamp),
                     ptraitvals[traitID,-c(1,2,3)], row.names = NULL)

# write.csv(poll.df, "test.csv", row.names = FALSE)
## Cut samples with 0 abundance
poll.df = poll.df %>%
  filter(val > 0 & !is.na(pvarname))

write.table(poll.df, "caranaPollTraitVals.csv", sep = ",", row.names = FALSE)
close(log_con)
