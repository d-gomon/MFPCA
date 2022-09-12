rm(list=ls())
dir <- "C:/Users/danie/OneDrive - Universiteit Leiden/Supervision_Daniel/Research/04_prediction_long+highdim/ADNI data/"
load(paste0(dir, "adni_cleaned_20220520.RData"))

library(funData)
#library(MFPCA)
load_all()
library(ggplot2)
dat <- df.long_censored_prc
#remove 1 person with follow up 5 days after entry
dat <- dat[-c(which(dat$VISCODE == "m0")),]
#use df.long_censored_prc

Npatients <- length(unique(df.long_censored_prc$RID))
times <- dat$VISCODE
times <- substring(times, first = 2)
times[which(times == "l")] <- "0"
times <- as.integer(times)
dat$times <- times

# RESCALE VARIANCE 1 ALL VARIABLES!!!!! MEAN 0 TODOTODOTODOTODO


#Write function to transform BioMarker to funData format.

WholeBrain <- create_x("WholeBrain")

Hippocampus = create_x("Hippocampus")

Fusiform <- create_x("Fusiform")

mdat <- multiFunData(WholeBrain, Hippocampus, Fusiform)

mfpca_test <- MFPCA::MFPCA(mdat,M = 4,
                           uniExpansions = list(list(type = "uFPCA", pve = 0.93),
                                                list(type = "uFPCA", pve = 0.93),
                                                list(type = "uFPCA", pve = 0.93)),
                           fit = TRUE)

scores_manual <- predict_mscores(mfpca_test, mFData_pred = mdat)

all.equal(mfpca_test$scores, scores_manual$mscores_test)

scores_manual2 <- predict_mscores(mfpca_test, mFData_pred = mdat[1:400,])

all.equal(mfpca_test$scores, scores_manual2$mscores_train)




