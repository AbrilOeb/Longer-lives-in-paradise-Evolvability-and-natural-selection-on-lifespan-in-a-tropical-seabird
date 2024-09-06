# R SCRIPT FOR RUNNING THE -ANIMAL MODEL-

# Model will be divided into shorter chains to speed up the process
# THIS SCRIPT ALLOWS ME TO SPLIT MY MODEL IN SHORTER CHAINS BY SUMBIMITTING MULTIPLE SHORTER JOBS


# Making the number process the argument, to use it when saving and save each job with a different name
args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("must give a number at the command line which will be used to create the filename", call.=FALSE)
} else {
  print("this script would write file number")
  print(args)
  # use argument to create filename
}

# The remote computer will run the model and give it back to me in a RDS format
# Since I am running multiple jobs with the same script I will get multiple RDS as my output

# MODEL: HERITABILITY OF LRS
#        Animal Model using a HURDLE distribution

#------------------------------------------------------------------------------
### PACKAGES ###
# Create a personal library to install packages into. 
# This has been created and ALL REQUIRED PROGRAMS have been installed using PATUNG
# dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)

# Add to the Path so remote computers can know where to look for my packages
.libPaths(Sys.getenv("R_LIBS_USER")) 

#------------------------------------------------------------------------------
### DATA FILES FOR THE MODEL ###

#Re-direct to the working directory to easily load files

setwd("/srv/home/abriloeb/wp1/h2/lrs/all")

#BFB_PEDIGREE

bfb_ped <- readRDS("bfb_ped.rds")

#BFB_DATABASE

bfb_all <- readRDS("bfb_all.rds")

#------------------------------------------------------------------------------
### CHOOSING AND USING MY DATA

# I can't give the whole database to the model because it slows down theprocess. 
# Therefore, I need to carefully check what information I am giving as data.

# BFB_PEDIGREE <- Will always use the complete and checked pedigree

# BFB_DATABASE <- Will always change depending on the variables I am using.

# THIS MODEL: 
# animal, sire and dam columns must ALWAYS be included
# variable of interest: lrs
# fixed effects: none
# random effects: nest, cohort

# Selecting only the relevant columns for the dataframe

library(tidyverse)

bfb_all <- select(bfb_all, animal, sire, dam, nest, cohort, lrs.fledged)

#------------------------------------------------------------------------------
### RUNNING THE MODEL ###

library(MCMCglmm)

#1.-lrs_hurdle

# Setting up the prior

prior_lrs_hurdle<- 
  list(G=list(G1=list(V=diag(2),nu=2),
              G2=list(V=diag(2),nu=2),
              G3=list(V=diag(2),nu=2)),
       R=list(V=diag(2),fix=2,nu=2))

# Running the model:

lrs_hurdle_batch<- MCMCglmm(lrs.fledged ~ trait-1 , 
                       random = ~ us(trait):animal + 
                         us(trait):nest +
                         us(trait):cohort,
                       rcov = ~ idh(trait):units ,
                       family = c("hupoisson"),
                       pedigree = bfb_ped ,
                       data = bfb_all ,
                       nitt = 55000 ,
                       thin = 2500 , 
                       burnin = 5000 ,
                       verbose = TRUE, 
                       prior = prior_lrs_hurdle)

#----Saving the chains in RDS format

#Within the remote server:

saveRDS(lrs_hurdle_batch,
        file = paste("lrs_hurdle_batch_", args, ".rds", sep=""))


#END...


