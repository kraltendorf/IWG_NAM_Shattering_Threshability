# Project: IWG_NAM_Shattering_Threshability
# Analysis - Heritability Calculations
# Author: Kayla R. Altendorf
# Date: 11/28/2020

# load required packages
library("dplyr")
library("tidyr")
library("stringr")
library("tibble")
library("lme4")
library("lmerTest")
library("emmeans")
library("gtools")
library("Hmisc")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"

# folder we're on
folder <- "/Phenotypic Data Analysis"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")

# read in the data 
dat <- read.csv(paste(dir, folder, "/output", "/data.csv", sep = ""), header = T)

# make important terms as factor
dat$year <- as.factor(dat$year)
dat$rep <- as.factor(dat$rep)


#### Step 1: Calculate Broad and Narrow Sense Heritability ####
# list traits
all_traits <- c("rachis_breaks_mean", "floret_score_mean", "threshability")

# set vectors for iterating through years and locations
loc <- c("STP", "TLI")

# create a data frame for the output
herit_df <- data.frame(env = c("STP", "TLI"), broad = NA, narrow = NA, trait = NA)
herit_df_list <- replicate(3, herit_df, simplify = FALSE) # how many traits are being calculated? 

for (i in 1:length(all_traits)) {
  for (j in 1:length(loc)) {
    
    # extract data 
    dat_loc_year <- dat[dat$loc == loc[j],]
    
    # calculate broad sense on a genet mean basis
    # if the trait is reproductive tiller number, it requires a transformation, if not, proceed with regular equation
    if (all_traits[i] == "reproductive_tiller_ct") {
      formula <- paste0("sqrt(reproductive_tiller_ct)", " ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)" , sep = "")
      model <- lmer(formula, data = dat_loc_year)
      }
    if (all_traits[i] == "threshability") {
      tran <- make.tran(type = "asin.sqrt", 100)
      model <- with(tran, lmer(linkfun(threshability) ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep), data = dat_loc_year))
      }
    if (! all_traits[i] %in% c("reproductive_tiller_ct", "threshability")) {
      formula <- paste0(all_traits[i], " ~ (1|famID:plantID3) + (1|year) + (1|famID:plantID3:year) + (1|rep) + (1|famID:plantID3:rep)", sep = "")
      model <- lmer(formula, data = dat_loc_year)
    }
    
    output <- as.data.frame(VarCorr(model))
    vg <- output$vcov[3] # extracting appropriate variance components
    vgy <- output$vcov[1]
    vgr <- output$vcov[6]
    heritability_broad <- vg / ((vg) + (vgy / 2) + (vgr / 4))
    
    # clean out variables before next iteration
    formula <- NA
    model <- NA
    output <- NA
    vg <- NA
    vgy <- NA
    vgr <- NA
      
    # calculate narrow sense according to falconer, again if reproductive tiller, requires trans
    # in this case, family and rep are random
    if (all_traits[i] == "reproductive_tiller_ct") {
      formula <- paste0("sqrt(reproductive_tiller_ct)", " ~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)" , sep = "")
      model <- lmer(formula, data = dat_loc_year)
      }
    if (all_traits[i] == "threshability") {
      tran <- make.tran(type = "asin.sqrt", 100)
      model <- with(tran, lmer(linkfun(threshability) ~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep), data = dat_loc_year))
      }
    if (! all_traits[i] %in% c("reproductive_tiller_ct", "threshability")) {
      formula <- paste0(all_traits[i], "~ (1|famID) + (1|year) + (1|famID:year) + (1|rep) + (1|famID:rep)" , sep = "")
      model <- lmer(formula, data = dat_loc_year)
    }
    
    output <- as.data.frame(VarCorr(model))
    vf <- output$vcov[3] # extracting appropriate variance components
    ve <- output$vcov[6]
    heritability_narrow <- (vf * 4) / ((vf * 4) + (ve))
      
    # clean out variables before next iteration
    model <- NA
    output <- NA
    vf <- NA
    ve <- NA

    # output result into heritability dataframe
    herit_df_list[[i]][j,2] <- heritability_broad
    herit_df_list[[i]][j,3] <- heritability_narrow
    herit_df_list[[i]][j,4] <- all_traits[i]
    
    # empty out these variables
    heritability_broad <- NA
    heritability_narrow <- NA
  }
}

herit_df <- do.call("rbind", herit_df_list)
herit_df_wide <- herit_df %>% 
  pivot_wider(values_from = c("broad", "narrow"), names_from = c("env")) %>%
  select(trait, broad_STP, narrow_STP, broad_TLI, narrow_TLI)

# calculate and append averages
herit_df_wide_round <- herit_df_wide %>% 
  mutate_if(is.numeric, funs(round(., 2)))

# write out result
write.table(herit_df_wide_round, paste(out_path, "/heritabilities.txt", sep = ""),  quote = F, row.names = F, sep = "\t")


#### Step 2: Edit TLI calculations ####
# these traits, shattering and threshability, were really affected by the drought in 2018 at TLI
# thus, this really brought down the estimates. So here we'll edit the calculations to exclude that year.

for (i in 1:length(all_traits)) {
    
    # extract data 
    dat_loc <- dat[dat$loc == "TLI",]
    dat_loc_year <- dat_loc[dat_loc$year == "2017",]
    
    # calculate broad sense on a genet mean basis
    if (all_traits[i] == "threshability") {
      tran <- make.tran(type = "asin.sqrt", 100)
      model <- with(tran, lmer(linkfun(threshability) ~ (1|famID:plantID3) + (1|rep), data = dat_loc_year))
    }
    if (! all_traits[i] %in% c("threshability")) {
      formula <- paste0(all_traits[i], " ~ (1|famID:plantID3) + (1|rep)", sep = "")
      model <- lmer(formula, data = dat_loc_year)
    }
    
    output <- as.data.frame(VarCorr(model))
    vg <- output$vcov[1] # extracting appropriate variance components
    vgr <- output$vcov[3]
    heritability_broad <- vg / ((vg) + (vgr / 2))
    
    herit_df_wide_round[i, 4] <- round(heritability_broad, 2)
    
    # clean out variables before next iteration
    formula <- NA
    model <- NA
    output <- NA
    vg <- NA
    vgy <- NA
    vgr <- NA
    
    # calculate narrow sense according to falconer, again if reproductive tiller, requires trans
    # in this case, family and rep are random
    
    if (all_traits[i] == "threshability") {
      tran <- make.tran(type = "asin.sqrt", 100)
      model <- with(tran, lmer(linkfun(threshability) ~ (1|famID) + (1|rep) + (1|famID:rep), data = dat_loc_year))
    }
    if (! all_traits[i] %in% c("threshability")) {
      formula <- paste0(all_traits[i], "~ (1|famID) + (1|rep) + (1|famID:rep)" , sep = "")
      model <- lmer(formula, data = dat_loc_year)
    }
    
    output <- as.data.frame(VarCorr(model))
    vf <- output$vcov[2] # extracting appropriate variance components
    ve <- output$vcov[4]
    heritability_narrow <- (vf * 4) / ((vf * 4) + (ve))
    
    # clean out variables before next iteration
    model <- NA
    output <- NA
    vf <- NA
    ve <- NA
    
    herit_df_wide_round[i, 5] <- round(heritability_narrow, 2)
}

# write out the results where I only include 2017 in the TLI estimates
write.table(herit_df_wide_round, paste(out_path, "/heritabilities_tli2017.txt", sep = ""),  quote = F, row.names = F, sep = "\t")

