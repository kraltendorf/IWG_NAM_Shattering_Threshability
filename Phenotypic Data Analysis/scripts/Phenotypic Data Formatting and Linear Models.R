# Project: IWG_NAM_Shattering_Threshability
# Analysis - Phenotypic Data Formatting and Linear Models
# Author: Kayla R. Altendorf
# Date: 11/24/2020

# load required packages
library("dplyr")
library("tidyr")
library("stringr")
library("tibble")
library("lme4")
library("lmerTest")
library("multcompView")
library("emmeans")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"

# folder we're on
folder <- "/Phenotypic Data Analysis"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")

# read in NAM data from the IWG database
dat <- read.table(paste(dir, "Phenotypic Data Analysis/data/NAM_Data.txt", sep = ""), header = T, sep = "\t")


#### Step 1: Format the Phenotypic Data #### 
# since shattering was measured on three spikes per plant, we need to calculate the mean of these subsamples
# and round it to the nearest integer

# NOTE: it's important to know that the floret shattering scale originally included a test where, if no
# florets broke off in the drop test method, the top floret in a random spikelet was bent away from
# the spike. If it could bend to a 45 degree angle without breaking off, it was given a 0, if it broke off
# it was given a 0.5. After the data was collected, I decided to ignore this part of the test because I found it to be somewhat 
# subjective and not informative. Thus here we will replace the 0.5 values with 0, 
# as well as replace two data entry errors -- where this trait had an incorrect value of 0.45 and 0.1. 

for (i in 1:nrow(dat)) {
  if (dat$trait_id[i] == "FLSHAT" & dat$phenotype_value[i] %in% c(0.5, 0.45, 0.1)) {
    dat$phenotype_value[i] <- 0 }
}

# create a dataframe for merging purposes
trait_names <- data.frame(trait_id = c("BRSHAT", "FLSHAT", "NKD"), 
                          trait_id_full = c("rachis_breaks_mean", "floret_score_mean", "threshability"))

dat1 <- dat %>% filter(trait_id %in% trait_names$trait_id) %>%
  left_join(trait_names, by = "trait_id") %>% # keep only the traits that are relevant for this analysis
  dplyr::rename(year = phenotype_year, # rename cols to match personal preference
                famID = family_name, 
                col = range) %>%
  mutate(loc = str_replace(substr(experiment_id, 4, 6), "SAL", "TLI"), # extract location, replacing SAL with TLI
         parent = substr(germplasm_id, 6, 6), # extract parent (C for common, D for donor, P for parent)
         merge_col = paste(germplasm_id, loc, year, rep, sep = "_"),
         trait_id_full = case_when(trait_id_full %in% c("floret_score_mean", "rachis_breaks_mean") ~ paste(trait_id_full, sample_number, sep = ""), # give each subsample a unique trait_id_full
                                   ! trait_id_full %in% c("floret_score_mean", "rachis_breaks_mean") ~ paste(trait_id_full))) %>%
 filter(parent != "P") %>% # exclude parents - hash this out to export that data with the parents (Line 75)
  select(famID, germplasm_id, loc, year, rep, trait_id_full, phenotype_value, plant_id, merge_col) %>% 
  pivot_wider(names_from = trait_id_full, values_from = phenotype_value)

# extract column names
colnames(dat1)
dat_cols <- colnames(dat1)[8:14]

# calculate the means and eliminate the original cols
dat2 <- dat1 %>% mutate_at(vars(dat_cols), funs(as.numeric(as.character(.)))) %>%
  mutate(floret_score_mean =  round(rowMeans(dplyr::select(., dat_cols[1:3]), na.rm = T)), 0, 
         rachis_breaks_mean = round(rowMeans(dplyr::select(., dat_cols[4:6]), na.rm = T)), 0) %>%
  dplyr::select(-dat_cols[1:6], -famID, -germplasm_id, -loc, -year, -rep, -plant_id, -'0')


# we'll need the parents later, so hash (##) out line 59 to eliminate that filtering step
# and export here
write.csv(dat2, paste(out_path, "data_including_parents_selfs_not_na.csv", sep = "/"), row.names = F)

# the database includes *all* entries, including parents, plants that were identified later as selfs, and so on.
# futhermore, the population itself is not balanced (e.g. unequal numbers of individuals within families, 
# entries per location and year), which causes problems with the ANOVA
# to address this, we will load in a backbone dataset, which is balanced using NA values,
# and we'll format the data to match. 

# read in the backbone csv
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/data/backbone.csv", sep = ""), header = T)

# left join with the data
dat3 <- left_join(backbone, dat2, by = "merge_col") %>% select(-merge_col)
dat3 %>% group_by(loc, year, rep, famID) %>% tally() # make sure it's still balanced

# before we change selfs to NA, export this
write.csv(dat3, paste(out_path, "data_selfs_not_na.csv", sep = "/"), row.names = F)

# change all selfs to NA
for (i in 1:nrow(dat3)) {
  if (! is.na(dat3$self[i])) { # if self is not NA (e.g. if it is outcross or self)
    dat3[i, 13:15] <- NA # change all phenotype data columns to NA
  }
}

colnames(dat3)

# write out the final dataset
write.csv(dat3, paste(out_path, "data.csv", sep = "/"), row.names = F)

# read it back in to convert everything to numerics
dat <- read.csv(paste(out_path, "data.csv", sep = "/")) %>%
  mutate_at(vars(plantID3, year, rep, famID), funs(as.factor(.))) %>%
  filter(parent != "P")


#### Step 2: Analysis of Variance Combined Across Environments ####

# threshability - where it requires a asin.sqrt using emmeans 'make.tran' and 'linkfun'
tran <- make.tran(type = "asin.sqrt", 100)
model <- with(tran, lmer(linkfun(threshability) ~ famID * loc * year + (1|loc:rep) + (1|loc:famID:plantID3), data = dat))
anova <- anova(model)
plot(model)
dir.create(paste(out_path, "threshability", sep = "/"))
emmeans_loc_year <- CLD(emmeans(model, ~  year * loc, type = "response"), Letters = c(LETTERS))
write.table(anova, paste(out_path, "/threshability/threshability_anova.txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
write.table(emmeans_loc_year, paste(out_path, "/threshability/threshability_emmeans_loc_year.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

#### Step 3: Analysis of Variance within Environments #### 
# filter data by location
stp17 <- filter(dat, loc == "STP" & year == "2017")
stp18 <- filter(dat, loc == "STP" & year == "2018")
tli17 <- filter(dat, loc == "TLI" & year == "2017")
tli18 <- filter(dat, loc == "TLI" & year == "2018")

# set some vectors for iteration
loc_list <- list(stp17, stp18, tli17, tli18)
loc_names <- c("stp17", "stp18", "tli17", "tli18")

# generalized linear models - for count data
# rachis_breaks_mean
dir.create(paste(out_path, "rachis_breaks_mean", sep = "/"))
for (j in 1:length(loc_list)) {
  model <- glm(rachis_breaks_mean ~ famID/plantID3 + rep, data = loc_list[[j]])
  anova <- anova(model)
  emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID), Letters = c(LETTERS)))
  emmeans_genet <- as.data.frame(emmeans(model,  ~ plantID3|famID))
  write.table(anova, paste(out_path, "/rachis_breaks_mean/rachis_breaks_mean_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
  write.table(emmeans_fam, paste(out_path, "/rachis_breaks_mean/rachis_breaks_mean_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(emmeans_genet, paste(out_path, "/rachis_breaks_mean/rachis_breaks_mean_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
}

# floret_score_mean
dir.create(paste(out_path, "floret_score_mean", sep = "/"))
for (j in 1:length(loc_list)) {
  model <- glm(floret_score_mean ~ famID/plantID3 + rep, data = loc_list[[j]])
  anova <- anova(model)
  emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID), Letters = c(LETTERS)))
  emmeans_genet <- as.data.frame(emmeans(model,  ~ plantID3|famID))
  write.table(anova, paste(out_path, "/floret_score_mean/floret_score_mean_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
  write.table(emmeans_fam, paste(out_path, "/floret_score_mean/floret_score_mean_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(emmeans_genet, paste(out_path, "/floret_score_mean/floret_score_mean_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
}

# threshability 
for (j in 1:length(loc_list)) {
  for (i in 1:length(scale)) {
    tran <- make.tran(type = "asin.sqrt", 100)
    model <- with(tran, lm(linkfun(threshability) ~ famID/plantID3 + rep, dat = loc_list[[j]]))
    anova <- anova(model)
    emmeans_fam <- as.data.frame(CLD(emmeans(model, ~ famID, type = "response"), Letters = c(LETTERS)))
    emmeans_genet <- as.data.frame(emmeans(model, ~ plantID3|famID, type = "response"))
    write.table(anova, paste(out_path, "/threshability/threshability_anova_", loc_names[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
    write.table(emmeans_fam, paste(out_path, "/threshability/threshability_emmeans_fam_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    if (scale[i] == "transformed") {
      emmeans_genet <- as.data.frame(emmeans(model, ~ plantID3|famID))
      write.table(emmeans_genet, paste(out_path, "/threshability/threshability_emmeans_genet_emmeans_genet_transformed_scale_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
      }
    if (scale[i] == "response") {
      emmeans_genet <- emmeans(model, ~ plantID3|famID, type = "response")
      write.table(emmeans_genet, paste(out_path, "/threshability/threshability_emmeans_genet_", loc_names[j], ".txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
}
