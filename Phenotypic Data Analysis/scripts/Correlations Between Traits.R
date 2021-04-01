# Project: IWG_NAM_Shattering_Threshability
# Analysis - Correlations Between Traits
# Author: Kayla R. Altendorf
# Date: 11/27/2020

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
library("corrplot")
library("ggplot2")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"

# folder we're on
folder <- "/Phenotypic Data Analysis"

# declare where you want the output to go
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")


#### Step 1: Read in Emmeans ####
# get the id_frame from data.csv - this is needed to convert the emmeans identity (which was required to be balanced for proper DF calcs) to their NAM identity (not balanced)
dat <- read.csv(paste(out_path, "/data.csv", sep = ""), header = T)
id_frame <- dat %>% 
  select(famID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  distinct() %>%
  filter(! is.na(longID))

# which traits are we dealing with? 
dirs <- list.dirs(path = out_path, full.names = TRUE, recursive = FALSE)

# set locations in order
env <- c("stp17", "stp18", "tli17", "tli18")

# create empty dataframe
emmeans_env <- data.frame(matrix(NA, nrow = 1295, ncol = 14))
emmeans_all <- list()

dirs <- dirs[-14]

for (j in 1:length(env)) {
  for (i in 1:length(dirs)) {
    files <- list.files(dirs[i], full.names = TRUE)
    files <- files[grepl(env[j], files)]
    file <- files[grepl("genet", files)]
    
    if (length(file) > 1) { # if there's more than one file that matches the "env" and genet characteristics, go with the transformed on for the sake of the correlation analysis
      file <- files[grepl("transformed", files)]
    }

    emmeans <- read.table(file, head = T)
    emmeans1 <- emmeans %>% mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% dplyr::select(8, 3)
    emmeans2 <- left_join(id_frame, emmeans1, by = "famID_plantID3")
    emmeans_env[,i] <- emmeans2[,5] 
    colnames(emmeans_env)[i] <- str_split(dirs[i], "/")[[1]][11] 
  }
  emmeans_all[[j]] <- emmeans_env
}


#### Step 1: Make a figure of Correlations ####
# set function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# order and rename columns
# change i from 1-4 to create different environments. export each individually and edit in PowerPoint.
df <- emmeans_all[[i]] %>% dplyr::select(rachis_breaks_mean, floret_score_mean, threshability, floret_site_utilization, height, 
                                            thousand_grain_weight, area, florets_per_spikelet, spike_length, head_weight_per_spike, 
                                            yield_per_spike, emergence_percent, anthesis_score)

colnames(df) <- c("Brittle Rachis", "Floret Shattering", "Threshability", "Floret Site Utilization", "Height", "Thousand Grain Weight", 
                    "Seed Area", "Florets Spikelet-1", "Spike Length", "Headweight", "Yield Spike-1", 
                    "Spike Emergence", "Anthesis")

M <- cor(df, y = df, use = "complete.obs", method = c("pearson"))
p.mat <- cor.mtest(df)
corrplot(M, method = "color", type = "upper", order = "original", addCoef.col = "black", tl.col= "black", p.mat = p.mat, sig.level = 0.01, insig = "blank", diag = FALSE)




# altnerative option 

#### Step 2: Run Correlation Analysis ####
# set the order of how we want the correlation analysis to feature them 
trait_order <- c("rachis_breaks_mean", "floret_score_mean", "reproductive_tiller_ct", "threshability", "floret_site_utilization", "height", "thousand_grain_weight", 
                 "area", "florets_per_spikelet", "spike_length", "head_weight_per_spike", "yield_per_spike", "emergence_percent", "anthesis_score")
# and their names as they should appear in the table
pub_names <- c("Brittle Rachis (ct)", "Floret Shattering (scale)", "Reproductive Tiller (ct)", "Threshability (%)", "Floret Site Utilization (%)", 
               "Height (cm)", "Thousand Grain Weight (g)", "Seed Area (mm2)", "Florets Spikelet-1 (ct)", "Spike Length (cm)", "Headweight (g)", "Yield Spike-1 (g)", 
               "Spike Emergence (%)", "Anthesis (1-10)")

cor_out <- list()

for (i in 1:length(emmeans_all)){
  # first select the order we want
  emmeans_all[[i]] <- emmeans_all[[i]] %>% select(trait_order)
  colnames(emmeans_all[[i]]) <- pub_names
  cor <- rcorr(as.matrix(emmeans_all[[i]]), type = "pearson")
  cor.p <- as.data.frame(cor$P)
  cor.r <- as.data.frame(cor$r)
  cor.r <- cor.r %>% mutate_all(funs(round(., 2))) # round r to two digits
  cor.p.stars <- cor.p
  for (k in 1:ncol(cor.p)) {
    for (j in 1:nrow(cor.p)) {
      cor.p.stars[k,j] <- stars.pval(cor.p[k,j])
    }
  }
  cor.r.p <- cor.r
  for (k in 1:ncol(cor.p)) {
    for (j in 1:nrow(cor.r)) {
      cor.r.p[k,j] <- paste(cor.r[k,j], cor.p.stars[k,j], sep = "")
    }
  }
  cor.r.p[upper.tri(cor.r.p)] <- NA
  row.names(cor.r.p) <- pub_names
  cor_out[[i]] <- cor.r.p
}

# create output directory
dir.create(paste(out_path, "/trait_correlations", sep = ""))

# write out output
for (j in 1:length(env)) {
  write.table(cor_out[[j]], paste(out_path, "/trait_correlations/cor_sig_", env[j], ".txt", sep = ""), quote = F, row.names = T, col.names = NA, sep = "\t")
}



