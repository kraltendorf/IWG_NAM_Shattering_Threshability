# Project: IWG_NAM_Shattering_Threshability
# Analysis - Variation Among Progeny, Maternal Effects
# Author: Kayla R. Altendorf
# Date: 11/24/2020

# location of github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"

# folder we're on
folder <- c("/Phenotypic Data Analysis")

# declare where you want the output to go 
out_path <- paste(dir, folder, "/output", sep = "")
in_path <- paste(dir, folder, "/data", sep = "")


#### Step 1: Read in the Data that Includes the Parents #### to create this see lines 59 and 72-73 of Phenotypic Data Formatting and Linear Models
dat <- read.csv(paste(dir, folder, "/output/data_including_parents_selfs_not_na.csv", sep = ""), header = T) %>%
  separate(merge_col, into = c("longID", "loc", "year", "rep"), sep = "_") %>%
  mutate(famID = substr(longID, 1, 5), 
         parent = substr(longID, 6, 6))


#### Step 2: Calculate Parental Means ####
# set vectors
# environments
loc <- c("STP", "STP", "TLI", "TLI")
year <- c("2017", "2018", "2017", "2018")
envs <- c("stp17", "stp18", "tli17", "tli18")
mother <- c("C", "D")

# families
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")

# list of traits we're dealing with 
trait <- c("rachis_breaks_mean", "floret_score_mean", "threshability")
# separate out common and donor parent means
common_parent_means <- list()
donor_parent_means <- list()

file <- paste(out_path, "/parental_means.txt", sep = "")

for (i in 1:length(trait)) {
  fam <- dat %>% filter(parent == "P") %>% select(famID, loc, year, trait[i])
  fam_sum <- fam %>% group_by(loc, year, famID) %>% summarise_all(funs(mean = mean), na.rm = T)
  
  # export fam_sum while we're here
  fam_sum_out <- fam_sum %>% pivot_wider(id_cols = famID, names_from = c("loc", "year"), values_from = "mean")
  fam_sum_out[,-1] <-round(fam_sum_out[,-1],2)
  
  cat(paste("\t", trait[i], sep = ""), file = file, append = T)
  cat("\n", file = file, append = T)
  write.table(fam_sum_out, file, row.names = T, col.names = T, sep = "\t", quote = F, append = T)
  cat("\n", file = file, append = T)
  ### continue
  
  common_parent_means[[i]] <- fam_sum %>% filter(famID == "WGN59")
  donor_parent_means[[i]] <- fam_sum %>% filter(famID != "WGN59")
  
  for (j in 1:length(donor_parent_means)) {
    donor_means <- donor_parent_means[[j]] %>% mutate(dif = NA)
    for (k in 1:nrow(donor_means)) {
      if (donor_means$loc[k] == "STP" & donor_means$year[k] == "2017") {
        donor_means$dif[k] <- abs(donor_means$mean[k] - common_parent_means[[j]][1, 4])[,1] }
      else if (donor_means$loc[k] == "STP" & donor_means$year[k] == "2018") {
        donor_means$dif[k] <- abs(donor_means$mean[k] - common_parent_means[[j]][2, 4])[,1] }
      else if (donor_means$loc[k] == "TLI" & donor_means$year[k] == "2017") {
        donor_means$dif[k] <- abs(donor_means$mean[k] - common_parent_means[[j]][3, 4])[,1] }
      else if (donor_means$loc[k] == "TLI" & donor_means$year[k] == "2018") {
        donor_means$dif[k] <- abs(donor_means$mean[k] - common_parent_means[[j]][4, 4])[,1] }
    }
    donor_parent_means[[j]] <- donor_means
  }
}



#### Step 3: Calculate Progeny Means ####
# uninstall plyr - or else the group_by summarise command fails
progeny_means <- list()

for (i in 1:length(trait)) {
  files <- list.files(paste(out_path, "/", trait[i], sep = ""), pattern = "emmeans_genet", full.names = T)
  files <- files[!grepl("transformed", files)]
  tables <- lapply(files, read.table, header = T, fill = TRUE)
  for (j in 1:length(tables)) {
    tables[[j]]$loc <- loc[j]
    tables[[j]]$year <- year[j]
  }
  tables_all <- do.call("rbind", tables)
  colnames(tables_all)[3] <- c("emmean")
  tables_end <- tables_all %>% 
    dplyr::select(famID, loc, year, emmean) %>%
    group_by(loc, year, famID) %>%
    summarise(trait_mean = mean(emmean, na.rm = T), trait_sd = sd(emmean), trait_min = min(emmean), trait_max = max(emmean))
  progeny_means[[i]] <- tables_end 
}

#### Step 4: Do Variable Parents Result in Variable Progeny? ####
# make a dataframe for the output
cor.result.df <- data.frame()
cor.estimate <- data.frame(matrix(NA, nrow = 3, ncol = 4))
cor.pval <- data.frame(matrix(NA, nrow = 3, ncol = 4))
colnames(cor.estimate) <- envs
colnames(cor.pval) <- envs

# iterate through traits and environments
for (i in 1:length(trait)) {
  means_bound <- bind_cols(donor_parent_means[[i]], progeny_means[[i]][,4:5]) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    droplevels() %>%
    mutate(loc = as.factor(loc),
           year = as.factor(year))
  for (j in 1:length(loc)) {
    means_loc <- means_bound[means_bound$loc == loc[j],]
    means_loc_year <- means_loc[means_loc$year == year[j],]
    result.df <- cor.test(means_loc_year$dif, means_loc_year$trait_sd)
    cor.estimate[i, j] <- round(result.df$estimate[[1]], 2)
    cor.pval[i, j] <- round(result.df$p.value[[1]], 2)
    
    plot(means_loc_year$dif ~ means_loc_year$trait_sd)
  }
}

# write out results
cor_bound <- list()

for (i in 1:length(envs)) {
  cor_bound[[i]] <- cbind(cor.estimate[,i], cor.pval[,i])
}

cor_out <- do.call("cbind", cor_bound)
colnames(cor_out) <- rep(c("est", "p-value"), 4)
rownames(cor_out) <- trait

file <- paste(out_path, "/variable_parents.txt", sep = "")
write.table(cor_out, file, row.names = T, col.names = T, sep = "\t", quote = F)


#### Step 5: Is there an Effect of Maternal Parent? ####
# read in the backbone to get the NAM names and parent identity (C or D)
# so we can attach this info to th emmeans, which utilize the plantID3 col, 
# this was balanced for the sake of the linear models
backbone <- read.csv(paste(in_path, "/backbone.csv", sep = ""), header = T) %>%
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>%
  select(famID_plantID3, longID) %>%
  distinct()

# maternal_means
# create an output dataframe
maternal_out <- data.frame(matrix(NA, nrow = 10, ncol = 4))
colnames(maternal_out) <- envs

mean_c <- maternal_out
mean_d <- maternal_out

maternal_out_list <- replicate(3, maternal_out, simplify = FALSE)
mean_c_list <- replicate(3, maternal_out, simplify = FALSE)
mean_d_list <- replicate(3, maternal_out, simplify = FALSE)

maternal_out_list <- list()

for (i in 1:length(trait)) {
  files <- list.files(paste(out_path, "/", trait[i], sep = ""), pattern = "emmeans_genet", full.names = T)
  files <- files[! grepl("transformed", files)]
  tables <- lapply(files, read.table, header = T, fill = TRUE)
  for (j in 1:length(tables)) {
    tables[[j]] <- tables[[j]] %>% 
      mutate(loc = loc[j], 
             year = year[j], 
             famID_plantID3 = paste(famID, plantID3, sep = "_"))
  }
  
  # add in the longID and the maternal parent
  tables_all <- do.call("rbind", tables)
  colnames(tables_all)[3] <- c("emmean")
  tables_join <- left_join(tables_all, backbone, by = "famID_plantID3") %>%
    mutate(maternal_parent = substr(longID, 6, 6))
  
  # iterate through locations and years, test for difference between maternal parent means
  for (k in 1:length(envs)) {
    for (l in 1:length(famID)) {
      table_loc_year <- tables_join[tables_join$loc == loc[k] & tables_join$year == year[k] & tables_join$famID == famID[l],]

      # separate by maternal parent
      c <- table_loc_year[table_loc_year$maternal_parent == "C",]
      d <- table_loc_year[table_loc_year$maternal_parent == "D",]
      
      # conduct t-test and output result
      t <- t.test(c$emmean, d$emmean)
      p <- t$p.value
      
      mean_c[l, k] <- round(mean(c$emmean, na.rm = T), 2)
      mean_d[l, k] <- round(mean(d$emmean, na.rm = T), 2)
      
      maternal_out[l, k] <- round(p, 2)
    }
  }
  
  mean_c_list[[i]] <- mean_c
  mean_d_list[[i]] <- mean_d
  
  maternal_out_list[[i]] <- maternal_out
}

# write out a table of these results
bind <- list()
trait_maternal_means_list <- list()

for (j in 1:length(trait)) {
  for (i in 1:length(envs)) {
    bind[[i]] <- cbind(mean_c_list[[j]][,i], mean_d_list[[j]][,i], maternal_out_list[[j]][,i])
  }
  cols_bound <- do.call("cbind", bind)
  colnames(cols_bound) <- rep(c("C", "D", "p-value"), 4)
  rownames(cols_bound) <- famID
  trait_maternal_means_list[[j]] <- cols_bound
}


file = paste(out_path, "/", "maternal_parents.txt", sep = "")
for (i in 1:length(trait)) {
  cat(paste("\t", trait[i], sep = ""), file = file, append = T)
  cat("\n", file = file, append = T)
  write.table(trait_maternal_means_list[[i]], file,
              col.names = NA, row.names = T, sep = "\t", quote = F, append = T)
  cat("\n", file = file, append = T)
}



  