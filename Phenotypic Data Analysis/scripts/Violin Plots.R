# Project: IWG_NAM_Shattering_Threshability
# Analysis - Violin Plots
# Author: Kayla R. Altendorf
# Date: 11/24/2020

# violin plots of variation for shattering, threshability

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
pub_names <- c("Brittle Rachis (Ct)", "Floret Shattering (scale)", "Threshability (%)")

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



for (i in 1:3) {
  
  i = 1
  files <- list.files(paste(out_path, "/", trait[i], sep = ""), pattern = "emmeans_genet", full.names = T)
  files <- files[! grepl("transformed", files)]
  tables <- lapply(files, read.table, header = T, fill = TRUE)
  for (j in 1:length(tables)) {
    tables[[j]] <- tables[[j]] %>% 
      mutate(loc = loc[j], 
             year = year[j],
             famID_plantID3 = paste(famID, plantID3, sep = "_"))
    colnames(tables[[j]])[3] <- "mean"
  }
  pheno_data_all <- do.call("rbind", tables)
  
  # declare order based on increasing progeny means for stp 2017
  
  pheno_order <- pheno_data_all %>% 
    filter(loc == "STP" & year == "2017") %>% 
    group_by(famID) %>% 
    summarise(trait_mean = mean(mean)) %>% 
    arrange(trait_mean) %>%
    as.data.frame() %>%
    select(famID)
  
  # arrange the data to follow this order
  pheno_data_all$famID <- factor(pheno_data_all$famID, levels = pheno_order$famID)
  
  p1 <- ggplot(pheno_data_all, aes(x = factor(famID), y = mean)) + 
    geom_violin() +
    stat_summary(fun=mean, geom="point", size=4, color = "black") + 
    facet_grid(loc~year, scales = "free_y") +
    geom_point(data = donor_parent_means[[i]],aes(color = famID), size = 4)  +
    scale_color_manual(values = c("#c1d42f", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#770026",  "#68bbaf", "#D55E00", "#0b3d4c", "#CC79A7"))+ 
    geom_hline(data = common_parent_means[[i]], aes(yintercept = mean), linetype = "dotted", size =2, color = "#616365") + # common parent line
    theme_bw() + 
    ylab(paste(pub_names[i])) + 
    xlab("Family ID")  + 
    theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size=20, face = "bold"), 
          strip.text.x = element_text(size = 20),
          strip.background.x = element_rect(fill = "white", colour = NA), 
          strip.text.y = element_text(size = 20),
          strip.background.y = element_rect(fill = "white", colour = NA), 
          legend.position = "none")
  
plot_grid(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"), label_size = 20)
          