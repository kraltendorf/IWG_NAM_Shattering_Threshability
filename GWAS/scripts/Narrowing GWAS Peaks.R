# Project: IWG_NAM_Shattering_Threshability
# Script 2 - Narrowing GWAS Peaks
# Author: Kayla R. Altendorf
# Date: 12/23/2020

# load required packages
library("sommer")
library('sommer')
library("vcfR")
library("dplyr")
library("ggplot2")
library("tibble")
library("rrBLUP")
library("tidyr")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"

# script we're on
folder <- c("/GWAS")

# create a vector of traits
traits <- c("rachis_breaks_mean", "floret_score_mean", "threshability")
year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")
env <- c("stp17", "stp18", "tli17", "tli18")

# set the base pair threshold - this was identified in the LD analysis
# LD decays below 0.2 at a median rate of 21 cM. Average distance bp distance per cM across chromosomes x 21 = 
bp_threshold <-  63093157


#### Step 1: Read in and Format Genotypic Data ####
# this is the imputed dataset which was used in GWAS
vcf <- read.vcfR(paste(dir, folder, "/data/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE)


# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode) to their sample names
# read in key
key <- read.table(paste(dir, folder, "/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes <- as.data.frame(t(genotypes)) %>% 
  rownames_to_column(var = "GATK_Sample") %>%
  left_join(., key, by = "GATK_Sample") %>%
  arrange(Sample) %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove any phasing information that might be present
genotypes[genotypes=="0|1"] <- "0/1"
genotypes[genotypes=="1|0"] <- "0/1"
genotypes[genotypes=="1|1"] <- "1/1"
genotypes[genotypes=="0|0"] <- "0/0"

# convert to marker matrix
genotypes_matrix <- genotypes # create a new dataframe for the output

for (i in 1:nrow(genotypes)) {
  gt <- unlist(genotypes[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 1
  gt1[gt == "0/1"] <- 0
  gt1[gt == "1/1"] <- -1
  gt1[gt == "./."] <- NA
  genotypes_matrix[i,] <- gt1 
}

genotype <- cbind(id_frame, genotypes_matrix)  # add the id_frame back
rownames(genotype) <- NULL # remove rownames


#### Step 2: Read in and Format Phenotypic Data ####
# read in backbone
backbone <- read.csv(paste(dir, "Phenotypic Data Analysis/Data/", "backbone.csv", sep = ""), header = T) %>% 
  dplyr::select(famID, parent, plantID, plantID3, longID) %>% 
  mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% 
  dplyr::select(-plantID3, -famID) %>% 
  distinct()


# read in emmeans files for rachis_breaks_mean
files <- list.files(paste(dir, "/Phenotypic Data Analysis/output/", traits[1], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
rachis_breaks_mean <- list()
for (j in 1:length(files)) {
  rachis_breaks_mean[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean)) %>%
    arrange(longID)
}

# read in emmeans files for floret_score_mean
files <- list.files(paste(dir, "/Phenotypic Data Analysis/output/", traits[2], sep = ""), pattern = "emmeans_genet", full.names = TRUE)
floret_score_mean <- list()
for (j in 1:length(files)) {
  floret_score_mean[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean)) %>%
    arrange(longID)
}

# read in emmeans files for threshability
files <- list.files(paste(dir, "/Phenotypic Data Analysis/output/", traits[3], sep = ""), pattern = "emmeans_genet_transformed", full.names = TRUE)
threshability <- list()
for (j in 1:length(files)) {
  threshability[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean)) %>%
    arrange(longID)
}

#### Step 3: Read in GWAS Results from GAPIT ####
rachis_breaks_mean_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[1], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
floret_score_mean_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[2], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)
threshability_files <- list.files(paste(dir, folder, "/output/GAPIT/", traits[3], sep = ""), pattern = "GWAS.Results.csv", full.names = TRUE)

# now iterate through each trait
# filter significant SNPs with -log10(p) 3.6 or p.value < 0.00025

# rachis_breaks_mean
rachis_breaks_mean_gwas <- list()
for (i in 1:length(rachis_breaks_mean_files)) {
  rachis_breaks_mean_gwas[[i]] <- read.csv(rachis_breaks_mean_files[i], header = T) %>% filter(P.value < 0.00025)
}

# floret_score_mean
floret_score_mean_gwas <- list()
for (i in 1:length(floret_score_mean_files)) {
  floret_score_mean_gwas[[i]] <- read.csv(floret_score_mean_files[i], header = T) %>% filter(P.value < 0.00025)
}

# threshability
threshability_gwas <- list()
for (i in 1:length(threshability_files)) {
  threshability_gwas[[i]] <- read.csv(threshability_files[i], header = T) %>% filter(P.value < 0.00025)
}

#### Step 3: Identify Priority SNPs ####
# identify SNPs that are commonly significant across traits and environments 
# sort by most significant pvalues
traits
gwas <- threshability_gwas

priority_snps_pval <- rbind(do.call("rbind", gwas)) %>%
  group_by(SNP) %>% 
  summarise(p_val_sum = sum(P.value)) %>%
  arrange(p_val_sum)

# sort by most frequent
priority_snps_n <- rbind(do.call("rbind", gwas)) %>%
  group_by(SNP) %>% 
  tally() %>% 
  arrange(-n) %>%
  filter(n >= 2) 

# merge togther
priority_snps <- left_join(priority_snps_n, priority_snps_pval, by = "SNP") %>%
  mutate(chrom = substr(SNP, 4, 5)) %>%
  arrange(chrom, -n, p_val_sum) 


#### Step 4: Run 'mmmer' in Sommer ####
# declare trait here - edit for each 
traits
phenotype <- threshability
gwas <- threshability_gwas

for (i in 1:4) { 
  
  # identify samples in common
  phen <- phenotype[[i]]$longID
  gen <- colnames(genotype)[-1:-5]
  common <- Reduce(intersect, list(phen, gen))
  
  # filter phenotype and genotype so they match 
  phenotype[[i]] <- phenotype[[i]] %>% filter(longID %in% common)
  genotype_common <- genotype[, colnames(genotype) %in% common]
  print(summary(phenotype[[i]]$longID == colnames(genotype_common))) 
  
  # create the genotype matrix
  genotype_matrix <- t(genotype_common)
  colnames(genotype_matrix) <- id_frame$ID
  genotype_matrix <- data.frame(apply(genotype_matrix, 2, function(x) as.numeric(as.character(x))))
  rownames(genotype_matrix) <- colnames(genotype_common)
  
  # extract markers to test
  snps <- gwas[[i]]$SNP
  
  # create phenotype_genotype dataset
  snps_to_test <- cbind(phenotype[[i]], genotype_matrix[,colnames(genotype_matrix) %in% snps])
  rownames(snps_to_test) <- NULL
  
  # calculate the A matrix
  A <- A.mat(genotype_matrix)
  
  # prepare model terms
  fixed <- as.formula(paste("emmean ~ ",  paste0(snps, collapse = " + "), sep = ""))
  random <- ~vs(longID, Gu=A)
  
  # run the model
  fit <- mmer(fixed = fixed, random = random, data = snps_to_test)
  
  # look at the anova 
  anova <- anova.mmer(fit)
  anova_snps <- anova %>% slice(3:nrow(anova) - 1) %>% # ignore first and last terms
    rownames_to_column("SNP") %>%
    mutate(chrom = substr(SNP, 4, 5)) %>%
    arrange(chrom)
  
  anova_snps_nest <- anova_snps %>% group_by(chrom) %>% nest()
  
  for (j in 1:length(anova_snps_nest$data)) {
    if (nrow(anova_snps_nest$data[[j]]) > 1) { # if there's more than one on a chromosome
      keep <- anova_snps_nest$data[[j]] %>% filter(`Pr(>F)` <= 0.00025) # keep the ones that are significant
      
      # but if none are significant to that point, then keep the most significant one
      if (nrow(keep) == 0) {
        keep <- anova_snps_nest$data[[j]] %>% arrange(`Pr(>F)`) %>% slice(1:1)
      }
      remove <- anova_snps_nest$data[[j]] %>% filter(`Pr(>F)` > 0.001) # and remove those that are not
      remove_priority <- remove %>% filter(SNP %in% priority_snps$SNP)
      
      if (length(remove_priority$SNP) > 1) {
        keep <- rbind(keep, remove_priority) # but if a priority SNP was not significant we'll keep it in for ease of reporting
      }
      
      # go through the keepers and test to see if any are priority SNPs
      for (k in 1:nrow(keep)) {
        if (keep$SNP[k] %in% priority_snps$SNP) { # if a SNP is in the priority SNP list
          priority_bp_pos <- as.numeric(substr(keep$SNP[k], 7, nchar(keep$SNP[k]))) # extract it's base pair position
          
          keep <- keep %>% mutate(bp_pos = as.numeric(substr(SNP, 7, nchar(SNP))), # calculate the base pair difference between the remaining SNPs
                                  bp_diff = abs(bp_pos - priority_bp_pos)) %>% # if they're less than the threshold (i.e. within the window)
            filter(! between(bp_diff, 1, bp_threshold)) # then remove them
        }
      }
      anova_snps_nest$data[[j]] <- anova_snps_nest$data[[j]] %>% filter(SNP %in% keep$SNP) # only keep the final SNPs
    }
  }
  
  # unnest the final snps
  anova_snps_final <- anova_snps_nest %>% unnest(cols = c(data))
  
  # re-run sommer for the final variance explained 
  # prepare model terms
  fixed <- as.formula(paste("emmean ~ ",  paste0(anova_snps_final$SNP, collapse = " + "), sep = ""))
  random <- ~vs(longID, Gu=A)
  
  # run the model
  fit <- mmer(fixed = fixed, random = random, data = snps_to_test)
  
  # look at the anova 
  anova <- anova.mmer(fit)
  
  # extract out variance explained
  fixed_beta <- fit$Beta$Estimate
  X <- model.matrix(fixed, snps_to_test)
  var_fixed <- var(matrix(data = fixed_beta, nrow = nrow(X), ncol = ncol(X), byrow = T) * X)
  snp_var <- diag(var_fixed)[-1]
  var <- snp_var / (snp_var + sum(unlist(fit$sigma)))
  var_percent <- as.data.frame(var*100) %>% rownames_to_column("SNP")
  colnames(var_percent)[2] <- "percent_variation_explained"
  
  
  gwas[[i]] <- gwas[[i]] %>% 
    filter(SNP %in% anova_snps_final$SNP) %>% 
    left_join(., var_percent, by = "SNP")

  print(paste(env[[i]], " is complete at ", Sys.time(), sep = ""))
}



## edit trait here
for (i in 1:length(gwas)) {
  gwas[[i]] <- gwas[[i]] %>% mutate(trait = "threshability",
                                    loc = loc[i], 
                                    year = year[i])
}

gwas_final <- do.call("rbind", gwas)

# for filtering by snps that are repeated at least 2 times
#keep <- gwas_final %>% group_by(SNP) %>% tally() %>% arrange(-n) #%>% filter(n >= 2) # must occur in more than 1 env
#gwas_final <- gwas_final %>% filter(SNP %in% keep$SNP)

View(gwas_final %>% group_by(SNP) %>% tally())

# change file name here for trait
write.csv(gwas_final, paste(dir, folder, "/output/final_gwas_threshability_all_sig_snps.csv", sep = ""), row.names = FALSE)



#### Step 5: Create a Table of Results ####
# read in the output - for SNPs detected in at least two environments
rachis_breaks_mean_output <- read.csv(paste(dir, folder, "/output/final_gwas_rachis_breaks_mean.csv", sep = ""), header = T)
floret_score_mean_output <- read.csv(paste(dir, folder, "/output/final_gwas_floret_score_mean.csv", sep = ""), header = T)
threshability_output <- read.csv(paste(dir, folder, "/output/final_gwas_threshability.csv", sep = ""), header = T)

# for ALL SNPS
rachis_breaks_mean_output <- read.csv(paste(dir, folder, "/output/final_gwas_rachis_breaks_mean_all_sig_snps.csv", sep = ""), header = T)
floret_score_mean_output <- read.csv(paste(dir, folder, "/output/final_gwas_floret_score_mean_all_sig_snps.csv", sep = ""), header = T)
threshability_output <- read.csv(paste(dir, folder, "/output/final_gwas_threshability_all_sig_snps.csv", sep = ""), header = T)

# bring in reference and alternate alleles
gwas_all <- rbind(rachis_breaks_mean_output, floret_score_mean_output, threshability_output)
colnames(id_frame)[3] <- "SNP"
gwas_all1 <- left_join(gwas_all, id_frame, by = "SNP") %>% 
  mutate(Alleles = paste(REF, ALT, sep = "/"))



# determine segregating famililes
# first create a vector of families
famID <- backbone %>% 
  mutate(famID = substr(famID_plantID3, 1, 5)) %>%
  dplyr::select(famID) %>%
  distinct()
famID <- as.vector(famID[,1])

# calculate segregation types for each family at each significant locus
output <- data.frame(matrix(NA, ncol=4, nrow=length(unique(gwas_all1$SNP))))
colnames(output) <- c("SNP", "het", "hom_alt", "hom_ref") # create an output dataframe/list
output_list <- replicate(10, output, simplify = FALSE) # one for each family 

# this SNP segregates at 0.5, but is not apparently segregating in any family
# Chr08_72972971

for (i in 1:length(famID)) {
    fam_genotype <- cbind(id_frame, genotype[, grepl(famID[i], names(genotype))]) %>% #iterating through families one at a time
      filter(SNP %in% gwas_all1$SNP) %>%
      arrange(SNP) # extract and arrange all the unique significant SNPs
  
  fam_id_frame <- fam_genotype %>% dplyr::select(CHROM:ALT) # extract the ID to bind back later  
  fam_genotype <- fam_genotype %>% dplyr::select(-CHROM:-ALT) # and the genotype data to analyze
  
  output_list[[i]]$SNP <- fam_id_frame$SNP
  output_list[[i]]$het <- rowSums(fam_genotype == "0") / ncol(fam_genotype)
  output_list[[i]]$hom_alt <- rowSums(fam_genotype == "1") / ncol(fam_genotype)
  output_list[[i]]$hom_ref <- rowSums(fam_genotype == "-1") / ncol(fam_genotype)
}

# now assess which families are segregating 
# create a dataframe for the output
seg_output <- data.frame(fam_id_frame)

for (i in 1:length(output_list)) {
  for (j in 1:nrow(output_list[[i]])) {
    seg_pattern <- output_list[[i]][j,] %>% dplyr::select(-SNP)
    
    if (length(seg_pattern[seg_pattern > 0.05]) > 2) { # if it's not segregating, one of the types would be 1, or 0.95 and 0.05, but if two are over 0.05, then we 
      output_list[[i]]$seg[j] <- substr(famID[i], 4, 5) # know it's segregating 
    }
    else (output_list[[i]]$seg[j] <- NA)
  }
  seg_output[,i] <- output_list[[i]]$seg
}

output_list
id_seg <- cbind(fam_id_frame, seg_output)
colnames(id_seg)[-1:-5] <- famID

# combine all the results
id_seg_frame <- id_seg %>% 
  unite("Segregating_Families", WGN07:WGN63, na.rm = TRUE, remove = FALSE, sep = ", ") %>% 
  dplyr::select(-WGN07:-WGN63)


# there's one QTL on chromsome 15 that's associated with floret socre mean but it technically does not segregate except only within family 39
# het: 0.975409836 hom_alt: 0.01639344 hom_ref: 0.008196721
# every other family has 100% hom_alt.
# just don't include the family seg pattern

only_once <- gwas_all1 %>% group_by(trait, SNP) %>% tally() %>% filter(n == 1)
only_once$SNP

gwas_all1 <- gwas_all1 %>% filter(SNP %in% only_once$SNP)

gwas_all1


# join it back with the gwas_all1 and arrange into a nice table
gwas_all2 <- left_join(gwas_all1, id_seg_frame, by = "SNP") %>%
  dplyr::select(trait, SNP, Alleles, loc, year, Segregating_Families, maf, P.value, effect, percent_variation_explained) %>%
  mutate(P.value = round(-log10(P.value), digits = 2),
         effect = round(effect, digits = 2),
         percent_variation_explained = as.numeric(round(percent_variation_explained, digits = 2)), 
         maf = round(maf, digits = 2)) %>%
  group_by(trait) %>%
  pivot_wider(names_from = c(loc, year), values_from = c("P.value", "effect", "percent_variation_explained")) %>%
  dplyr::select(trait, SNP, Alleles, maf, Segregating_Families, P.value_STP_2017, effect_STP_2017, percent_variation_explained_STP_2017,
         P.value_STP_2018, effect_STP_2018, percent_variation_explained_STP_2018,
         P.value_TLI_2017, effect_TLI_2017, percent_variation_explained_TLI_2017,
         P.value_TLI_2018, effect_TLI_2018, percent_variation_explained_TLI_2018) %>%
  separate(SNP, into = c("Chromosome", "Position"), sep = "_") %>%
  arrange(trait, Chromosome, as.numeric(Position)) %>%
  mutate(SNP = paste(Chromosome, Position, sep = "_")) %>%
  dplyr::select(-Chromosome, -Position) 


#### not including segregation pattern ####
gwas_all1 %>% group_by(trait, SNP) %>% tally()
head(gwas_all1)

# join it back with the gwas_all1 and arrange into a nice table
gwas_all2 <- gwas_all1 %>%
  dplyr::select(trait, SNP, Alleles, loc, year, maf, P.value, effect, percent_variation_explained) %>%
  mutate(percent_variation_explained = percent_variation_explained/100) %>%
  mutate(P.value = round(-log10(P.value), digits = 2),
         effect = round(effect, digits = 2),
         percent_variation_explained = as.numeric(round(percent_variation_explained, digits = 2)), 
         maf = round(maf, digits = 2)) %>%
  group_by(trait) %>%
  pivot_wider(names_from = c(loc, year), values_from = c("P.value", "effect", "percent_variation_explained")) %>%
  separate(SNP, into = c("Chromosome", "Position"), sep = "_") %>%
  dplyr::select(trait, Chromosome, Position, Alleles, maf, P.value_STP_2017, effect_STP_2017, percent_variation_explained_STP_2017,
                P.value_STP_2018, effect_STP_2018, percent_variation_explained_STP_2018,
                P.value_TLI_2017, effect_TLI_2017, percent_variation_explained_TLI_2017,
                P.value_TLI_2018, effect_TLI_2018, percent_variation_explained_TLI_2018) %>%
  arrange(trait, Chromosome, as.numeric(Position)) %>%
  mutate(SNP = paste(Chromosome, Position, sep = "_")) %>%
  dplyr::select(-SNP)
 

gwas_all2$trait = factor(gwas_all2$trait, levels = traits)
gwas_all2 <- gwas_all2[order(gwas_all2$trait), ]

write.csv(gwas_all2, paste(dir, folder, "/output/gwas_final_table_only_once.csv", sep = ""), row.names = F)


head(gwas_all2)
gwas_all2 %>% group_by(Position) %>% tally() %>% arrange(-n)
### EXTRA ###

View(gwas_all2 %>% filter(Position == 72972971))

# calculate some summary statistics for the results section
gwas_all1 %>% select(SNP) %>% tally() # 45 total SNPs detected
gwas_all1 %>% group_by(trait) %>% tally() # 26 and 19 per trait
gwas_all1 %>% select(SNP) %>% group_by(SNP) %>% tally() %>% arrange(-n) # 7 detected in 2 environments, 2 in three or more environment/trait combinations
gwas_all1 %>% group_by(trait) %>% summarise(mean_pve = median(percent_variation_explained)) # 1.8 and 1.7
gwas_all1 %>% group_by(trait) %>% summarise(mean_effect = (median(abs(effect)))) # 0.324, 0.0647
gwas_all1 %>% select(SNP, Chromosome) %>% distinct() %>% group_by(Chromosome) %>% tally() %>% arrange(-n) # most detected on 5, 2, 6, 21, 9 and 16
gwas_all1 %>% group_by(loc, year) %>% tally() %>% arrange(-n)

# looking at the pddh1 locus
gwas_all1 %>% filter(SNP == "Chr06_27073072") %>% summarise(mean_pve = mean(percent_variation_explained))
# looking at the QTL on Chr 17
gwas_all1 %>% filter(Chromosome == 17)

# iterate through them
for (i in 1:length(snps)) {
  locus <- genotype_data %>% filter(rs == snps[i])
  
  alleles[[i]] <- locus$alleles
  loc[[i]] <- locus$rs
  
  # extract one family at a time
  for (j in 1:length(famID)) {
    fam <- locus[, grepl(famID[j], colnames(locus))]
    output_list[[1]][i, j] <- rowSums(fam == "0") / ncol(fam)
    output_list[[2]][i, j] <- rowSums(fam == "1") / ncol(fam)
    output_list[[3]][i, j] <- rowSums(fam == "-1") / ncol(fam)
    output_list[[4]][i, j] <- rowSums(is.na(fam)) / ncol(fam)
  }
}


##### Step 6: these are some extra bits of code to use for calculating LD if necessary ####   
# how about LD
check <- keep$SNP
ld_test <- genotype_matrix[, names(genotype_matrix) %in% check]
cor_out <- cor(ld_test)^2

# change the diagonals to NA
diag(cor_out) <- 0



