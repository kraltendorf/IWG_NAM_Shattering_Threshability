# Project: IWG_NAM_Shattering_Threshability
# Analysis - GWAS
# Author: Kayla R. Altendorf
# Date: 12/1/2020

# load packages
library("sommer")
library("vcfR")
library("dplyr")
library("tibble")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler") #this library is already installed in R
library("scatterplot3d")
library("cowplot")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c("multtest"))

# download gapit
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# location of the github directory
dir <- "/users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github/"
# folder we're on
folder <- "/GWAS"

#### Step 1: Load in Emmeans ####
# create a vector of traits

traits <- c("rachis_breaks_mean", "floret_score_mean", "threshability")
year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")
env <- c("stp17", "stp18", "tli17", "tli18")


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
           ! is.na(emmean))
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
           ! is.na(emmean))
}

# read in emmeans files for threshability
files <- list.files(paste(dir, "/Phenotypic Data Analysis/output/", traits[4], sep = ""), pattern = "emmeans_genet_transformed_scale", full.names = TRUE)
threshability <- list()
for (j in 1:length(files)) {
  threshability[[j]] <- read.table(files[j], header = T) %>% 
    mutate(year = year[j], loc = loc[j]) %>% 
    mutate(famID_plantID3 = paste(famID, plantID3, sep = "_")) %>% # create a column for joining
    left_join(backbone, ., by = "famID_plantID3") %>%
    dplyr::select(longID, emmean) %>%
    filter(! is.na(longID), 
           ! is.na(emmean))
}

##### Step 2: Prepare SNP Data in Numeric Format #####
# prepare imputed genotype file in numeric format
vcf <- read.vcfR(paste(dir, folder, "/data/NAM_GATK_imputed.vcf", sep = ""), convertNA = TRUE)

# extract genotypes and id_frame
genotypes <- extract.gt(vcf, convertNA = FALSE)
id_frame <- as.data.frame(vcf@fix[,1:5]) # extract a dataframe with SNP ids -- this will come in handy later
id_frame <- id_frame %>% mutate(ID = paste(CHROM, "_", POS, sep = ""))

head(id_frame)
test <- filter(id_frame, CHROM == "Chr21") %>%
  mutate(POS = as.numeric(POS))

hist(test$POS)

View(sommer::manhattan())

-log10(0.001)


View(manhattan)# update sample names from their ID in variant calling (e.g. flowcell, lane, barcode)
# to their sample names
# read in key
key <- read.table(paste(dir, folder, "/data/new_key.txt", sep = ""), header = T) %>% 
  mutate(GATK_Sample = paste(Flowcell, "_", Lane, "_", Barcode_ID, "_", sep = "")) %>%
  dplyr::select(GATK_Sample, Sample)

# prepare genotypes
genotypes1 <- t(genotypes) # transpose
genotypes2 <- as.data.frame(genotypes1) %>% rownames_to_column(var = "GATK_Sample")

genotypes3 <- left_join(genotypes2, key, by = "GATK_Sample") %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-GATK_Sample) %>%
  t()

# change all | to / to remove phasing information, if any
genotypes3[genotypes3=="0|1"] <- "0/1"
genotypes3[genotypes3=="1|0"] <- "0/1"
genotypes3[genotypes3=="1|1"] <- "1/1"
genotypes3[genotypes3=="0|0"] <- "0/0"

genotypes_numeric <- genotypes3

for (i in 1:nrow(genotypes3)) {
  gt <- unlist(genotypes3[i,])
  gt1 <- gt
  gt1[gt == "0/0"] <- 0
  gt1[gt == "0/1"] <- 1
  gt1[gt == "1/1"] <- 2
  gt1[gt == "./."] <- NA
  genotypes_numeric[i,] <- gt1 
}

genotypes_numeric <- as.data.frame(t(genotypes_numeric)) %>% rownames_to_column(var = "taxa")

# extract snp information as a separate file
mdp_SNP_information <- id_frame %>% 
  dplyr::select(ID, CHROM, POS) %>%
  dplyr::rename(Name = ID, 
                Chromosome = CHROM, 
                Position = POS) %>%
  mutate(Chromosome = as.numeric(substr(Chromosome, 4, 5)))

# write out result
write.table(genotypes_numeric, paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(mdp_SNP_information, paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

# read it back in 
myGD <- read.table(paste(dir, folder, "/output/GATK_NAM_snp_matrix_imputed.table.txt", sep = ""), head = T) 
myGM <- read.table(paste(dir, folder, "/output/mdp_SNP_information.txt", sep = ""), head = T) 


#### Step 3: Run GAPIT ####
# edit here only
traits # to reference trait order
trait <- floret_score_mean
trait_name <- c("floret_score_mean")

# set directory
dir.create(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))
setwd(paste(dir, folder, "/output/GAPIT/", trait_name, sep = ""))

# run

for (i in 1:length(trait)) {
  phenotype <- trait[[i]]
  
  # choose only samples that are in common
  P <- phenotype$longID
  G <- myGD$taxa
  common <- Reduce(intersect, list(G, P))
  
  # filter phenotype data
  P2 <- filter(phenotype, longID %in% common) %>% droplevels() %>% arrange(longID)
  
  myGD2 <- myGD %>% 
    filter(taxa %in% common) %>%
    arrange(taxa)
  
  # make sure sample names align
  print(summary(myGD2$taxa == P2$longID))
  
  # rename phenotype header
  colnames(P2)[1] <- "Taxa"
  
  # rename phenotype to correct trait
  colnames(P2)[2] <- paste(trait_name, env[i], sep = "_")
  
  # gapit command
  myGAPIT <- GAPIT(
    Y=P2,
    GD=myGD2,
    GM=myGM,
    SNP.MAF=0.005, # inserting a MAF here since imputation introduced a few low freq alleles post vcftools filtering
    PCA.total = 0)
}



