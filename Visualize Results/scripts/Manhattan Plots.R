# Project: IWG_NAM_Shattering_Threshability
# Visualize Results
# Author: Kayla R. Altendorf
# Date: 03/02/21

# load packages
library('sommer')

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Shattering and Threshability/Scripts for Github"

# script we're on
folder <- c("/Visualize Results")
traits <- c("rachis_breaks_mean", "flroet_score_mean", "threshability")
traits <- c("emergence_percent", "anthesis_score")

year <- c("2017", "2018", "2017", "2018")
loc <- c("STP", "STP", "TLI", "TLI")

#### Step 1: Read in GWAS Results ####
files <- list.files(paste(dir, "GWAS/output/GAPIT/anthesis_score/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
anthesis_score <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/rachis_breaks_mean/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
rachis_breaks_mean <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/floret_score_mean/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
floret_score_mean <- lapply(files, read.csv)

files <- list.files(paste(dir, "GWAS/output/GAPIT/threshability/", sep = "/"), pattern = "GWAS.Results.csv", full.names = T)
threshability <- lapply(files, read.csv)

#### Step 2: Make Manhattan Plots ####
# edit trait and save each plot separately
trait <- emergence_percent
  
for (i in 1:length(trait)) {
  
  trait[[i]] <- trait[[i]] %>% dplyr::select(Chromosome, Position, P.value)
  colnames(trait[[i]]) <- c("Chrom", "Position", "p.val")
  
  trait[[i]] <- trait[[i]] %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(Position)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(trait[[i]], ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, Position) %>%
  mutate( BPcum=Position+tot) %>%
  
  # add loc and year info to each plot
  mutate(loc = loc[i], 
         year = year[i])
}

dat <- do.call("rbind", trait)
axisdf = dat %>% group_by(Chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#Ready to make the plot using ggplot2:
  
  ggplot(dat, aes(x=BPcum, y=-log10(p.val))) +
    facet_grid(rows = loc ~ year) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  geom_hline(yintercept = 3.6, linetype = "dashed", color = "grey", size = 0.5) + 
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Chromosome") + 
    
  # Custom the theme:
  theme_bw() +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size=20, face = "bold"), 
          strip.text.x = element_text(size = 20),
          strip.background.x = element_rect(fill = "white", colour = NA), 
          strip.text.y = element_text(size = 20),
          strip.background.y = element_rect(fill = "white", colour = NA), 
          legend.position = "none", 
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
          )


    