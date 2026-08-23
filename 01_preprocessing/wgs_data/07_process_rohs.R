---
title: 'Process FROH '
author: "B. K. Montero"
date: "`r Sys.Date()`"
output: html_document

---


## LOAD LIBRARIES

```{r libraries, message=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=45), paged.print=FALSE}

rm(list = ls())

# basics
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(ggplot2)

```


## WORKING DIRECTORIES
```{r source, message=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=45), paged.print=FALSE}

source("~/Dropbox/FringillaMicrobiome/05_Scripts/workingDirectories.R")
source("~/Dropbox/FringillaMicrobiome/05_Scripts/CleanScripts/graphics_param.R")

```


## LOAD DATA


```{r rohs, message=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=45), paged.print=FALSE}

# import the bcftools output 
roh_raw <- fread(paste0(dataDir, "/microbiome_subset_roh_rg.txt"), 
             col.names=c("state", "file_id", "chr", "start", "end", "length", "nsnp", "qual"))


# Froh is based on autosomes only:
# keep only assembled chromosomes (OY7407*) and remove Z chromosome (OY740727)
roh_autoscaf <- roh_raw %>%
  filter(grepl("^OY7407", chr)) %>%   # removes CAUPSF scaffolds
  filter(chr != "OY740727")            # removes Z chromosome

# Filter based on quality (minimum of 30), ROH length (minimum of 10kb) and number of SNPs within the ROH (minimum of 100) 
roh_clean <- subset(roh_autoscaf, qual > 30 & length >= 10000 & < 100000 & nsnp >= 100) 

# Calculate Froh (sum the total length of all ROHs for all individuals and divide it by the total autosomal genome length)
#Autosomal genome length: 1005935441 bp

froh <- roh_clean %>% group_by(file_id) %>% summarise(total_length_bp = sum(length))
froh$froh <- froh$total_length_bp/1005935441

#as.data.frame(froh)

# Calculate the number of ROHs per individual, regardless of their length
n_roh <- roh_clean %>% group_by(file_id) %>% count()

# Quick plot
metadata <- read_delim(paste0(dataDir, "/WGS_Fringilla_samples.tsv"), delim = ",")
head(as.data.frame(metadata))
metadata$file_id <- metadata$TipID_microbiome

roh_meta <- left_join(froh, metadata, by ="file_id")
ggplot(roh_meta, aes(froh)) + 
        geom_histogram(fill = "#0f86a9", colour = "black", alpha = 0.3) +
        facet_wrap(. ~ Species , scales = "free") +
        theme_bw()

# Write to file

write.csv(froh, paste0(dataDir, "/fringilla_microbiome_froh_100kb_1Mb.csv"), row.names=FALSE)



