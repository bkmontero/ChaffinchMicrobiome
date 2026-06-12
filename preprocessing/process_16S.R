#!/usr/bin/Rscript
########################################################
##### ------------    FRINGILLA 16S   ------------ #####
#####           		 Dada2			           #####
########################################################

#options(warn=-1)
rm(list = ls())

## ------ LIBRARIES ------


library(devtools)
library(dada2); packageVersion("dada2")
library(ShortRead)
library(fastqcr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(data.table)
library(stringr)
library(phyloseq)
library(Biostrings)
library(microbiome)
library(seqateurs)

## ------ WORKING DIRECTORIES ------
source("~/Dropbox/FringillaMicrobiome/05_Scripts/workingDirectories.R")

## -----------------------------------------------------------------------------
## ------ I. LOAD & CHECK DATA ------

## ------   1. metadata --------

fringilla <- read.csv(file.path(dataDir, "metadata_fringilla_16S.tsv"))
sample_report <- read.csv(file.path(dataDir, "sample_report_all.csv"))


path <- fastaDir 
list.files(path)


## -----------------------------------------------------------------------------
## ------ II. RUN DADA2 -----

## ------   1. select files with at least 500 reads and 1400 read length --------

select_samples_1400_500reads <- read.csv(paste0(dataDir, "/select_samples_1400_500reads.csv"))
#331 samples


# update path
original_path <- list.files(path)

# Define pattern
keep_samples <- paste(select_samples_1400_500reads$SampleID, sep = "", collapse = "_|")

str(keep_samples)

# Subset list 
update_path <- list.files(fastaDir, pattern = keep_samples, ignore.case = TRUE, full.names = TRUE)


## ------  2. Remove primers

# - Enter primers. 

FWD <- "AGAGTTTGATCMTGGCTCAG" # Loop 16S forward primer
REV <- "TACCTTGTTACGACTT" # Loop 16S reverse primer


#Remove the primers and any flanking sequence from the reads, and filter out reads that don’t contain both primers:

remove_primers <- file.path(path, "remove_primers", basename(update_path))

track_remove_primers <- removePrimers(update_path, remove_primers, FWD, rc(REV), max.mismatch = 4, verbose=TRUE)

## ------  3. Filtering and Trimming

filt <- file.path(path, "filteredTrimmed", basename(update_path))

track_filter_trim <- filterAndTrim(remove_primers, filt, maxEE=6, truncQ = 8, minLen=1400, maxLen=1600, verbose=TRUE, multithread = FALSE)

track_reads <- cbind(SampleFile=track_remove_primers[,1], raw=track_remove_primers[,2], primers=track_remove_primers[,3], filtered=track_filter_trim[,2])

## ------  3. Denoising

filtpath <- paste0(fastaDir,"/filteredTrimmed") # CHANGE ME to the directory containing your filtered fastq files
filts <- list.files(filtpath, pattern=".fq", full.names=TRUE) # CHANGE if different file extensions


## ------  4. Dereplication


#Dereplicating the data collapses together reads that encode the same sequence this ends up saving computational time in later stages. (see section 4 https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html)

filts_pl1 <- sort(list.files(filtpath, pattern="3792_sample_A", full.names = TRUE))
filts_pl2 <- sort(list.files(filtpath, pattern="3784_sample", full.names = TRUE))
filts_pl3 <- sort(list.files(filtpath, pattern="3785_sample", full.names = TRUE))
filts_pl4 <- sort(list.files(filtpath, pattern="3786_sample", full.names = TRUE))
filts_pl5 <- sort(list.files(filtpath, pattern="3787_sample", full.names = TRUE))



sample_names_filts_pl1 <- sapply(strsplit(basename(filts_pl1), "_contig"), `[`, 1) 
names(filts_pl1) <- sample_names_filts_pl1
sample_names_filts_pl2 <- sapply(strsplit(basename(filts_pl2), "_contig"), `[`, 1) 
names(filts_pl2) <- sample_names_filts_pl2
sample_names_filts_pl3 <- sapply(strsplit(basename(filts_pl3), "_contig"), `[`, 1) 
names(filts_pl3) <- sample_names_filts_pl3
sample_names_filts_pl4 <- sapply(strsplit(basename(filts_pl4), "_contig"), `[`, 1) 
names(filts_pl4) <- sample_names_filts_pl4
sample_names_filts_pl5 <- sapply(strsplit(basename(filts_pl5), "_contig"), `[`, 1) 
names(filts_pl5) <- sample_names_filts_pl5


derep_pl1 <- derepFastq(filts_pl1, # file paths to fastq files
                          n = 1e+06, # maximum number of reads to parse and dereplicate at any one time.
                          verbose = TRUE) # outputs final status of the dereplication

names(derep_pl1) <- sample_names_filts_pl1

derep_pl2 <- derepFastq(filts_pl2, n = 1e+06,verbose = TRUE) 
names(derep_pl2) <- sample_names_filts_pl2

derep_pl3 <- derepFastq(filts_pl3, n = 1e+06,verbose = TRUE) 
names(derep_pl3) <- sample_names_filts_pl3

derep_pl4 <- derepFastq(filts_pl4, n = 1e+06,verbose = TRUE) 
names(derep_pl4) <- sample_names_filts_pl4

derep_pl5 <- derepFastq(filts_pl5, n = 1e+06,verbose = TRUE) 
names(derep_pl5) <- sample_names_filts_pl5

## ------  5. Error rates
# Binned error rates troubleshoot: https://github.com/benjjneb/dada2/issues/1307


set.seed(56456)


# Option 1: Alters the weights and span in loess, also enforce monotonicity

loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# Option 2: Only enforce monotonicity

loessErrfun_mod2 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}



# Option 3: Only alter loess weights and also enforce monotonicity


loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


# Option 4: Alter loess function arguments (weights, span, and degree) also enforce monotonicity.

loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# keep track of how long this step takes to run
tic(msg = NULL, quiet = TRUE)

# check what this looks like
errF_1 <- learnErrors(
  derep_pl1,
  multithread = TRUE,
  randomize = TRUE,
  nbases = 1e08,
  MAX_CONSIST = 12, 
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

# option 1 
saveRDS(errF_1, paste0(outputDir, "/02_Dada2/err_pl1_mod1.rds"))
#err <- readRDS(paste0(outputDir, "/02_Dada2/err.rds"))

plotErrors(errF_1, nominalQ=TRUE)
ggsave(paste(outputDir, "01_Exploration/plotErrors_pl1_mod1.png", sep = "/"))

write.csv(errF_1, paste0(outputDir, "/02_Dada2/err_pl1_mod1.csv"))

# check what this looks like
errF_2 <- learnErrors(
  derep_pl1,
  multithread = TRUE,
  randomize = TRUE,
  nbases = 1e08,
  MAX_CONSIST = 12, 
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)


# option 2
saveRDS(errF_2, paste0(outputDir, "/02_Dada2/err_pl1_mod2.rds"))
#err <- readRDS(paste0(outputDir, "/02_Dada2/err.rds"))

plotErrors(errF_2, nominalQ=TRUE)
ggsave(paste(outputDir, "01_Exploration/plotErrors_pl1_mod2.png", sep = "/"))

write.csv(errF_2, paste0(outputDir, "/02_Dada2/err_pl1_mod2.csv"))


# check what this looks like
errF_3 <- learnErrors(
  derep_pl1,
  multithread = TRUE,
  randomize = TRUE,
  nbases = 1e08,
  MAX_CONSIST = 12, 
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

# option 3
saveRDS(errF_3, paste0(outputDir, "/02_Dada2/err_pl1_mod3.rds"))
#err <- readRDS(paste0(outputDir, "/02_Dada2/err.rds"))

plotErrors(errF_3, nominalQ=TRUE)
ggsave(paste(outputDir, "01_Exploration/plotErrors_pl1_mod3.png", sep = "/"))

write.csv(errF_3, paste0(outputDir, "/02_Dada2/err_pl1_mod3.csv"))


# check what this looks like
errF_4 <- learnErrors(
  derep_pl1,
  multithread = TRUE,
  randomize = TRUE,
  nbases = 1e08,
  MAX_CONSIST = 12, 
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

# SAFE and PLOT

# option 4
saveRDS(errF_4, paste0(outputDir, "/02_Dada2/err_pl1_mod4.rds"))
#err <- readRDS(paste0(outputDir, "/02_Dada2/err.rds"))

plotErrors(errF_4, nominalQ=TRUE)
ggsave(paste(outputDir, "01_Exploration/plotErrors_pl1_mod4.png", sep = "/"))

write.csv(errF_4, paste0(outputDir, "/02_Dada2/err_pl1_mod4.csv"))


# So, the total time to get the error models from the data is
errmodel <- toc(log = FALSE, quiet = TRUE)

time_errormodel <- errmodel$toc - errmodel$tic
time_errormodel/3600

# use model 4 for this dataset

## ------  6. run Dada2

#Dereplicating the data collapses together reads that encode the same sequence this ends up saving computational time in later stages. (see section 4 https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html)


err_pl1<-readRDS(paste0(outputDir, "/02_Dada2/err_pl1_mod4.rds")) 
err_pl2<-readRDS(paste0(outputDir, "/02_Dada2/err_pl2_mod4.rds")) 
err_pl3<-readRDS(paste0(outputDir, "/02_Dada2/err_pl3_mod4.rds")) 
err_pl4<-readRDS(paste0(outputDir, "/02_Dada2/err_pl4_mod4.rds")) 
err_pl5<-readRDS(paste0(outputDir, "/02_Dada2/err_pl5_mod4.rds")) 


# ------- PLATE 1
# Infer sequence variants per plate


dds_1 <- vector("list", length(sample_names_filts_pl1))
names(dds_1) <- sample_names_filts_pl1

for(sam in sample_names_filts_pl1) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts_pl1[[sam]])
  dds_1[[sam]] <- dada(derep, err=err_pl1, multithread=TRUE,  pool=TRUE)
}

# Construct sequence table and write to disk
seqtab_pl1 <- makeSequenceTable(dds_1)
saveRDS(seqtab_pl1,paste0(outputDir, "/02_Dada2/seqtab_pl1_mod4_pooled.rds")) # CHANGE ME to where you want sequence table saved


# ------- PLATE 2
# Infer sequence variants per plate


dds_2 <- vector("list", length(sample_names_filts_pl2))
names(dds_2) <- sample_names_filts_pl2

for(sam in sample_names_filts_pl2) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts_pl2[[sam]])
  dds_2[[sam]] <- dada(derep, err=err_pl2, multithread=TRUE,  pool=TRUE)
}

# Construct sequence table and write to disk
seqtab_pl2 <- makeSequenceTable(dds_2)
saveRDS(seqtab_pl2,paste0(outputDir, "/02_Dada2/seqtab_pl2_mod4_pooled.rds")) # CHANGE ME to where you want sequence table saved


# ------- PLATE 3
# Infer sequence variants per plate


dds_3 <- vector("list", length(sample_names_filts_pl3))
names(dds_3) <- sample_names_filts_pl3

for(sam in sample_names_filts_pl3) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts_pl3[[sam]])
  dds_3[[sam]] <- dada(derep, err=err_pl3, multithread=TRUE,  pool=TRUE)
}

# Construct sequence table and write to disk
seqtab_pl3 <- makeSequenceTable(dds_3)
saveRDS(seqtab_pl3,paste0(outputDir, "/02_Dada2/seqtab_pl3_mod4_pooled.rds")) # CHANGE ME to where you want sequence table saved


# ------- PLATE 4
# Infer sequence variants per plate


dds_4 <- vector("list", length(sample_names_filts_pl4))
names(dds_4) <- sample_names_filts_pl4

for(sam in sample_names_filts_pl4) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts_pl4[[sam]])
  dds_4[[sam]] <- dada(derep, err=err_pl4, multithread=TRUE,  pool=TRUE)
}

# Construct sequence table and write to disk
seqtab_pl4 <- makeSequenceTable(dds_4)
saveRDS(seqtab_pl4,paste0(outputDir, "/02_Dada2/seqtab_pl4_mod4_pooled.rds")) # CHANGE ME to where you want sequence table saved


# ------- PLATE 5
# Infer sequence variants per plate


dds_5 <- vector("list", length(sample_names_filts_pl5))
names(dds_5) <- sample_names_filts_pl5

for(sam in sample_names_filts_pl5) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts_pl5[[sam]])
  dds_5[[sam]] <- dada(derep, err=err_pl5, multithread=TRUE,  pool=TRUE)
}

# Construct sequence table and write to disk
seqtab_pl5 <- makeSequenceTable(dds_5)
saveRDS(seqtab_pl5,paste0(outputDir, "/02_Dada2/seqtab_pl5_mod4_pooled.rds")) # CHANGE ME to where you want sequence table saved


# Merge multiple runs (if necessary)
seqtab_pl1 <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_pl1_pooled.rds"))
seqtab_pl2 <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_pl2_pooled.rds"))
seqtab_pl3 <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_pl3_pooled.rds"))
seqtab_pl4 <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_pl4_pooled.rds"))
seqtab_pl5 <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_pl5_pooled.rds"))
seqtab_all_pool <- mergeSequenceTables(seqtab_pl1, seqtab_pl2, seqtab_pl3, seqtab_pl4, seqtab_pl5)


## ------  7. Chimera removal

seqtab <- removeBimeraDenovo(seqtab_all_pool, method="pooled", minFoldParentOverAbundance = 8, multithread=TRUE)

saveRDS(seqtab,paste0(outputDir, "/02_Dada2/seqtab_nochim_all_mod4_pool.rds")) 

# ------  8. Assign taxonomy

seqtab <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_nochim_all_mod4_pool.rds"))


# assign taxonomy in batches of 50,000
to_split <- seq(1, ncol(seqtab), by = 50000) #change to 50000 with full dataset
to_split2 <- c(to_split[2:length(to_split)]-1, ncol(seqtab))

taxtab = NULL


for(i in 1:length(to_split)){
  seqtab2 <- seqtab[, to_split[i]:to_split2[i]]
  taxtab2 <- assignTaxonomy(seqtab2, paste0(dataDir,"/GTDB_bac120_arc53_ssu_r214_fullTaxo.fa.gz"), multithread = TRUE, minBoot = 80)
  #if(!is.null(paste0(dataDir,"/silva_species_assignment_v138.fa.gz"))){taxtab2 <- addSpecies(taxtab2, refFasta = paste0(dataDir,"/silva_species_assignment_v138.fa.gz"), verbose = TRUE)}
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}


saveRDS(taxtab, paste0(outputDir, "/02_Dada2/taxtab_GTDB2023_fulldata.rds"))
write.csv(taxtab, paste0(outputDir, "/02_Dada2/taxtab_GTDB2023_fulldata.csv"))


## ------  9. Clean ASV table

# Load 16S dataset using GTDB database taxonomy assignment and create phyloseq object 
fringilla <- import_qiime_sample_data(mapfilename = paste0(dataDir, "/metadata_fringilla_16S.tsv"))

seqtab <- readRDS(paste0(outputDir, "/02_Dada2/seqtab_nochim_all_mod4_pool.rds"))
taxtab <- readRDS(paste0(outputDir, "/02_Dada2/taxtab_GTDB2023_fulldata.rds"))

ps_GTDB2023 <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(fringilla), 
               tax_table(taxtab))

# Add a DNAStringSet object 
sequences <- Biostrings::DNAStringSet(taxa_names(ps_GTDB2023))
names(sequences) <- taxa_names(ps_GTDB2023)
ps_GTDB2023 <- merge_phyloseq(ps_GTDB2023, sequences)

# Rename ASVs
taxa_names(ps_GTDB2023) <- paste0("ASV", seq(ntaxa(ps_GTDB2023)))

# Clean data

tax <- data.frame(tax_table(ps_GTDB2023))
tax %>% head

tax.clean <- data.frame(row.names = row.names(tax),Kingdom = str_replace(tax[,1], "d__",""),Phylum = str_replace(tax[,2], "p__",""),Class = str_replace(tax[,3], "c__",""),Order = str_replace(tax[,4], "o__",""),Family = str_replace(tax[,5], "f__",""),Genus = str_replace(tax[,6], "g__",""),Species = str_replace(tax[,7], "s__",""),stringsAsFactors = FALSE)

tax_table(ps_GTDB2023) <- as.matrix(tax.clean)

# add ASV and ASV_SEQ columns to taxonomy table
tax_table(ps_GTDB2023) <- cbind(tax_table(ps_GTDB2023), 
    rownames(tax_table(ps_GTDB2023)), names(sequences))

colnames(tax_table(ps_GTDB2023)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",  "ASV_ID", "ASV_SEQ")

# remove unclassified phyla, mitochondria and chloroplast
ps_GTDB2023_clean <- ps_GTDB2023 %>%
  subset_taxa(Kingdom == "Bacteria" & Phylum != "NA" & !Phylum %in% c("", "uncharacterized"))

ps_GTDB2023_clean <- ps_GTDB2023_clean %>% subset_taxa(Order != "Chloroplast" & Order != "Mitochondria" & Family != "Mitochondria" & Family != "Chloroplast" & Genus != "Chloroplast" & Genus != "Mitochondria"| is.na(Order) | is.na(Family) | is.na(Genus))

#filter ASV with a sequencing depth of less than 10 on all samples
ps_GTDB2023_clean <- prune_taxa(taxa_sums(ps_GTDB2023_clean) > 10, ps_GTDB2023_clean)
ps_GTDB2023_clean 

# Add sequencing depth column 
sample_data(ps_GTDB2023_clean)$Seq_depth <- sample_sums(ps_GTDB2023_clean)

# Update metadata Fringilla project
metadata <- data.frame(sample_data(ps_GTDB2023_clean))

df <-  data.frame(sample_data(ps_GTDB2023_clean)) %>%  select("sampleid", "Lab_Code", "Plate", "Seq_depth", "Ring","Genus" , "Species", "Subspecies", "Species_Label", "Month", "Year", "Locality", "Habitat", "Xdec", "Ydec", "Archipelago_Mainland", "Island_Mainland",  "Area_km2", "Distance_mainland", "Age_ma", "Age", "Sex", , "Weight", "Wing", "X8_P", "Tail", "Nar", "P_Cr", "Height", "Width", "Cran", "Tars")

sample_data(ps_GTDB2023_clean) <- df

# Replace sampleid with Lab_Code to match sample names across datasets 

sample_names(ps_GTDB2023_clean) <- paste0("X",sample_data(ps_GTDB2023_clean)$Lab_Code)

# create Age (juvenile and adult) category
sample_data(ps_GTDB2023_clean)$Age_cat <- ifelse(sample_data(ps_GTDB2023_clean)$Age > 3, "adult", "fledging")

# recode Sex covariate 
sample_data(ps_GTDB2023_clean)$Sex <- ifelse(sample_data(ps_GTDB2023_clean)$Sex == 1, "male", "female")


# Fasta file for all ASVs in the filtered data set

ps_to_fasta(ps_GTDB2023_clean, out.file = paste0(dataDir, "/ps_clean_asv_GTDB2023.fasta"), width = 1800, seqnames = "ASV_ID")

# align in cesga

#scp -C csddepla@ft3.cesga.es:/mnt/netapp1/Store_csddepla/FringillaProject/02_Data/ps_clean_asv_GTDB2023_aln.fasta /home/fringilla/Dropbox/FringillaMicrobiome/07_Data/

# Save phyloseq object, sequence and taxonomy tables

saveRDS(ps_GTDB2023_clean, paste0(dataDir, "/ps_clean_16S_GTDB2023.rds"))
