# Environmental and genetic drivers differentially shaped resident and transient gut microbiota during the Macaronesian chaffinch radiation

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20126521-blue)](https://doi.org/10.5281/zenodo.20126521) ![Code License](https://img.shields.io/badge/code%20license-MIT-green) ![Data License](https://img.shields.io/badge/data%20license-CC%20BY%204.0-blue)

## Description
This repository contains the data-processing and analysis code for a study of gut microbiome assembly across the Macaronesian chaffinch (Fringilla) radiation. It integrates full-length long-read 16S rRNA microbiome profiling, rbcL and COI diet metabarcoding, and host whole-genome sequencing from 302 birds representing six Fringilla species across 11 populations in the Azores, Madeira, the Canary Islands and the Iberian mainland, using this replicated island system to separate spatial from host-evolutionary drivers of microbiome structure. The scripts cover the full workflow, from raw amplicon processing and OTU clustering through to the community-assembly tests, co-phylogenetic analyses and figures presented in the manuscript.

## Abstract
Disentangling whether gut microbiome diversification reflects host evolutionary history or spatial processes remains challenging because these factors are often confounded. The Macaronesian chaffinch (Fringilla) radiation provides a replicated natural experiment to separate these effects. Using long-read 16S sequencing, diet metabarcoding, and whole-genome data, we evaluated competing models of community assembly, testing host-mediated dispersal versus neutral assembly, and environmental filtering versus host selection. We show that microbiome alpha diversity across host species was largely shaped by island size and isolation, aligning with Island Biogeography Theory, rather than host colonisation history. This indicates that microbial dispersal occurs independently of their hosts. However, spatial and host factors influenced different aspects of the microbiome: geography and diet primarily predicted which taxa were present, whereas host genetic relatedness and heterozygosity were associated with the abundance of established residents. Notably, we show that inbreeding compromises a host's ability to regulate microbial abundance, leading to highly variable microbiomes. While we detected phylosymbiosis (closely related hosts harboring similar microbiomes), we found no evidence of co-diversification between birds and microbes. Instead, associations likely arise from convergent ecological filtering, where shared host traits filter similar bacteria from the local environment . Together, our results show that avian microbiomes are shaped by distinct spatial and host-mediated processes, and that separating microbial presence from abundance is important for understanding the mechanisms driving microbiome assembly.

## Instructions
Scripts are numbered by stage below; run them in this order, as later
  scripts read `.rds` objects written by earlier ones.

  **Requirements:** R ≥ 4.3 with Bioconductor (`dada2`, `phyloseq`, `DECIPHER`,
  `Biostrings`, `ANCOMBC`, `ggtree`, `treeio`), plus CRAN packages including
  `vegan`, `ape`, `phytools`, `paco`, `betapart`, `glmmTMB`, `tidyverse`.
  External tools: cutadapt, QIIME 2 (2024.5), PLINK 2, bcftools/samtools,
  and eMPRess. Taxonomy assignment requires the GTDB SSU reference
  (`bac120_arc53_ssu_r214`).

  Download the data from Zenodo (10.5281/zenodo.20126521) into `data/`, or
  point `FRINGILLA_DATA` at wherever you unpacked it.

  ### 1. Preprocessing (`preprocessing/`)
  1. `process_16S.R` — filter, denoise and assign taxonomy to long-read 16S
  2. `process_rbcL.Rmd` / `process_COI.Rmd` — diet marker processing (needs cutadapt)
  3. `rarefy_16S_diet.rmd` — rarefaction of 16S and diet data
  4. `OTU_clustering_16S.Rmd` / `OTU_clustering_diet.Rmd` — 99% OTU clustering
     (16S tree building uses QIIME 2)
  5. `taxFix_16S.rmd` — taxonomic curation
  6. `vennSampleSize.rmd` — final overlapping sample sets

  ### 2. Host genomic data (`analysis/`)
  Run on the HPC cluster; both shell scripts contain cluster-specific paths
  that must be edited.
  1. `plink_PCA_noLD.sh` → `pca_WGS.rmd` — population structure
  2. `bcftools_rohs.sh` → `process_rohs.R` — runs of homozygosity and F_ROH
  3. `draw_host_mt_nuclear_tree.rmd` — host mitochondrial and nuclear trees

  ### 3. Analyses and figures (`analysis/`)
  - `nmds_16S.rmd`, `nmds_diet.rmd` — ordination and PERMANOVA
  - `MRM_mantel_tests_mtHostTree.rmd` — Mantel and MRM tests
  - `test_fst_froh.rmd`, `test_fst_heteroz_clean.rmd` — host genetic effects
  - `PNM_tests.rmd`, `test_neutrality_allTaxa.rmd` — neutral model fitting
  - `list_topFamilyFix_overallData.rmd` — core family subsets
  - `PACo_global_tests.rmd`, `PACo_core_families.rmd` — cophylogenetic tests
  - `prepare_eMPRess_input.rmd` — writes input for eMPRess reconciliation
  - `SharedTaxa.rmd`, `SharedTaxa_UniFrac.rmd` — shared taxa analyses
  - `Figure3_plots.Rmd`, `Figure4_plots.rmd`, `Figure5_plots.rmd`,
    `Suppl_triangle.Rmd` — manuscript figuresNo 

## Authors
- B. Karina Montero (0000-0003-4246-6004)
- Mark A. F. Gillingham

## Affiliations
- Biodiversity Research Institute, Consejo Superior de Investigaciones Científicas (CSIC) and Oviedo University–Principality of Asturias, University of Oviedo, Campus of Mieres, Mieres E-33600, Spain.
- Centre for the Synthesis and Analysis of Biodiversity, Foundation for Research on Biodiversity, Montpellier 34000, France
- Center for Evolutionary and Functional Ecology (CEFE), Centre National de la Recherche Scientifique (CNRS),  Montpellier 34000, France
- Department of Ornithology, Max Planck Institute for Biological Intelligence, Eberhard Gwinner Straße, 82319 Seewiesen, Germany

## Contact
📧 b.karina.montero@gmail.com

## Funding
- Spanish Ministry of Science, Innovation and Universities and the European Regional Development Fund (PGC2018-097575-B-I00; PID2022-140091NB-I00)
- EECG research award granted to B.K.M. and M.A.F.G by the American Genetic Association.

## Citation
If you use this work please cite it using the DOI(s) above.

> Montero B K, Gillingham M A F, Ravinet M, Illera J C (under review) Environmental and genetic drivers differentially shaped resident and transient gut microbiota during the Macaronesian chaffinch radiation. Nature Communications

## License
- **Code:** licensed under MIT.
- **Data:** licensed under CC BY 4.0.

## Directory Structure
```text
ChaffinchMicrobiome/
├── 01_preprocessing/
│   ├── 01_process_16S.R
│   ├── 02_process_rbcL.Rmd
│   ├── 03_process_COI.Rmd
│   ├── 04_rarefy_16S_diet.rmd
│   ├── 05_OTU_clustering_16S.Rmd
│   ├── 06_OTU_clustering_diet.Rmd
│   ├── 07_taxFix_16S.rmd
│   ├── 08_vennSampleSize.html
│   ├── 08_vennSampleSize.rmd
│   └── MacaronesianChaffinches.png
└── 02_analysis/
    ├── 09_plink_PCA_noLD.sh
    ├── 10_pca_WGS.html
    ├── 10_pca_WGS.rmd
    ├── 11_bcftools_rohs.sh
    ├── 12_process_rohs.R
    ├── 13_draw_host_mt_nuclear_tree.rmd
    ├── 14_nmds_16S.html
    ├── 14_nmds_16S.rmd
    ├── 15_nmds_diet.html
    ├── 15_nmds_diet.rmd
    ├── 16_SharedTaxa.html
    ├── 16_SharedTaxa.rmd
    ├── 17_SharedTaxa_UniFrac.rmd
    ├── 18_Suppl_triangle.html
    ├── 18_Suppl_triangle.Rmd
    ├── 19_MRM_mantel_tests_mtHostTree.html
    ├── 19_MRM_mantel_tests_mtHostTree.rmd
    ├── 20_test_fst_heteroz_clean.rmd
    ├── 21_test_fst_froh.html
    ├── 21_test_fst_froh.rmd
    ├── 22_PNM_tests.rmd
    ├── 23_test_neutrality_allTaxa.rmd
    ├── 24_list_topFamilyFix_overallData.rmd
    ├── 25_PACo_global_tests.rmd
    ├── 26_PACo_core_families.rmd
    ├── 27_prepare_eMPRess_input.rmd
    ├── 28_Figure3_plots.html
    ├── 28_Figure3_plots.Rmd
    ├── 29_Figure4_plots.html
    ├── 29_Figure4_plots.rmd
    └── 30_Figure5_plots.rmd
```

## Other Files

### `01_preprocessing/08_vennSampleSize.html`

Resolves duplicate and re-extracted samples, generates the final corrected phyloseq objects, and summarises sample sizes and overlap across the 16S, diet and WGS datasets (Figure 2C).

### `01_preprocessing/MacaronesianChaffinches.png`

Icon image for GitHub

### `02_analysis/09_plink_PCA_noLD.sh`

Runs a PLINK 2 principal component analysis on the filtered SNP set without LD-pruning, writing eigenvector and eigenvalue files. Cluster script (SLURM).

### `02_analysis/10_pca_WGS.html`

Reads the PLINK output and plots host population genetic structure from whole-genome data (Figure 2F).

### `02_analysis/11_bcftools_rohs.sh`

Calls runs of homozygosity with bcftools, extracts RG records, and computes autosomal genome length excluding the Z chromosome. Cluster script (SLURM).

### `02_analysis/14_nmds_16S.html`

NMDS ordination of gut microbiome composition (Bray-Curtis) with PERMANOVA across taxonomic levels, exporting the PERMANOVA results table (Figure 2D).

### `02_analysis/15_nmds_diet.html`

NMDS ordination of diet composition (Raup-Crick) and PERMANOVA testing the effects of sex and age (Figure 2E).

### `02_analysis/16_SharedTaxa.html`

Tests whether sympatric chaffinches share more microbial taxa, using pairwise Jaccard distances across taxonomic levels with permutation t-tests.

### `02_analysis/18_Suppl_triangle.html`

Partitions Bray-Curtis and Jaccard beta diversity into turnover and nestedness components across taxonomic levels and draws the triangle plots (supplementary).

### `02_analysis/19_MRM_mantel_tests_mtHostTree.html`

### `02_analysis/21_test_fst_froh.html`

Tests whether host inbreeding (F_ROH) and heterozygosity predict microbiome composition and among-individual variability across taxonomic levels.

### `02_analysis/28_Figure3_plots.html`

Tests colonisation history against island biogeography using alpha diversity boxplots, Faith’s PD, and GLMs from phylum to OTU99, and exports the GLM and biogeography tables (Figure 3).

### `02_analysis/29_Figure4_plots.html`

Draws variance-partitioning Euler and bar plots alongside F_ST and F_ROH relationships with microbiome distance, using population-level subsampling and geographic centroids (Figure 4).

## Code
Scripts should be run in the following order:


1. **`01_preprocessing/01_process_16S.R`**
   
   Processes raw synthetic long-read 16S rRNA data with DADA2: selects reads passing depth and length thresholds, removes primers, filters, denoises, merges, removes chimeras, assigns taxonomy against GTDB r214, and builds the initial phyloseq object.

2. **`01_preprocessing/02_process_rbcL.Rmd`**
   
   Processes raw rbcL amplicon reads (plant component of the diet): primer removal with cutadapt, filtering and trimming, DADA2 denoising and merging, taxonomic assignment, and phyloseq assembly.

3. **`01_preprocessing/03_process_COI.Rmd`**
   
   Processes raw COI amplicon reads (invertebrate component of the diet): primer removal with cutadapt, filtering and trimming, DADA2 denoising and merging, taxonomic assignment, and phyloseq assembly.

4. **`01_preprocessing/04_rarefy_16S_diet.rmd`**
   
   Rarefies the 16S and merged diet datasets to even sequencing depth, calculates alpha diversity, defines the overlapping microbiome/diet sample sets, and exports rarefied ASV sequences as FASTA.

5. **`01_preprocessing/05_OTU_clustering_16S.Rmd`**
   
   Clusters rarefied 16S ASVs into 99% OTUs (DECIPHER alignment and clustering), rebuilds phyloseq objects with alpha diversity, exports the OTU FASTA and imports the rooted phylogeny built in QIIME 2.

6. **`01_preprocessing/06_OTU_clustering_diet.Rmd`**
   
   Performs the equivalent 99% OTU clustering for the rbcL and COI diet datasets.

7. **`01_preprocessing/07_taxFix_16S.rmd`**
   
   Curates the 16S taxonomy by removing chloroplast sequences misassigned to Rhizobiaceae and correcting other misassignments, producing the _taxFix phyloseq objects used by all downstream analyses.

8. **`01_preprocessing/08_vennSampleSize.rmd`**
   
   Resolves duplicate and re-extracted samples, generates the final corrected phyloseq objects, and summarises sample sizes and overlap across the 16S, diet and WGS datasets (Figure 2C).

9. **`02_analysis/10_pca_WGS.rmd`**
   
   Reads the PLINK output and plots host population genetic structure from whole-genome data (Figure 2F).

10. **`02_analysis/12_process_rohs.R`**
   
   Filters ROH calls on quality, length and SNP count, then calculates per-individual inbreeding coefficients (F_ROH) over the autosomal genome and ROH counts per bird.

11. **`02_analysis/13_draw_host_mt_nuclear_tree.rmd`**
   
   Builds and plots the host nuclear and mitogenome phylogenies, rooted on brambling and subset to individuals with microbiome data, showing bootstrap support >=70 (Figure 2A).

12. **`02_analysis/14_nmds_16S.rmd`**
   
   NMDS ordination of gut microbiome composition (Bray-Curtis) with PERMANOVA across taxonomic levels, exporting the PERMANOVA results table (Figure 2D).

13. **`02_analysis/15_nmds_diet.rmd`**
   
   NMDS ordination of diet composition (Raup-Crick) and PERMANOVA testing the effects of sex and age (Figure 2E).

14. **`02_analysis/16_SharedTaxa.rmd`**
   
   Tests whether sympatric chaffinches share more microbial taxa, using pairwise Jaccard distances across taxonomic levels with permutation t-tests.

15. **`02_analysis/17_SharedTaxa_UniFrac.rmd`**
   
   Repeats the shared-taxa comparison with phylogenetically informed UniFrac distances, from OTU to phylum level.

16. **`02_analysis/18_Suppl_triangle.Rmd`**
   
   Partitions Bray-Curtis and Jaccard beta diversity into turnover and nestedness components across taxonomic levels and draws the triangle plots (supplementary).

17. **`02_analysis/19_MRM_mantel_tests_mtHostTree.rmd`**
   
   Constructs geographic, microbiome, diet and host phylogenetic distance matrices, then runs Mantel tests, multiple regression on distance matrices (MRM) and variance partitioning.

18. **`02_analysis/20_test_fst_heteroz_clean.rmd`**
   
   Tests associations between host genetic differentiation (F_ST) and heterozygosity and microbiome dissimilarity at each taxonomic level.

19. **`02_analysis/21_test_fst_froh.rmd`**
   
   Tests whether host inbreeding (F_ROH) and heterozygosity predict microbiome composition and among-individual variability across taxonomic levels.

20. **`02_analysis/22_PNM_tests.rmd`**
   
   Fits Sloan’s neutral community model to the overlap dataset and classifies OTUs as neutral or non-neutral, saving each subset for downstream tests.

21. **`02_analysis/23_test_neutrality_allTaxa.rmd`**
   
   Applies neutral model fitting to the full dataset and identifies core taxa shared across populations.

22. **`02_analysis/24_list_topFamilyFix_overallData.rmd`**
   
   Builds lists of phyloseq objects for the core bacterial families, agglomerated at OTU, species and genus level, for the family-level cophylogenetic tests.

23. **`02_analysis/25_PACo_global_tests.rmd`**
   
   Runs Procrustean Approach to Cophylogeny (PACo) tests across the whole community and the neutral and non-neutral OTU subsets, testing host-microbe phylogenetic congruence.

24. **`02_analysis/26_PACo_core_families.rmd`**
   
   Restricts the PACo tests to the core bacterial families, reporting goodness-of-fit and ProTest p-values per family.

25. **`02_analysis/27_prepare_eMPRess_input.rmd`**
   
   Merges samples by island, balances island representation, and writes the host and symbiont trees plus tip mappings for event-based reconciliation in eMPRess.

26. **`02_analysis/28_Figure3_plots.Rmd`**
   
   Tests colonisation history against island biogeography using alpha diversity boxplots, Faith’s PD, and GLMs from phylum to OTU99, and exports the GLM and biogeography tables (Figure 3).

27. **`02_analysis/29_Figure4_plots.rmd`**
   
   Draws variance-partitioning Euler and bar plots alongside F_ST and F_ROH relationships with microbiome distance, using population-level subsampling and geographic centroids (Figure 4).

28. **`02_analysis/30_Figure5_plots.rmd`**
   
   Assembles the neutral model panel and associated microbiome visualisations into the final multi-panel figure (Figure 5).

## R Environment
**R version:** 4.6.1

| Package | Version |
| :------ | :------ |
| `adespatial` | 0.3.29 |
| `ANCOMBC` | 2.14.0 |
| `ape` | 5.8.1 |
| `betapart` | 1.6.1 |
| `biohelper` | 0.0.24.0 |
| `biomeUtils` | 0.22 |
| `Biostrings` | 2.80.1 |
| `car` | 3.1.5 |
| `cowplot` | 1.2.0 |
| `dada2` | 1.40.0 |
| `data.table` | 1.18.4 |
| `DECIPHER` | 3.8.0 |
| `decontam` | 1.32.0 |
| `details` | 0.4.0 |
| `devtools` | 2.5.2 |
| `doParallel` | 1.0.17 |
| `dplyr` | 1.2.1 |
| `ecodist` | 2.1.3 |
| `effects` | 4.2.5 |
| `eulerr` | 7.1.0 |
| `fantaxtic` | 0.2.1 |
| `fastqcr` | 0.1.3 |
| `future` | 1.70.0 |
| `future.apply` | 1.20.2 |
| `gdata` | 3.0.1 |
| `geosphere` | 1.6.8 |
| `ggeffects` | 2.3.2 |
| `ggforce` | 0.5.0 |
| `ggh4x` | 0.3.1 |
| `ggmagnify` | 0.4.2 |
| `ggnewscale` | 0.5.2 |
| `ggpattern` | 1.3.1 |
| `ggplot2` | 4.0.3 |
| `ggpubr` | 1.0.0 |
| `ggrepel` | 0.9.8 |
| `ggstar` | 1.0.6 |
| `ggtern` | 4.0.0 |
| `ggtext` | 0.1.2 |
| `ggtree` | 4.2.0 |
| `ggtreeExtra` | 1.22.0 |
| `glmmTMB` | 1.1.14 |
| `glue` | 1.8.1 |
| `gplots` | 3.3.0 |
| `gridExtra` | 2.3.1 |
| `gtable` | 0.3.6 |
| `hilldiv` | 1.5.3 |
| `Hmisc` | 5.2.6 |
| `iNEXT` | 3.0.2 |
| `kableExtra` | 1.4.0 |
| `knitr` | 1.51 |
| `magrittr` | 2.0.5 |
| `MASS` | 7.3.66 |
| `metagMisc` | 0.6.0.9000 |
| `MicEco` | 0.10.0 |
| `microbiome` | 1.31.3 |
| `microViz` | 0.13.1 |
| `minpack.lm` | 1.2.4 |
| `mirlyn` | 1.4.2 |
| `MKinfer` | 1.3 |
| `modEvA` | 3.45 |
| `multcomp` | 1.4.30 |
| `MuMIn` | 1.48.19 |
| `paco` | 0.5.0 |
| `patchwork` | 1.3.2 |
| `phangorn` | 2.12.1 |
| `phylobase` | 0.8.12 |
| `phyloseq` | 1.56.0 |
| `phylosmith` | 1.0.8 |
| `phytools` | 2.5.2 |
| `plyr` | 1.8.9 |
| `png` | 0.1.9 |
| `primers` | 1.2.1 |
| `purrr` | 1.2.2 |
| `RColorBrewer` | 1.1.3 |
| `Rcpp` | 1.1.2 |
| `readxl` | 1.5.0 |
| `RInSp` | 1.2.5 |
| `rlist` | 0.4.6.2 |
| `rmdformats` | 1.0.4 |
| `RVAideMemoire` | 0.9.83.12 |
| `scales` | 1.4.0 |
| `seqateurs` | 0.0.0.9000 |
| `sessioninfo` | 1.2.3 |
| `sf` | 1.1.1 |
| `ShortRead` | 1.70.0 |
| `speedyseq` | 0.5.3.9021 |
| `stringr` | 1.6.0 |
| `TDbook` | 0.0.6 |
| `tibble` | 3.3.1 |
| `tidyr` | 1.3.2 |
| `tidytree` | 0.4.7 |
| `tidyverse` | 2.0.0 |
| `treedata.table` | 0.1.1 |
| `treedataverse` | 0.0.1 |
| `treeio` | 1.36.1 |
| `tsiR` | 0.4.3 |
| `tyRa` | 0.1.0 |
| `vegan` | 2.7.5 |
| `VennDiagram` | 1.8.2 |

---
*README generated with READMEBuilder on 15 August 2026.*