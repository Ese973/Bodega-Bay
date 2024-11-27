# Libraries 
ibrary(tidyverse)
library(qiime2R)
library(phyloseq)
library(ggpubr)
library(randomcoloR)
library(cowplot)
library(RColorBrewer)

# Importing and cleaning data
# Imported taxonomy, dada2 table, rooted tree, and metadata from local storage
taxonomy_16s <- qiime2R::read_qza("16s/QIIME2 Outputs/16s-rep-sequences-taxonomy.qza")

table_16s <- qiime2R::read_qza("16s/QIIME2 Outputs/16s-dada2-table.qza")

metadata_16s <- read.csv("16s/Raw Data/16s_metadata.csv",
                     sep = ",", header = T, row.names = 1, quote = "")

tree_16s <- qiime2R::read_qza("16s/QIIME2 Outputs/rooted-16s-tree.qza")

# Extracted taxonomy data from taxonomy file
taxonomy_info_16s <- taxonomy_16s$data

# Parsed the taxonomy into 6 different taxonomic levels
parse_taxonomy_16s <- taxonomy_info_16s %>% separate(Taxon, sep =";", into = c("V1","V2","V3","V4","V5","V6", "V7"))

# Made Feature IDs the row names and removed the Feature.ID column 
rownames(parse_taxonomy_16s) <- parse_taxonomy_16s$Feature.ID 

parse_taxonomy_16s$Feature.ID <- NULL 

parse_taxonomy_16s$Consensus <- NULL 

# Removed any digits or characters from the beginning of taxon names
# Empty cells were classified as NA
# Replaced all NAs with "Unassigned"
parse_taxonomy_16s[] <- lapply(parse_taxonomy_16s, function(x) gsub("D_\\d+__", "", x))

parse_taxonomy_16s[parse_taxonomy_16s==""] <- NA

parse_taxonomy_16s <- parse_taxonomy_16s %>% replace(is.na(.), "Unassigned")

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_16s <- phyloseq::tax_table(as.matrix(parse_taxonomy_16s))

OTUMAT_16s <- table_16s$data

OTU_16s <- otu_table(OTUMAT_16s, taxa_are_rows = TRUE, replace_na(0))

TREE_16s <- tree_16s$data

# Creating phyloseq object
phylo_16s <- phyloseq(OTU_16s, sample_data(metadata_16s), TAX_16s, TREE_16s)

phylo_16s <- prune_samples(sample_sums(phylo_16s) >= 1, phylo_16s)

# Subset to exclude Mitochondria and Chloroplasts
phylo_16s <- subset_taxa(phylo_16s, V4 != "Chloroplast")

phylo_16s <- subset_taxa(phylo_16s, V5 != "Mitochondria")

phylo_16s <- subset_taxa(phylo_16s, V1 != "Unassigned")

# Subset phyloseq by Habitat
# This removed any mocks, blanks, or controls that were left in the original phyloseq
# Normalized the phyloseq for relative abundance 
phylo_16s <- phylo_16s %>% subset_samples(Habitat %in% c("Mason's Marina", "Campbell Cove", "West Park"))

phylo_normalized_16s <- phylo_normalized_16s <- transform_sample_counts(phylo_16s, function(x) x/sum(x))

# Created data frame from the phyloseq object
df_phylo_16s <- psmelt(phylo_16s)

