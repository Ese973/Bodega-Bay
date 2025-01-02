# Libraries 
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(ggpubr)
library(randomcoloR)
library(cowplot)
library(RColorBrewer)

# Importing and cleaning data
# Imported taxonomy, dada2 table, rooted tree, and metadata from local storage
taxonomy_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-rep-sequences-taxonomy.qza") 

table_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-dada2-table.qza")

tree_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/rooted-18S-tree.qza") 

metadata_18s <- read.csv("Raw Data/seagrasses_metadata_habitat_location.csv",
                     sep = ",", header = T, row.names = 1, quote = "")

# Extracted taxonomy data from taxonomy file
taxonomy_info_18s <- taxonomy_18s$data

# Parsed the taxonomy into 22 different taxonomic levels
# V1 is Domain, V2 is Kingdom, V14 is Phylum, V15 is Class, V17 is Order, V21 is Family, and V22 is Genus
parse_taxonomy_18s <- taxonomy_info_18s %>% separate(Taxon, sep =";", into = c("V1","V2","V3","V4","V5","V6","V7","V8",
                                                                                       "V9","V10","V11","V12","V13","V14","V15","V16","V17",
                                                                                       "V18","V19","V20","V21","V22")) 
# Made Feature IDs the row names and removed the Feature.ID column 
rownames(parse_taxonomy_18s) <- parse_taxonomy_18s$Feature.ID 

parse_taxonomy_18s$Feature.ID <- NULL 

# Removed any digits or characters from the beginning of taxon names
# Empty cells were classified as NA
# Replaced all NAs with "Unassigned"
parse_taxonomy_18s[] <- lapply(parse_taxonomy_18s, function(x) gsub("D_\\d+__", "", x))

parse_taxonomy_18s[parse_taxonomy_18s==""] <- NA

parse_taxonomy_18s <- parse_taxonomy_18s %>% replace(is.na(.), "Unassigned")

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_18s <- phyloseq::tax_table(as.matrix(parse_taxonomy_18s))

OTUMAT_18s <- table_18s$data

OTU_18s <- otu_table(OTUMAT_18s, taxa_are_rows = TRUE, replace_na(0))

TREE_18s <- tree_18s$data

# Created phyloseq and removed any samples that had less than 1 read
phylo_18s <- phyloseq(OTU_18s, sample_data(metadata_18s), TAX_18s, TREE_18s) 

phylo_18s <- prune_samples(sample_sums(phylo_18s) >= 1, phylo_18s)

# Removed the group Charophyta (sea grasses)
phylo_18s <- subset_taxa(phylo_18s, V4 != "Charophyta")

# Subset phyloseq by Habitat
# This removed any mocks, blanks, or controls that were left in the original phyloseq
# Normalized the phyloseq for relative abundance 
phylo_18s <- phylo_18s %>% subset_samples(Habitat %in% c("Mason's Marina", "Campbell Cove", "West Park"))

phylo_normalized_18s <- transform_sample_counts(phylo_18s, function(x) x/sum(x))

# Subset the samples that were Raw Sediment sample type
# Normalized the Raw Sediment samples for relative abundance 
phylo_raw_18s <- phylo_18s %>% subset_samples(SampleType %in% c("RawSediment"))

phylo_normalized_raw_18s <- transform_sample_counts(phylo_raw_18s, function(x) x/sum(x))

# Subset the samples that were Ludox sample type
# Normalized the Ludox samples for relative abundance 
phylo_ludox_18s <- phylo_18s %>% subset_samples(SampleType %in% c("Ludox"))

phylo_normalized_ludox_18s <- transform_sample_counts(phylo_ludox_18s, function(x) x/sum(x))

# Subset the samples for only Nematodes from Ludox samples
phylo_nematoda_18s <- subset_taxa(phylo_ludox_18s, V14 %in% c("Nematoda"))

phylo_normalized_nematoda_18s <- transform_sample_counts(phylo_nematoda_18s, function(x) x/sum(x))




# Plotting for Relative Abundance
 
# Raw Sediment
# Agglomerated taxa down to V4 rank and made the phyloseq object a data frame
glom_raw_18s <- tax_glom(phylo_normalized_raw_18s, taxrank = "V4")
data_raw_18s <- psmelt(glom_raw_18s)

# Obtained any taxa in the V4 rank with less than 5% abundance into "Other"
data_raw_18s$V4 <- as.character(data_raw_18s$V4)
data_raw_18s$V4[data_raw_18s$Abundance < 0.05] <- "Other (< 5% abund)"  

# Created taxonomy bar graph
raw_taxonomy_bar_18s <- ggplot(data_raw_18s, aes(x = Sample, y = Abundance, fill = V4)) +
  theme_bw() +
  geom_bar(stat = "identity") + 
  scale_color_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5", "#C7EAE5", "#80CDC1","#35978F", "#01665E", "salmon")) + 
  scale_fill_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "salmon"),
                    name = "Raw Sediment Rank", breaks = c("Apicomplexa", "Cercozoa","Chlorophyta",
                                                                  "Dinoflagellata",
                                                                  "Labyrinthulomycetes","MAST-1", "Ochrophyta", "Opisthokonta",
                                                                  "Peronosporomycetes", "Unassigned", "Other (< 5% abund)")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  facet_nested(~Habitat + Location, scales = "free", drop = TRUE) 


# Ludox
# Repeated the same process that was done for the raw sediment method
glom_ludox_18s <- tax_glom(phylo_normalized_ludox_18s, taxrank = "V4")
data_ludox_18s <- psmelt(glom_ludox_18s)
data_ludox_18s$V4 <- as.character(data_ludox_18s$V4)
data_ludox_18s$V4[data_ludox_18s$Abundance < 0.05] <- "Other (< 5% abund)"  

ludox_taxonomy_bar_18s <- ggplot(data_ludox_18s, aes(x = Sample, y = Abundance, fill = V4)) +
  theme_bw() +
  geom_bar(stat = "identity") + 
  scale_color_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1","#35978F", "#01665E", "salmon")) + 
  scale_fill_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#C7EAE5", "#80CDC1", "#35978F", "#01665E", "salmon"),
                    name = "Ludox Rank", breaks = c("Apicomplexa", "Cercozoa","Chlorophyta",
                                                                  "Dinoflagellata","Florideophycidae","Ochrophyta",
                                                                  "Opisthokonta",
                                                                  "Peronosporomycetes", "Unassigned", "Other (< 5% abund)")) +
  facet_nested(~Habitat + Location, scales = "free", drop = TRUE) +
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

# Combining the graphs into  one plot
plot_grid(raw_taxonomy_bar_18s, ludox_taxonomy_bar_18s, labels = c('A', 'B'), nrow = 2, label_size = 12)