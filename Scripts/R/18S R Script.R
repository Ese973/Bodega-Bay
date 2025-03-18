##################### Libraries #### 
library(tidyverse) #
library(qiime2R) #
library(phyloseq) #
library(tibble)
library(dplyr)
library(ggpubr) #
library(ggplot2)
library(randomcoloR) #
library(cowplot) #
library(RColorBrewer)
library(readxl)
library(gt)
library(ggh4x)
library(gridExtra)
library(fantaxtic) # 
library(MicrobiotaProcess)
library(reltools)
library(Biostrings)

##################### Importing and cleaning data ####

# Imported taxonomy, dada2 table, rooted tree, and metadata from local storage
taxonomy_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-rep-sequences-taxonomy.qza") 

table_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-dada2-table.qza")

tree_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/rooted-18S-tree.qza") 

metadata_18s <- read.csv("Raw Data/seagrasses_metadata_habitat_location.csv",
                     sep = ",", header = T, row.names = 1, quote = "")

refseq <- qiime2R::read_qza("18s/QIIME2 Outputs/denoised-rep-sequences.qza") 

refseq_df <- refseq$data

physeq_seq <- refseq(refseq_df)



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

write.csv(parse_taxonomy_18s,"18s/Clean Data/taxonomy.csv")

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_18s <- phyloseq::tax_table(as.matrix(parse_taxonomy_18s))

OTUMAT_18s <- table_18s$data

OTU_18s <- otu_table(OTUMAT_18s, taxa_are_rows = TRUE, replace_na(0))

TREE_18s <- tree_18s$data

# Created phyloseq and removed any samples that had less than 1 read
phylo_18s <- phyloseq(OTU_18s, sample_data(metadata_18s), TAX_18s, TREE_18s) 

phylo_18s <- prune_samples(sample_sums(phylo_18s) >= 1, phylo_18s)

phylo_18s <- merge_phyloseq(phylo_18s, physeq_seq)

# Removed the group Charophyta (sea grasses)
phylo_18s <- subset_taxa(phylo_18s, V4 != "Charophyta")


##################### Using Decontam for filering out contaminants ####

sample_data(phylo_18s)$is.neg <- sample_data(phylo_16s)$Location == "NegtCtrl" | sample_data(phylo_18s)$Location == "Mock" |
  sample_data(phylo_18s)$Location == "Blank" # create a sample-variable for contaminants
phylo_contaminants_18s <- isContaminant(phylo_18s, method = "prevalence", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance
table(phylo_contaminants_18s$contaminant) # check number of ASVs that are contaminents

# Make phyloseq object of presence-absence in negative controls and true samples
phylo_contaminants.pa_18s <- transform_sample_counts(phylo_18s, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence

ps.pa.neg_18s <- prune_samples(sample_data(phylo_contaminants.pa_18s)$Location == "Control" | sample_data(phylo_contaminants.pa_18s)$Location == "Blank" | 
                                 sample_data(phylo_contaminants.pa_18s)$Location == "Mock",phylo_contaminants.pa_18s) # identify controls

ps.pa.pos_18s <- prune_samples(sample_data(phylo_contaminants.pa_18s)$Location != "Control" | sample_data(phylo_contaminants.pa_18s)$Location != "Blank" |
                                 sample_data(phylo_contaminants.pa_18s)$Location != "Mock", phylo_contaminants.pa_18s) # identify samples

df.pa_18s <- data.frame(pa.pos=taxa_sums(ps.pa.pos_18s), pa.neg=taxa_sums(ps.pa.neg_18s), contaminant=phylo_contaminants_18s$contaminant) # convert into a dataframe

# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa_18s, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

phylo_18s <- prune_taxa(!phylo_contaminants_18s$contaminant, phylo_18s) # remove ASVs identified as decontaminants from the dataset

df_phylo_18s <- psmelt(phylo_18s)


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

##################### Plotting for Relative Abundance ####

# Raw Sediment
# Agglomerated taxa down to V4 rank and made the phyloseq object a data fram
glom_raw_18s <- tax_glom(phylo_normalized_raw_18s, taxrank = "V4")
data_raw_18s <- psmelt(glom_raw_18s)

# Obtained any taxa in the V4 rank with less than 5% abundance into "Other"
data_raw_18s$V4 <- as.character(data_raw_18s$V4)
data_raw_18s$V4[data_raw_18s$Abundance < 0.05] <- "Other (< 5% abund)"  


# Created taxonomy bar graph
raw_taxonomy_bar_18s <- ggplot(data_raw_18s, aes(x = Sample, y = Abundance, fill = V4)) +
  theme_bw() +
  geom_bar(stat = "identity") + 
  scale_color_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5", "#C7EAE5", "#80CDC1","#35978F","gray", "gray45")) + 
  scale_fill_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "grey","gray45"),
                    name = "Raw Sediment Rank", breaks = c("Apicomplexa", "Cercozoa","Chlorophyta",
                                                                  "Dinoflagellata",
                                                                  "Labyrinthulomycetes","MAST-1", "Ochrophyta", "Opisthokonta",
                                                                  "Peronosporomycetes", "Unassigned", "Other (< 5% abund)")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) 

raw_taxonomy_bar_18s


# Ludox
# Repeated the same process that was done for the raw sediment method
glom_ludox_18s <- tax_glom(phylo_normalized_ludox_18s, taxrank = "V4")
data_ludox_18s <- psmelt(glom_ludox_18s)
data_ludox_18s$V4 <- as.character(data_ludox_18s$V4)
data_ludox_18s$V4[data_ludox_18s$Abundance < 0.05] <- "Other (< 5% abund)"  

# Reorder Habitat 

ludox_taxonomy_bar_18s <- ggplot(data_ludox_18s, aes(x = Sample, y = Abundance, fill = V4)) +
  theme_bw() +
  geom_bar(stat = "identity") + 
  scale_color_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5", "#C7EAE5", "#80CDC1","#35978F","gray", "grey45")) + 
  scale_fill_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5","#C7EAE5", "#80CDC1", "#35978F","gray", "grey45"),
                    name = "Ludox Rank", breaks = c("Apicomplexa", "Cercozoa","Chlorophyta",
                                                                  "Dinoflagellata","Florideophycidae","Incertae Sedis", "Ochrophyta",
                                                                  "Opisthokonta",
                                                                  "Peronosporomycetes", "Unassigned", "Other (< 5% abund)")) +
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

ludox_taxonomy_bar_18s

# Combining the graphs into  one plot
plot_grid(raw_taxonomy_bar_18s, ludox_taxonomy_bar_18s, labels = c('A', 'B'), nrow = 2, label_size = 12) 

# Nematodes % vs Others

glom_ludox_18s_v14 <- tax_glom(phylo_normalized_ludox_18s, taxrank = "V14")
#glom_ludox_18s_v14 <- subset_taxa(glom_ludox_18s_v14, V14 != "Unassigned")
data_ludox_18s_v14 <- psmelt(glom_ludox_18s_v14)
data_ludox_18s_v14$V14 <- as.character(data_ludox_18s_v14$V14)
#data_ludox_18s_v14[, 28][data_ludox_18s_v14[, 28] != "Nematoda" | "Unassigned"] <- "Other"

data_ludox_18s_v14[, 28][!(data_ludox_18s_v14[, 28] == "Nematoda" | data_ludox_18s_v14[, 28] == "Unassigned")] <- "Other"


ggplot(data_ludox_18s_v14, aes(x = Sample, y = Abundance, fill = V14)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) +
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) 


##################### Looking for specific nematode with known symbionts ####

symbionts_18s <- subset_taxa(phylo_18s, V22 == "Robbea" | V22 == "Astomonema" | V22 == "Eubostrichus")
symbionts_18s <- prune_samples(sample_sums(symbionts_18s) >= 1, symbionts_18s)
symbionts_df_18s <- psmelt(symbionts_18s)

symbionts_df_18s <- symbionts_df_18s[-(26:136),]
write.csv(symbionts_df_18s, "18s_sym.csv")

palette <- distinctColorPalette(8)

symbionts_bar_18s <- ggplot(symbionts_df_18s, aes(x = Sample, y = Abundance, fill = V22)) +
  geom_bar(stat = "identity", na.rm = TRUE) +
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette, name = "Genus")

symbionts_bar_18s

symbiont_data <- read_tsv("Raw Data/blast_hit.tsv")

symbionts_plot <- grid.arrange(symbionts_bar_18s, tableGrob(symbiont_data), nrow = 1)

##################### Beta Diversity ####

phylo_18s_raw_ord <- ordinate(phylo_normalized_raw_18s, "PCoA", "bray")
phylo_18s_ord_raw_plot <- plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord, type="samples", color="Habitat", shape="Location", title="18S Raw Sediment ASVs") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("#2166AC","#B2182B","#F4A582"),breaks = c("Mason's Marina", "West Park", "Campbell Cove")) +
  #ylim(-0.45,0.45) +
  #xlim(-0.6,0.6) +
  theme(legend.position = "none")

phylo_18s_ord_raw_plot
        
phylo_18s_nematoda_ord <- ordinate(phylo_normalized_nematoda_18s, "PCoA", "bray")
phylo_18s_ord_nematoda_plot <- plot_ordination(phylo_normalized_nematoda_18s, phylo_18s_nematoda_ord, type="samples", color="Habitat", shape="Location", title="18S Nematoda ASVs (Ludox)") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("#2166AC","#B2182B","#F4A582"),breaks = c("Mason's Marina", "West Park", "Campbell Cove")) +
  #ylim(-0.45,0.45) +
  #xlim(-0.6,0.6) +
  theme(legend.position = "none")

phylo_18s_ord_nematoda_plot

phylo_18s_ludox_ord <- ordinate(phylo_normalized_ludox_18s, "PCoA", "bray")
phylo_18s_ord_ludox_plot <- plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord, type="samples", color="Habitat", shape="Location", title="18S Ludox ASVs") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("#2166AC","#B2182B","#F4A582"),breaks = c("Mason's Marina", "West Park", "Campbell Cove")) +
  #ylim(-0.45,0.45) +
  #xlim(-0.6,0.6) +
  theme(legend.position = "none")

phylo_18s_ord_ludox_plot


phylo_18s_ord_ludox_plot_2 <- plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord, type="samples", color="Location", shape="Habitat", title="18S Ludox ASVs") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("saddlebrown","seagreen3"),breaks = c("Bare Sediment", "Sea Grass")) +
  #ylim(-0.45,0.45) +
  #xlim(-0.6,0.6) +
  theme(legend.position = "none")

phylo_18s_ord_ludox_plot_2

phylo_18s_ord_raw_plot_2 <- plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord, type="samples", color="Location", shape="Habitat", title="18S Raw ASVs") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("saddlebrown","seagreen3"),breaks = c("Bare Sediment", "Sea Grass")) +
  #ylim(-0.45,0.45) +
  #xlim(-0.6,0.6) +
  theme(legend.position = "none")

phylo_18s_ord_raw_plot_2


# PERMANOVA
# Calculate bray curtis distance matrix
bray_18s_raw <- phyloseq::distance(phylo_normalized_raw_18s, method = "bray")
bray_18s_ludox <- phyloseq::distance(phylo_normalized_ludox_18s, method = "bray")

# make a data frame from the sample_data
sampledf_18s_raw <- data.frame(sample_data(phylo_raw_18s))
sampledf_18s_ludox <- data.frame(sample_data(phylo_ludox_18s))

# Adonis test
adonis_results_18s_raw_habitat <- adonis2(bray_18s_raw ~ Habitat, data = sampledf_18s_raw)
adonis_results_18s_raw_habitat

adonis_results_18s_raw_location <- adonis2(bray_18s_raw ~ Location, data = sampledf_18s_raw)
adonis_results_18s_raw_location

adonis_results_18s_ludox_habitat <- adonis2(bray_18s_ludox ~ Habitat, data = sampledf_18s_ludox)
adonis_results_18s_ludox_habitat

adonis_results_18s_ludox_location <- adonis2(bray_18s_ludox ~ Location, data = sampledf_18s_ludox)
adonis_results_18s_ludox_location


# Homogeneity of dispersion test
beta_18s_raw_results_habitat <- betadisper(bray_18s_raw, sampledf_18s_raw$Habitat)
permutest(beta_18s_raw_results_habitat)

beta_18s_raw_results_location <- betadisper(bray_18s_raw, sampledf_18s_raw$Location)
permutest(beta_18s_raw_results_location)

beta_18s_ludox_results_habitat <- betadisper(bray_18s_ludox, sampledf_18s_ludox$Habitat)
permutest(beta_18s_ludox_results_habitat)

beta_18s_ludox_results_location <- betadisper(bray_18s_ludox, sampledf_18s_ludox$Location)
permutest(beta_18s_ludox_results_location)

  
TukeyHSD(beta_18s_raw_results_habitat)
##################### Nematode taxonomy bar charts ####


nematode_palette <- brewer.pal(n = 11, name = 'RdYlBu')
nematode_palette <- gsub("#4575B4", "grey", nematode_palette)
nematode_palette <- gsub("#313695", "grey45", nematode_palette)


top_order_taxa <- fantaxtic::top_taxa(phylo_normalized_nematoda_18s, n_taxa=10, tax_level = "V17", include_na_taxa = T)
order_bar <- plot_bar(top_order_taxa$ps_obj, fill = "V17") + 
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE)  + 
  theme_bw() +
  geom_bar(aes(color = V17, fill = V17), stat = "identity", position = "stack", show.legend = F) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  scale_color_manual(values = nematode_palette) + 
  scale_fill_manual(values = nematode_palette, name = "Order") 

order_bar

top_fam_taxa <- fantaxtic::top_taxa(phylo_normalized_nematoda_18s, n_taxa=10, tax_level = "V21", include_na_taxa = T)
family_bar <- plot_bar(top_fam_taxa$ps_obj, fill = "V21")+ facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE)  + 
  theme_bw() +
  geom_bar(aes(color = V21, fill = V21), stat = "identity", position = "stack", show.legend = F) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  scale_color_manual(values = nematode_palette, breaks = c('Anticomidae','Chromadoridae','Comesomatidae',
                                                  'Cyatholaimidae','Desmodoridae','Enchelidiidae',
                                                  'Linhomoeidae','Oncholaimidae','Xyalidae','Unassigned', 'Other')) +
  scale_fill_manual(values = nematode_palette, name = "Family", breaks = c('Anticomidae','Chromadoridae','Comesomatidae',
                                                                  'Cyatholaimidae','Desmodoridae','Enchelidiidae',
                                                                  'Linhomoeidae','Oncholaimidae','Xyalidae','Unassigned', 'Other'))  

family_bar

top_gen_taxa <- fantaxtic::top_taxa(phylo_normalized_nematoda_18s, n_taxa=10, tax_level = "V22", include_na_taxa = T)
genus_bar <- plot_bar(top_gen_taxa$ps_obj, fill="V22") + facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) + 
  theme_bw() +
  geom_bar(aes(color = V22, fill = V22), stat = "identity", position = "stack", show.legend = F) + 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))  + 
  scale_color_manual(values = nematode_palette, breaks = c('Anticoma','Calyptronema','Daptonema','
                                                  Desmolaimus','Metalinhomoeus','Ptycholaimellus',
                                                  'Sabatieria','Terschellingia','Viscosia','Unassigned','Other')) + 
  scale_fill_manual(values = nematode_palette, name = "Genus", breaks = c('Anticoma','Calyptronema','Daptonema',
                                                                 'Desmolaimus','Metalinhomoeus','Ptycholaimellus',
                                                                 'Sabatieria','Terschellingia','Viscosia','Unassigned','Other'))

genus_bar




##################### Alpha Div #######
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset

alpha_div_18S_ludox <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "InvSimpson"),
  "Location" = phyloseq::sample_data(phylo_ludox_18s)$Location,
  "Habitat" = phyloseq::sample_data(phylo_ludox_18s)$Habitat)

head(alpha_div_18S_ludox)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_18S_ludox <- alpha_div_18S_ludox %>%
  rename(Simpson = InvSimpson)
head(alpha_div_18S_ludox)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_ludox <- alpha_div_18S_ludox[, c( 5, 4, 1, 2, 3)]
head(alpha_div_18S_ludox)


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by location, column 2
kw_alpha_18S_ludox <- as.data.frame(sapply(3:5, function(x) kruskal.test(alpha_div_18S_ludox[,x],
                                                                       alpha_div_18S_ludox[,2])))

# Rename columns with the proper variable names
kw_alpha_18S_ludox <- kw_alpha_18S_ludox %>%
  rename(Observed = V1,
         Shannon = V2,
         Simpson = V3)

kw_alpha_18S_ludox <- t(kw_alpha_18S_ludox) # transpose
kw_alpha_18S_ludox <- as_tibble(kw_alpha_18S_ludox, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_ludox <- kw_alpha_18S_ludox[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_ludox) # checking object class


# Plot alpha diversity measures
# Change names
soil.labs <- c("Bare Sediment", "Sea Grass")
names(soil.labs) <- c("BS", "SG")


habitat.labs <- c("MM", "CC", "WP")
names(habitat.labs) <- c("Mason's Marina", "Campbell Cove", "West Park")

alpha_color_habitat <- c("#F4A582", "#B2182B", "#2166AC")

alpha_div_18S_ludox <- transform(alpha_div_18S_ludox, Habitat=factor(Habitat,levels=c("Campbell Cove","West Park","Mason's Marina")))

ad_18s_ludox_location <- alpha_div_18S_ludox %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(comparisons = location_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Location, scales = "free") +
  scale_color_manual(values = alpha_color_habitat) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins 
  scale_x_discrete(
    labels = c("CC", "WP", "MM"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(legend.position="none")  

ad_18s_ludox_location


ad_18s_ludox_habitat <- alpha_div_18S_ludox %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Location, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Location), height = 0, width = .2) +
  facet_nested(metric ~ Habitat, scales = "free") +
  scale_color_manual(values = alpha_color_location) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("BS", "SG"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")

ad_18s_ludox_habitat


#### Raw
alpha_div_18S_raw <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_raw_18s, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_raw_18s, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_raw_18s, measures = "InvSimpson"),
  "Location" = phyloseq::sample_data(phylo_raw_18s)$Location,
  "Habitat" = phyloseq::sample_data(phylo_raw_18s)$Habitat)

head(alpha_div_18S_raw)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_18S_raw <- alpha_div_18S_raw %>%
  rename(Simpson = InvSimpson)
head(alpha_div_18S_raw)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_raw <- alpha_div_18S_raw[, c( 5, 4, 1, 2, 3)]
head(alpha_div_18S_raw)


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by location, column 2
kw_alpha_18S_raw <- as.data.frame(sapply(3:5, function(x) kruskal.test(alpha_div_18S_raw[,x],
                                                                       alpha_div_18S_raw[,2])))

# Rename columns with the proper variable names
kw_alpha_18S_raw <- kw_alpha_18S_raw %>%
  rename(Observed = V1,
         Shannon = V2,
         Simpson = V3)

kw_alpha_18S_raw <- t(kw_alpha_18S_raw) # transpose
kw_alpha_18S_raw <- as_tibble(kw_alpha_18S_raw, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_raw <- kw_alpha_18S_raw[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_raw) # checking object class


# Plot alpha diversity measures
# Change names
soil.labs <- c("Bare Sediment", "Sea Grass")
names(soil.labs) <- c("BS", "SG")


habitat.labs <- c("CC", "WP", "MM")
names(habitat.labs) <- c("Mason's Marina", "Campbell Cove", "West Park")

alpha_color_habitat <- c("#F4A582", "#B2182B", "#2166AC")

alpha_div_18S_raw <- transform(alpha_div_18S_raw, Habitat=factor(Habitat,levels=c("Campbell Cove","West Park","Mason's Marina")))



ad_18s_raw_location <- alpha_div_18S_raw %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(comparisons = location_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Location, scales = "free") +
  scale_color_manual(values = alpha_color_habitat) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins +
  scale_x_discrete(
    labels = c("CC", "WP", "MM"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(legend.position="none") 

ad_18s_raw_location

ad_18s_raw_habitat <- alpha_div_18S_raw %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Location, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Location), height = 0, width = .2) +
  facet_nested(metric ~ Habitat, scales = "free") +
  scale_color_manual(values = alpha_color_location) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("BS", "SG"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")

ad_18s_raw_habitat
                                         




