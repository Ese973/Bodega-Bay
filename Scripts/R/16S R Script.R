###################### Libraries ####
library(tidyverse) 
library(qiime2R) #Integrating QIIME2 and R for data visualization and analysis using qiime2R.
library(phyloseq) #Handling and analysis of high-throughput microbiome census data.
library(ggh4x) #ggh4x package is a ggplot2 extension package.
library(cowplot) #Simple add-on to ggplot2; provides additional features to improve graphic quality.  
library(zCompositions)
library(RColorBrewer) #Provides color schemes for maps/graphics designed by Cynthia Brewer.
library(colorRamp2)
library(randomcoloR)
library(ComplexHeatmap) # Package for Complex heatmaps reveal patterns and correlations in multidimensional genomic data.
library(ggpubr) #Provides some easy-to-use functions to customize ggplot2 graphics.
library(ALDEx2)
library(decontam)
library(fantaxtic) # 
library(microbiome)
library(gridExtra)
library(microViz)
library(wesanderson)
library(vegan)


###################### Importing and cleaning data ####

# Imported taxonomy, dada2 table, rooted tree, and metadata from local storage
taxonomy_16s <- qiime2R::read_qza("16s/QIIME2 Outputs/updated-16s-rep-sequences-taxonomy.qza")

table_16s <- qiime2R::read_qza("16s/QIIME2 Outputs/16s-dada2-table.qza")

metadata_16s <- read.csv("16s/Raw Data/16s_metadata.csv",
                     sep = ",", header = T, row = 1, quote = "")

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

write.csv(parse_taxonomy_16s,"16s/Raw Data/taxonomy.csv")

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_16s <- phyloseq::tax_table(as.matrix(parse_taxonomy_16s))

OTUMAT_16s <- table_16s$data

OTU_16s <- otu_table(OTUMAT_16s, taxa_are_rows = TRUE, replace_na(0))

TREE_16s <- tree_16s$data

# Creating phyloseq object
phylo_16s <- phyloseq(OTU_16s, sample_data(metadata_16s), TAX_16s, TREE_16s)

phylo_16s <- prune_samples(sample_sums(phylo_16s) >= 1, phylo_16s)

to_remove <- c("16S-bodega-bay-fecal-experiment_CC.SG.3.2")

phylo_16s <- prune_samples(!(sample_names(phylo_16s) %in% to_remove), phylo_16s)

# Subset to exclude Mitochondria and Chloroplasts
phylo_16s <- subset_taxa(phylo_16s, V4 != "Chloroplast")

phylo_16s <- subset_taxa(phylo_16s, V5 != "Mitochondria")

phylo_16s <- subset_taxa(phylo_16s, V1 != "Unassigned")




###################### Using Decontam for filering out contaminants ####

sample_data(phylo_16s)$is.neg <- sample_data(phylo_16s)$Location == "Control" | sample_data(phylo_16s)$Location == "Blank" # create a sample-variable for contaminants
phylo_contaminants <- isContaminant(phylo_16s, method = "prevalence", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance
table(phylo_contaminants$contaminant) # check number of ASVs that are contaminents

# Make phyloseq object of presence-absence in negative controls and true samples
phylo_contaminants.pa <- transform_sample_counts(phylo_16s, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence
ps.pa.neg <- prune_samples(sample_data(phylo_contaminants.pa)$Location == "Control" | sample_data(phylo_contaminants.pa)$Location == "Blank", phylo_contaminants.pa) # identify controls
ps.pa.pos <- prune_samples(sample_data(phylo_contaminants.pa)$Location != "Control" | sample_data(phylo_contaminants.pa)$Location != "Blank", phylo_contaminants.pa) # identify samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=phylo_contaminants$contaminant) # convert into a dataframe

# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

phylo_16s <- prune_taxa(!phylo_contaminants$contaminant, phylo_16s) # remove ASVs identified as decontaminants from the dataset



# Subset phyloseq by Habitat
# This removed any mocks, blanks, or controls that were left in the original phyloseq
# Normalized the phyloseq for relative abundance 
phylo_16s <- phylo_16s %>% subset_samples(Habitat %in% c("Mason's Marina", "Campbell Cove", "West Park"))

phylo_normalized_16s <- phylo_normalized_16s <- transform_sample_counts(phylo_16s, function(x) x/sum(x))

# Created data frame from the phyloseq object
df_phylo_16s <- psmelt(phylo_16s)
df_normalized_16s <- psmelt(phylo_normalized_16s)


###################### Plotting for Relative Abundance ####

# Raw Sediment
# Agglomerated taxa down to V3 rank and made the phyloseq object a data fram
glom_16s <- tax_glom(phylo_normalized_16s, taxrank = "V3")
data_16s <- psmelt(glom_16s)

# Obtained any taxa in the V4 rank with less than 5% abundance into "Other"
data_16s$V3 <- as.character(data_16s$V3)
data_16s$V3[data_16s$Abundance < 0.05] <- "Other (< 5% abund)" 

# Reorder Habitat levels

data_16s <- transform(data_16s, Habitat=factor(Habitat,levels=c("Campbell Cove","West Park","Mason's Marina")))

# Created taxonomy bar graph
taxonomy_bar_16s <-  ggplot(data_16s, aes(x = Sample, y = Abundance, fill = V3)) +
  geom_bar(stat = "identity") +
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) +
  scale_color_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5","#C7EAE5", "#80CDC1", "#35978F","#003C30", "grey45")) + 
  scale_fill_manual(values = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3","#F5F5F5","#C7EAE5", "#80CDC1", "#35978F","#003C30", "grey45"),
                    name = "Bacteria Rank", breaks = c("Acidimicrobiia", "Alphaproteobacteria","Anaerolineae",
                                                    "Bacteroidia","Campylobacteria","Desulfobacteria",
                                                    "Desulfobulbia",
                                                    "Gammaproteobacteria", "Planctomycetes","Verrucomicrobiae", "Other (< 5% abund)")) +
 scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) 

taxonomy_bar_16s


###################### Beta Diversity ####

phylo_16s_ord <- ordinate(phylo_normalized_16s, "PCoA", "bray")

phylo_16s_ord_plot <- plot_ordination(phylo_normalized_16s, phylo_16s_ord, type="samples", color="Habitat", shape="Location", title="16s All ASVs") + 
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values = c("#2166AC","#B2182B","#F4A582"),breaks = c("Mason's Marina", "West Park", "Campbell Cove")) +
  #ylim(-0.45,0.45) +
 # xlim(-0.6,0.6) +
  theme(legend.position = "none")



phylo_16s_ord_plot


phylo_16s_ord_plot_2 <- plot_ordination(phylo_normalized_16s, phylo_16s_ord, type="samples", color="Location", shape="Habitat", title="16S All ASVs") + 
  scale_color_manual(values = c("saddlebrown","seagreen3"),breaks = c("Bare Sediment", "Sea Grass")) +
  #ylim(-0.45,0.45) +
  # xlim(-0.6,0.6) +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(aes(color =  Location, shape = Habitat), alpha = 0.7, size = 4) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_16s_ord_plot_2

# PERMANOVA
# Calculate bray curtis distance matrix
bray_16s <- phyloseq::distance(phylo_normalized_16s, method = "bray")

# make a data frame from the sample_data
sampledf_16s <- data.frame(sample_data(phylo_16s))

# Adonis test
adonis_results_16s_habitat <- adonis2(bray_16s ~ Habitat, data = sampledf_16s)
adonis_results_16s_habitat

adonis_results_16s_location <- adonis2(bray_16s ~ Location, data = sampledf_16s)
adonis_results_16s_location

# Homogeneity of dispersion test
beta_16s_results_habitat <- betadisper(bray_16s, sampledf_16s$Habitat)
permutest(beta_16s_results_habitat)

beta_16s_results_location <- betadisper(bray_16s, sampledf_16s$Location)
permutest(beta_16s_results_location)


###################### Looking for sulfur-utilizing symbionts ####

phylo_ss <- phylo_normalized_16s

#phylo_ss <- prune_samples(sample_sums(phylo_ss) >= 1, phylo_ss)

phylo_ss <- tax_glom(phylo_ss, taxrank = "V6")

phylo_ss <- tax_fix(phylo_ss, unknowns = c("Unassigned", "uncultured", "Unkown"))
 
df_ss <- psmelt(phylo_ss)

sulfur_symbionts_16s <- subset_taxa(phylo_ss, V6 == "Pseudoalteromonas" | 
                                      V6 == "Candidatus Thiodiazotropha" | 
                                      V6 =="Thiobacillus" |
                                      V6 == "Paracoccus" |
                                      V6 == "Pseudomonas" |
                                      V5 == "Chromatiaceae" |
                                      V6 == "Candidatus Thiosymbion" |
                                      V5 == "Ectothiorhodospiraceae")

sulfur_symbionts_16s <- prune_samples(sample_sums(sulfur_symbionts_16s) >= 1, sulfur_symbionts_16s)

sulfur_symbionts_df_16s <- psmelt(sulfur_symbionts_16s)

ss_palette <- wes_palette("Darjeeling2", n = 13, type = "continuous")

sulfur_symbionts_bar_16s <- ggplot(sulfur_symbionts_df_16s, aes(x = Sample, y = Abundance, fill = V6)) +
  geom_bar(stat = "identity", na.rm = TRUE) +
  theme_bw() +
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  scale_color_manual(values = ss_palette) +
  scale_fill_manual(values = ss_palette, name = "Genus") +
  facet_nested(~factor(Habitat, levels = c("Campbell Cove","West Park","Mason's Marina")) + Location, scales = "free", drop = TRUE) 



sulfur_symbionts_bar_16s






###################### Aldex2 analysis on phyloseq object ASV level ####
aldex2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phylo_16s)),
                                     phyloseq::sample_data(phylo_16s)$Habitat,
                                     mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_asv_table <- data.frame(phyloseq::otu_table(phylo_16s))
aldex2_asv_table <- rownames_to_column(aldex2_asv_table, var = "OTU")
write.csv(aldex2_asv_table, "16s/Tables:Results/aldex2_asv_table.csv")

# Write aldex2 results to file as csv
write.csv(aldex2, "16s/Tables:Results/aldex2_results.csv")

# Import aldex2 results and rename the X variable by OTU
aldex2_asv_result <- read.csv("16s/Tables:Results/aldex2_results.csv")
colnames(aldex2_asv_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

#Clean up presentation
taxa_info <- data.frame(tax_table(phylo_16s))
taxa_info <- taxa_info %>%
  rownames_to_column(var = "OTU")
sample_tab <- data.frame(sample_data(phylo_16s))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_result <- aldex2_asv_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_result <- left_join(sig_aldex2_result, taxa_info)
write.csv(sig_aldex2_result, "16s/Tables:Results/sig_aldex2_result.csv") 

# Create clr objects by using OTU ids and original otu table.
# All significant OTUs
sig_aldex2_result_count <- left_join(sig_aldex2_result, aldex2_asv_table)
write.csv(sig_aldex2_result, "16s/Tables:Results/sig_aldex2_result.csv")

clr_asv <- sig_aldex2_result_count[, -(2:10)]
rownames(clr_asv) <- clr_asv$OTU
clr_asv <- clr_asv[, -1]

# Adjusting zeros on the matrix and applying log transformation
clr_asv_czm <- cmultRepl(t(clr_asv),  label=0, method="CZM")
clr_asv_czm_tv <- t(apply(clr_asv_czm, 1, function(x){log(x) - mean(log(x))}))
clr_asv_czm <- (apply(clr_asv, 1, function(x){log(x+1) - mean(log(x+1))}))

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_matrix4 <- colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red"))

# Combine the heatmap and the annotation
Z.Score.clr_asv <- scale(t(clr_asv_czm))
heatmap(Z.Score.clr_asv) # simple heatmap

ha_top = HeatmapAnnotation(
  Location = as.vector(sample_tab$Location),
  Habitat = as.vector(sample_tab$Habitat),
  col = list(
    Location = c("Bare Sediment" = "#674422", "Sea Grass" = "#004B00"),
    Habitat = c("Mason's Marina" = "#52392F", "West Park" = "#C29B6C",
                "Campbell Cove" = "#397A4C")
  ),
  annotation_legend_param = list(
    Location = list(
      title = "Location",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    ),
    Habitat = list(
      title = "Habitat",
      at = c("Mason's Marina", "West Park", "Campbell Cove"),
      labels = c ("MM", "WP", "CC")
    )
  ))

# Organize total abundance and taxa name
# Taxa name
Z.Score.clr_asv_name <- as.data.frame(Z.Score.clr_asv)
str(Z.Score.clr_asv)
Z.Score.clr_asv_name <- rownames_to_column(Z.Score.clr_asv_name, var = "OTU")
Z.Score.clr_asv_name <- left_join(Z.Score.clr_asv_name, taxa_info)
head(Z.Score.clr_asv_name)
#Z.Score.clr_asv_name <- Z.Score.clr_asv_name[, -(85:91)]

# Taxa abundance
aldex2_asv_table_total <- as.data.frame(aldex2_asv_table)
aldex2_asv_table_total$Total <- rowSums(aldex2_asv_table[, -1])
head(aldex2_asv_table_total)
aldex2_asv_table_total <- aldex2_asv_table_total[, -(2:84)]

Z.Score.clr_asv_count_total <- left_join(Z.Score.clr_asv_name, aldex2_asv_table_total)
head(Z.Score.clr_asv_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 2000, 4000, 6000, 8000),
                               c("white", "#c7eae5", "#80cdc1","#35978f","#01665e"))

ha_right_phy = rowAnnotation(
  Abundance = Z.Score.clr_asv_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_phy = Z.Score.clr_asv_count_total$V5

#row_labels_phy = Z.Score.clr_asv_count_total$OTU


# Plot heatmap at the phylum level
hm_asv <- Heatmap(Z.Score.clr_asv, name = "Z-score, CLR", col = col_matrix,
                       column_title = "16S rRNA microbiome", 
                       column_title_gp = gpar(fontface = "bold", fontsize = 16),
                       column_split = as.vector(as.vector(sample_tab$Location)),
                       border = TRUE,
                       top_annotation = ha_top,
                       right_annotation = ha_right_phy,
                       row_title = "ASV",
                       row_labels = row_labels_phy,
                       row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 6),
                       row_order = order(row_labels_phy),
                       rect_gp = gpar(col = "white", lwd = 1),
                       show_column_names = FALSE,
                       show_heatmap_legend = TRUE) 
  
hm_asv










###################### Alpha diversity measures #########################
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset

alpha_div_16S <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_16s, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_16s, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_16s, measures = "InvSimpson"),
  "Location" = phyloseq::sample_data(phylo_16s)$Location,
  "Habitat" = phyloseq::sample_data(phylo_16s)$Habitat)

head(alpha_div_16S)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S <- alpha_div_16S %>%
  rename(Simpson = InvSimpson)
head(alpha_div_16S)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S <- alpha_div_16S[, c( 5, 4, 1, 2, 3)]
head(alpha_div_16S)

#Summarize alpha diversity measures by location
summary_alpha_16S_location <- alpha_div_16S %>%
  group_by(Location) %>%
  summarise(count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_Simpson = mean(Simpson),
    sd_Simpson = sd(Simpson))

write_csv(summary_alpha_16S_location, "16s/Tables:Results/alpha_div_results_16s_location.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_16S_habitat <- alpha_div_16S %>%
  group_by(Habitat) %>%
  summarise(count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_Simpson = mean(Simpson),
    sd_Simpson = sd(Simpson))

write_csv(summary_alpha_16S_habitat, "16s/Tables:Results/alpha_div_results_16s_habitat.csv")

# KW analysis alpha_div_16S on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by location, column 2
kw_alpha_16S_location <- as.data.frame(sapply(3:5, function(x) kruskal.test(alpha_div_16S[,x], alpha_div_16S[,2])))

# Rename columns with the proper variable names
kw_alpha_16S_location <- kw_alpha_16S_location %>%
  rename(Observed = V1,
         Shannon = V2,
         Simpson = V3)

kw_alpha_16S_location <- t(kw_alpha_16S_location) # transpose
kw_alpha_16S_location <- as_tibble(kw_alpha_16S_location, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_location <- kw_alpha_16S_location[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_location) # checking object class


# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_location, "16s/Tables:Results/kw_ad_16s_location.csv")

# KW analysis alpha_div_16S on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by habitat, column 1
kw_alpha_16S_habitat <- as.data.frame(sapply(3:5, function(x) kruskal.test(alpha_div_16S[,x], alpha_div_16S[,1])))

# Rename columns with the proper variable names
kw_alpha_16S_habitat <- kw_alpha_16S_habitat %>%
  rename(Observed = V1,
         Shannon = V2,
         Simpson = V3)

kw_alpha_16S_habitat <- t(kw_alpha_16S_habitat) # transpose
kw_alpha_16S_habitat <- as_tibble(kw_alpha_16S_habitat, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_habitat <- kw_alpha_16S_habitat[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_habitat) # checking object class


# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_habitat, "16s/Tables:Results/kw_ad_16s_habitat.csv")

# Plot alpha diversity measures

location_comparisons <- list( c("West Park", "Mason's Marina"), c("West Park", "Campbell Cove"), c("Campbell Cove", "Mason's Marina"))

alpha_color_habitat <- c("#2166AC", "#B2182B", "#F4A582")

alpha_div_16S <- transform(alpha_div_16S, Habitat=factor(Habitat,levels=c("Campbell Cove","West Park","Mason's Marina")))

ad_16s_location <- alpha_div_16S %>%
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
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("CC", "WP", "MM"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")

ad_16s_location

alpha_color_location <- c("saddlebrown","seagreen3")

habitat_comparisons <- list(c("Bare Sediment", "Sea Grass"))

ad_16s_habitat <- alpha_div_16S %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Location, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label =..p.signif..), size = 4, hide.ns = FALSE) +
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

ad_16s_habitat









###################### Aldex2 analysis at Genus level for habitat comparison ########


phylo_16s_genus <- tax_glom(phylo_16s, taxrank = "V6")

aldex2_gen <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phylo_16s_genus)),
                                phyloseq::sample_data(phylo_16s_genus)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
                                
gen_otu_table <- data.frame(phyloseq::otu_table(phylo_16s_genus))
gen_otu_table <- rownames_to_column(gen_otu_table, var = "OTU")
write.csv(gen_otu_table, "16s/Tables:Results//gen_otu_table.csv")

write.csv(aldex2_gen, "16s/Tables:Results//aldex2_gen.csv")

aldex2_gen_result <- read.csv("16s/Tables:Results//aldex2_gen.csv")
colnames(aldex2_gen_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex_taxa_info <- data.frame(tax_table(phylo_16s))
aldex_taxa_info <- aldex_taxa_info %>%
  rownames_to_column(var = "OTU")
sample_tab <- data.frame(sample_data(phylo_16s))

write.csv(aldex_taxa_info, "16s/Tables:Results//taxa_info.csv")

aldex_taxa_info <- read.csv("16s/Tables:Results//taxa_info.csv")

# filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_gen_result <- aldex2_gen_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_gen_result <- left_join(sig_aldex2_gen_result, aldex_taxa_info)
write.csv(sig_aldex2_gen_result, "16s/Tables:Results//sig_aldex2_gen_result.csv")

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 30
sig_aldex2_gen_result_count <- left_join(sig_aldex2_gen_result, gen_otu_table)
write.csv(sig_aldex2_gen_result_count, "16s/Tables:Results//sig_aldex2_gen_result_count.csv")

clr_gen <- sig_aldex2_gen_result_count[, -(2:11)] ###
rownames(clr_gen) <- clr_gen$OTU
clr_gen <- clr_gen[, -1] ##

# Adjusting zeros on the matrix and applying log transformation
shsk_gen_czm <- cmultRepl(t(clr_gen),  label=0, method="CZM")
shsk_gen_czm_tv <- t(apply(shsk_gen_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_gen_czm <- (apply(clr_gen, 1, function(x){log(x+1) - mean(log(x+1))}))


# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)


# Combine the heatmap and the annotation
Z.Score.gen <- scale(t(shsk_gen_czm))

# Define colors for each level of qualitative variables, i.e. soil type and habitats
# Create the heatmap annotation
ha_top_gen = HeatmapAnnotation(
  Location = as.vector(sample_tab$Location),
  Habitat = as.vector(sample_tab$Habitat),
  col = list(
    Location = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "seagreen3"),
    Habitat = c("Mason's Marina" = "#2166AC", "West Park" = "#B2182B",
                "Campbell Cove" = "#F4A582")
  ),
  annotation_legend_param = list(
    Location = list(
      title = "Location",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    ),
    Habitat = list(
      title = "Habitat",
      at = c("Mason's Marina", "West Park", "Campbell Cove"),
      labels = c ("MM", "WP", "CC")
    )
  ))

# Organize total abundance and taxa name
Z.Score.gen_name <- as.data.frame(Z.Score.gen)
str(Z.Score.gen_name)
Z.Score.gen_name <- rownames_to_column(Z.Score.gen_name, var = "OTU")
Z.Score.gen_name <- left_join(Z.Score.gen_name, aldex_taxa_info)
head(Z.Score.gen_name)
#Z.Score.gen_name <- Z.Score.gen_name[, -(84:90)] ###

gen_otu_table_total <- as.data.frame(gen_otu_table)
gen_otu_table_total$Total <- rowSums(gen_otu_table[, -1])
head(gen_otu_table_total)
gen_otu_table_total <- gen_otu_table_total[, -(2:84)] ##

Z.Score.gen_count_total <- left_join(Z.Score.gen_name, gen_otu_table_total)
head(Z.Score.gen_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 15000, 25000, 50000),
                               c("#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))

ha_right_gen = rowAnnotation(
  Abundance = Z.Score.gen_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_gen = Z.Score.gen_count_total$V6

#row_labels_gen = Z.Score.gen_count_total$V6

# Plot heatmap at the phylum level
hm_gen <- Heatmap(Z.Score.gen, name = "Z-score, CLR", col = col_matrix,
                         column_title = "16S rRNA microbiome", 
                         column_title_gp = gpar(fontface = "bold", fontsize = 14),
                         column_split = as.vector(as.vector(sample_tab$Habitat)),
                         #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                         border = TRUE,
                         top_annotation = ha_top_gen,
                         right_annotation = ha_right_gen,
                         row_title = "Genus",
                         row_labels = row_labels_gen,
                         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         #row_order = order(row_labels_gen),
                         rect_gp = gpar(col = "white", lwd = 1),
                         show_column_names = FALSE,
                         show_heatmap_legend = TRUE)
hm_gen

###################### Aldex2 analysis at Genus level for location comparison ######

# Agglomerate taxa at genus level 
phylo_16s_genus <- tax_glom(phylo_16s, taxrank = "V6")

# Run ALEDx2
aldex2_gen_location <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phylo_16s_genus)),
                                     phyloseq::sample_data(phylo_16s_genus)$Location,
                                     mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
# Create data frame and format
gen_otu_table_location <- data.frame(phyloseq::otu_table(phylo_16s_genus))
gen_otu_table_location <- rownames_to_column(gen_otu_table_location, var = "OTU")

write.csv(aldex2_gen_location, "16s/Tables:Results//aldex2_gen_l.csv")

# Import results CSV and format
aldex2_gen_result_location <- read.csv("16s/Tables:Results//aldex2_gen_l.csv")
colnames(aldex2_gen_result_location) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info <- data.frame(tax_table(phylo_16s))
aldex_taxa_info <- aldex_taxa_info %>%
  rownames_to_column(var = "OTU")

# Export CSV and edit "Unassigned" names at V6 level 
write.csv(aldex_taxa_info, "16s/Tables:Results//taxa_info.csv")

# Import edited taxonomy CSV
aldex_taxa_info <- read.csv("16s/Tables:Results//taxa_info.csv")


# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_gen_result_location <- aldex2_gen_result_location %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_gen_result_location <- left_join(sig_aldex2_gen_result_location, aldex_taxa_info)
write.csv(sig_aldex2_gen_result_location, "16s/Tables:Results//sig_aldex2_gen_result_location.csv")

# Made changes to taxonomy names if needed 
sig_aldex2_gen_result_location <- read.csv("16s/Tables:Results//sig_aldex2_gen_result_location.csv")

# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_gen_result_count_location <- left_join(sig_aldex2_gen_result_location, gen_otu_table_location)
sig_aldex2_gen_result_count_location <- sig_aldex2_gen_result_count_location[, -1]
write.csv(sig_aldex2_gen_result_count_location, "16s/Tables:Results//sig_aldex2_gen_result_count_location.csv")

# Dropped taxonomic rank and formatted columns 
clr_gen_location <- sig_aldex2_gen_result_count_location[, -(2:11)] 
rownames(clr_gen_location) <- clr_gen_location$OTU
clr_gen_location <- clr_gen_location[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_gen_czm_location <- cmultRepl(t(clr_gen_location),  label=0, method="CZM")
shsk_gen_czm_tv_location <- t(apply(shsk_gen_czm_location, 1, function(x){log(x) - mean(log(x))}))
shsk_gen_czm_location <- (apply(clr_gen_location, 1, function(x){log(x+1) - mean(log(x+1))}))


# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)


# Combine the heatmap and the annotation
Z.Score.gen_location <- scale(t(shsk_gen_czm_location))

# Define colors for each level of qualitative variables, i.e. locations and habitats
# Create the heatmap annotation
ha_top_gen = HeatmapAnnotation(
  Location = as.vector(sample_tab$Location),
  Habitat = as.vector(sample_tab$Habitat),
  col = list(
    Location = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "seagreen3"),
    Habitat = c("Mason's Marina" = "#2166AC", "West Park" = "#B2182B",
                "Campbell Cove" = "#F4A582")
  ),
  annotation_legend_param = list(
    Location = list(
      title = "Location",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    ),
    Habitat = list(
      title = "Habitat",
      at = c("Mason's Marina", "West Park", "Campbell Cove"),
      labels = c ("MM", "WP", "CC")
    )
  ))


# Organize total abundance and taxa name
# Taxa name
Z.Score.gen_name_location <- as.data.frame(Z.Score.gen_location)
str(Z.Score.gen_name_location)
Z.Score.gen_name_location <- rownames_to_column(Z.Score.gen_name_location, var = "OTU")
Z.Score.gen_name_location <- left_join(Z.Score.gen_name_location, aldex_taxa_info)
head(Z.Score.gen_name_location)
#Z.Score.gen_name <- Z.Score.gen_name[, -(84:90)] #

# Taxa abundance 
gen_otu_table_total_location <- as.data.frame(gen_otu_table_location)
gen_otu_table_total_location$Total <- rowSums(gen_otu_table_location[, -1])
head(gen_otu_table_total_location)
gen_otu_table_total_location <- gen_otu_table_total_location[, -(2:84)] 

Z.Score.gen_count_total_location <- left_join(Z.Score.gen_name_location, gen_otu_table_total_location)
head(Z.Score.gen_count_total_location)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 15000, 25000, 50000),
                               c("#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))

hb_right_gen = rowAnnotation(
  Abundance = Z.Score.gen_count_total_location$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_gen_location = Z.Score.gen_count_total_location$V6


# Plot heatmap at the phylum level
hm_gen_location <- Heatmap(Z.Score.gen_location, name = "Z-score, CLR", col = col_matrix,
                           column_title  = "16S rRNA microbiome", 
                           column_title_gp = gpar(fontface = "bold", fontsize = 14),
                           column_split = as.vector(as.vector(sample_tab$Location)),
                           #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                           border = TRUE,
                           top_annotation = ha_top_gen,
                           right_annotation = hb_right_gen,
                           row_title = "Genus",
                           row_labels = row_labels_gen_location,
                           row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                           row_names_gp = gpar(fontsize = 6),
                           column_names_gp = gpar(fontsize = 6),
                           #row_order = order(row_labels_gen_location),
                           rect_gp = gpar(col = "white", lwd = 1),
                           show_column_names = FALSE,
                           show_heatmap_legend = TRUE)


hm_gen_location

