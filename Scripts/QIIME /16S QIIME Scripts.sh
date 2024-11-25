# Step 1 Importing Data
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/16s/16s_manifest.txt
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences.qza

# Script
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $INPUT \
--output-path $OUTPUT \
--input-format PairedEndFastqManifestPhred33V2



# Step 2 Summarize Demultiplexed Sequences
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/16s/analysis/data-sequences.qza
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences.qzv

# Script
qiime demux summarize \
--i-data $INPUT \
--o-visualization $OUTPUT



# Step 3 Denoise with DADA2
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis/data-sequences.qza
OUTPUT=/home/gas51591/boba/16s/analysis/16s-denoised-rep-sequences.qza

# Script
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $INPUT \
--p-trim-left-f 0 \
--p-trunc-len-f 280 \
--p-trim-left-r 0 \
--p-trunc-len-r 135 \
--o-representative-sequences $OUTPUT \
--o-table /home/gas51591/boba/16s/analysis/16s-dada2-table.qza \
--o-denoising-stats /home/gas51591/boba/16s/analysis/16s-dada2-stats.qza



# Step 4 Exporting DADA2 Stats
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime metadata tabulate \
--m-input-file /home/gas51591/boba/16s/analysis/16s-dada2-stats.qza \
--o-visualization /home/gas51591/boba/16s/analysis/16s-dada2-stats.qzv

qiime feature-table summarize \
--i-table /home/gas51591/boba/16s/analysis/16s-dada2-table.qza \
--o-visualization /home/gas51591/boba/16s/analysis/16s-dada2-table.qzv \
--m-sample-metadata-file /home/gas51591/boba/16s/metadata/16s_metadata.txt

qiime feature-table tabulate-seqs \
--i-data /home/gas51591/boba/16s/analysis/16s-denoised-rep-sequences.qza \
--o-visualization /home/gas51591/boba/16s/analysis/16s-denoised-rep-sequences.qzv



# Step 5 Import Taxonomy Files
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/tax
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime tools import \
--type FeatureData[Sequence] \
--input-path /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.fasta \
--output-path /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.qza

qiime tools import \
--type FeatureData[Taxonomy] \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path /home/gas51591/boba/16s/tax/SILVA138-nr99_taxonomy_16S.txt \
--output-path /home/gas51591/boba/16s/tax/SILVA138-nr99_ref_taxonomy_16S.qza



# Step 6 BLAST Taxonomy 
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/tax
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime feature-classifier classify-consensus-blast \
--i-query /home/gas51591/boba/16s/analysis/16s-denoised-rep-sequences.qza \
--i-reference-reads /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.qza \
--i-reference-taxonomy /home/gas51591/boba/16s/tax/SILVA138-nr99_ref_taxonomy_16S.qza \
--p-maxaccepts 1 \
--p-perc-identity .90 \
--o-classification /home/gas51591/boba/16s/analysis/16s-rep-sequences-taxonomy.qza \
--o-search-results /home/gas51591/boba/16s/analysis/16s-rep-sequences-tax-search-results.qza



# Step 7 Taxonomic Visualization 
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis/16s-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/16s/analysis/16s-rep-sequences-taxonomy.qzv

# Script
qiime metadata tabulate \
--m-input-file $INPUT \
--o-visualization $OUTPUT



# Step 8 Create Tree for Phylogentic Diversity Analysis 
# Activate QIIME2-2024.2 in Miniconda3 Enviorment 
module load Miniconda3
source activate /home/gas51591/conda-envs/envs/QIIME2-2024.2
INPUT=/home/gas51591/boba/16s/analysis/18S-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/16s/analysis/

# Script
qiime phylogeny  align-to-tree-mafft-fasttree \
--i-sequences /home/gas51591/boba/16s/analysis/16s-denoised-rep-sequences.qza \
--o-alignment  /home/gas51591/boba/16s/analysis/aligned-16S-rep-seqs.qza \
--o-masked-alignment /home/gas51591/boba/16s/analysis/16s-masked-aligned-rep-seqs.qza \
--o-tree /home/gas51591/boba/16s/analysis/unrooted-16s-tree.qza \
--o-rooted-tree /home/gas51591/boba/16s/analysis/rooted-16s-tree.qza