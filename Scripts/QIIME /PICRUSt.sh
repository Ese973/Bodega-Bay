# Step 1 export reference sequences and feature table from QIIME2 artifacts
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2-amplicon
REFSEQ=/home/gas51591/boba/16s/analysis/data-sequences-ca-trimmed-dada2.qza
REFTABLE=/home/gas51591/boba/16s/analysis/16s-dada2-table-ca.qza
SEQOUT=/home/gas51591/boba/16s/picrust2/output/refseq/
TABLEOUT=/home/gas51591/boba/16s/picrust2/output/biomtable/

# Script
qiime tools export \
  --input-path ${REFSEQ} \
  --output-path ${SEQOUT}

qiime tools export \
  --input-path ${REFTABLE} \
  --output-path ${TABLEOUT}



# Step 2
# Activate Picrust2 and add Path Variables
module load PICRUSt2/2.6.2-foss-2023a
SEQ=/home/gas51591/boba/16s/picrust2/output/refseq/dna-sequences.fasta
BIOM=/home/gas51591/boba/16s/picrust2/output/biomtable/feature-table.biom
OUT=/home/gas51591/boba/16s/picrust2/final_outputs2/

# Script
picrust2_pipeline.py \
  -s ${SEQ} \
  -i ${BIOM} \
  -o ${OUT} -p 12 \
  --stratified \
  --in_traits EC,KO



# Step 3 add desctiptions to each EC number
# Activate Picrust2
module load PICRUSt2/2.6.2-foss-2023a

# Script
add_descriptions.py -i /home/gas51591/boba/16s/picrust2/final_outputs2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o /home/gas51591/boba/16s/picrust2/final_outputs2/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i /home/gas51591/boba/16s/picrust2/final_outputs2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o /home/gas51591/boba/16s/picrust2/final_outputs2/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i /home/gas51591/boba/16s/picrust2/final_outputs2/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o /home/gas51591/boba/16s/picrust2/final_outputs2/pathways_out/path_abun_unstrat_descrip.tsv.gz
