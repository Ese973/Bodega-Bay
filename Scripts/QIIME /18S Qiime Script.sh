# Step 1 Importing Data

# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/data-clean
OUTPUT=/home/gas51591/boba/updated-database/analysis/data-sequences.qza

# Script
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $INPUT \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path $OUTPUT


