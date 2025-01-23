###############################################
# 1. Trimming and filtering raw reads #
###############################################

#!/bin/sh
# -c 8
# --mem-per-cpu 6G

module load bbmap/39.01 seqkit
sample_list=$(cat sample_list.txt)
for i in $sample_list
  do
    cd "$i"
      bbduk.sh in1=${i}_1.fastq.gz in2=${i/_1/_2}_2.fastq.gz out1=/$(basename ${i/_1/_1})_trim_1.fq.gz out2=/$(basename ${i/_1/_2})_trim_2.fq.gz ref=adapters,phix k=23 hdist=1 tbo tpe;

      bbduk.sh in1=/$(basename ${i/_1/_1})_trim_1.fq.gz in2=/$(basename ${i/_1/_2})_trim_2.fq.gz out1=/$(basename ${i/_1/_1})_clean_1.fq.gz o
ut2=/$(basename ${i/_1/_2})_clean_2.fq.gz qtrim=rl trimq=12 minlength=51 maxns=2 maq=3;

      seqkit rmdup /$(basename ${i/_1/_1})_clean_1.fq.gz -s -o /$(basename ${i/_1/_1})_cleandp_1.fq.gz;
      seqkit rmdup /$(basename ${i/_1/_2})_clean_2.fq.gz -s -o /$(basename ${i/_1/_2})_cleandp_2.fq.gz;
    cd ..
  done

###############################################
# 2. Host removal #
###############################################

# --mem=600gb
#cpu=20 

mamba activate kraken2/2.1.2 pigz/2.6.0

i=$(cat names_list | sed -n "${SLURM_ARRAY_TASK_ID}p")
threads=20

path/conda/kraken2-2.1.2/bin/kraken2 \
-db /path/kraken2_RefSeqV205_Complete_500GB \
--paired /${i}_v_clean_1.fq.gz /${i}_v_clean_2.fq.gz \
--output /${i}_v_clean_12_krk2_500.kraken2 \
--report /${i}_v_clean_12_krk2_500.kraken2.report \
--use-names --threads $threads \

pigz -p $threads /${i}_v_clean_12_krk2_500.kraken2 /${i}_v_clean_12_krk2_500.kraken2.report \

/path/bin/extract_kraken_readsBB.py -k \
/path/${i}_v_clean_12_krk2_500.kraken2.gz \
-s1 /${i}_v_clean_1.fq.gz \
-s2 /${i}_v_clean_2.fq.gz \
-o /${i}_v_noEUK_1.fq \
-o2 /${i}_v_noEUK_2.fq \
--fastq-output  --exclude --taxid 2759 --include-children  -r //${i}_v_clean_12_krk2_500.kraken2.report.gz \

pigz -p $threads /${i}_v_noEUK_1.fq.gz /${i}_v_noEUK_2.fq.gz




###############################################
# 3.1. script for coassembly #
###############################################

#!/bin/sh
#-c 20
#mem=60gb

source ~/.bashrc
mamba activate spades

input="input_path/"
out="/output_path/"

threads=20

spades.py --meta --only-assembler -k 25,33,55,77,95,101,127 \
-1 $input/forward_all_virome.fq.gz \
-2 $input/reverse_all_virome.fq.gz \
-t $threads -m 60 \
-o $out/coassembly_contigs




###############################################
# 3.2. script for single sample assembly #
###############################################

#!/bin/sh
# -c 20
# --mem=60gb

source ~/.bashrc
mamba activate spades

input="input_path/"
out="/output_path/"

i=$(cat $input/sample_list | sed -n "${SLURM_ARRAY_TASK_ID}p")
threads=20

spades.py --meta --only-assembler -k 25,33,55,77,95,101,127 \
-1 $input/${i}_1.fq.gz \
-2 $input/${i}_2.fq.gz \
-t $threads -m 60 \
-o $out/${i}_contigs

###############################################
# 4. Virus identification #
###############################################

#!/bin/sh
# -c 5
# --mem=30gb

source ~/.bashrc
mamba activate genomad

input="input_path/"
out="/output_path/"
db="/genomad_db"

i=$(cat $out/sample_list | sed -n "${SLURM_ARRAY_TASK_ID}p")

genomad end-to-end --cleanup $input/${i}_contigs/contigs.fasta $out/${i}_geNomad_output $db

#!/bin/sh
# -c 5
# --mem=30gb

source ~/.bashrc
mamba activate genomad

input="input_path/"
out="/output_path/"
db="/genomad_db"

genomad end-to-end --cleanup $input/coassembly_contigs/contigs.fasta $out/coassembly_geNomad_output $db

## concatenate the fasta files from the predicted viruses

cat ./*_geNomad_output/contigs_summary/contigs_virus.fna > contigs_virus_all_ssANDco.fna

####################################################
# 5. dereplication using CheckV anicalc.py and aniclust.py #
####################################################

#!/bin/sh
# -c 3
# --mem=10gb

script="/path/" 
data="/path/"
wd="/path/"
ani="/path/" 

python $script/ani_clustering_script.py -wd $wd -i $data/contigs_virus_all_ssANDco.fna -t 3 -ani $ani/

#python script for running checkV anicalc.py and aniclust.py 
#! /usr/bin/env python

import argparse
import os
import subprocess
import pandas as pd

# Flag setup
parser = argparse.ArgumentParser(description='bl4stn - 1st pyth0n')
parser.add_argument('-i', '--input', help='Path to fasta file.', required=True)
parser.add_argument('-wd', help='Path to working_directory.', required=True)
parser.add_argument('-o', '--output', help='Name of the output directory.', required=True)
parser.add_argument('-t', '--threads', type=int, default=1)
parser.add_argument('-ani', help='Path to the anicalc.py and aniclust.py scripts.', required=True)
parser.add_argument('--min_ani', type=int, default=95, help='Minimum ANI value for clustering (default = 95)')
parser.add_argument('--min_tcov', type=int, default=85, help='Minimum coverage of the target sequence (default = 85)')
parser.add_argument('--min_qcov', type=int, default=0, help='Minimum coverage of the query sequence (default = 0)')

args = parser.parse_args()

# Create a directory for the output
set_folder = os.path.join(args.wd, args.output)
os.makedirs(set_folder, exist_ok=True)

# Make path to anicalc.py and aniclust.py
anicalc_path = os.path.join(args.ani, 'anicalc.py')
aniclust_path = os.path.join(args.ani, 'aniclust.py')

# Run makeblastdb fasta file
blastdb_output = os.path.join(set_folder, 'virus_contigs_blast_DB')
subprocess.run(['makeblastdb', '-in', args.input, '-dbtype', 'nucl', '-out', blastdb_output])

# Run blastn on the concatenated FASTA file
blastn_output = os.path.join(set_folder, 'viral_contigs_ready_for_cluster.tsv')
subprocess.run(['blastn', '-query', args.input, '-db', blastdb_output, '-outfmt', '6 std qlen slen', '-max_target_seqs', '10000', '-out', blastn_output, '-num_threads', str(args.threads)])

# Run anicalc.py on the blastn output
ani_output = os.path.join(set_folder, 'viral_contig_ANI_calculation.tsv')
subprocess.run(['python', anicalc_path, '-i', blastn_output, '-o', ani_output])

# Run aniclust.py on the anicalc.py output
cluster_output = os.path.join(set_folder, 'viral_contig_ANI_clusters.tsv')
subprocess.run(['python', aniclust_path, '--fna', args.input, '--ani', ani_output, '--out', cluster_output, '--min_ani', str(args.min_ani), '--min_tcov', str(args.min_tcov), '--min_qcov', str(args.min_qcov)])

# Add column names to the output file
df = pd.read_csv(cluster_output, sep='\t', header=None)
df.columns = ['Representative_Sequence', 'Sequences']
contig_ids_to_keep = set(df['Representative_Sequence'].tolist())
df.to_csv(cluster_output, sep='\t', index=False)

# Find the directory containing the viruses and proviruses fasta files for this sample
fasta_lines = []
with open(args.input, 'r') as f:
    fasta_lines = f.readlines()

# Filter the concatenated file
representative_fasta_file = os.path.join(set_folder, 'clust_virus_contigs.fna')
with open(representative_fasta_file, 'w') as f:
    for line in fasta_lines:
        if line.startswith('>'):
            contig_id = line.strip()[1:]
            if contig_id in contig_ids_to_keep:
                write_line = True
                f.write(line)
            else:
                write_line = False
        else:
            if write_line:
                f.write(line)



###################################################
# 5.2 rerun geNomad to get the phylogeny, length etc. #
###################################################

### filter out contigs larger than 3kb which have a geNomad fdr score of min. 0.8

seqtk seq <your_input.fasta> | awk 'BEGIN {RS=">"; ORS=""} length($0) > 3000 {print ">"$0}' > virus_3kb.fasta

#
vim filter_script.sh
#!/bin/bash

# Define the score threshold
threshold=0.8

# Create a list of contigs with score > threshold
awk -v threshold="$threshold" 'BEGIN {FS="\t"} $7 > threshold {print $1}' virus_3kb.tsv > contigs_to_keep.txt

# Filter the .fna file
while read contig; do
  awk -v contig="$contig" '/^>/{p = ($0 ~ contig)} p' virus_3kb.fasta
done < contigs_to_keep.txt > filtered_0.8_virus_3kb.fna

###################################################
# 6. profilling #
###################################################

vim profilling_cf.sh

#!/bin/sh                       
# -c 5
# --mem=5gb

source ~/.bashrc
mamba activate bowtie2_samtools_coverm

# Define variables
READS="/path/"
vOTUs="/path/"
OUT="/path/"
bowtie2-build -f $vOTUs/filtered_0.8_virus_3kb.fna $vOTUs/dbname
sample_list=$(cat $READS/sample_list)

for file in $sample_list
do
  echo "$file"
  output_sam="${OUT}/${file}.sam"
  output_bam="${OUT}/${file}.bam"
  sorted_bam="${OUT}/${file}_sorted.bam"
  output_coverm="${OUT}/${file}_CoverM.tsv"
  output_stats="${OUT}/${file}.stats"

  # Run Bowtie2
  bowtie2 -x $vOTUs/dbname -1 $READS/"$file"_v_noEUK_1.fq.gz -2 $READS/"$file"_v_noEUK_2.fq.gz -p $SLURM_CPUS_PER_TASK | samtools view -bS > "${output_sam}"

  # Convert SAM to BAM
  samtools view -S -b "${output_sam}" -o "${output_bam}" -@ $SLURM_CPUS_PER_TASK

  # Sort BAM files
  samtools sort "${output_bam}" -o "${sorted_bam}" -@ $SLURM_CPUS_PER_TASK

  # Run CoverM
  coverm contig -b "${sorted_bam}" -m mean trimmed_mean covered_bases covered_fraction variance length count reads_per_base rpkm tpm -o "${output_coverm}" --output-format sparse

  # Generate CoverM stats files
  echo -e 'vOTU\t'"${file}"'' > "${output_coverm_stats}"
  cut -f2,11 "${output_coverm}" | sed '1d' >> "${output_coverm_stats}" #workes wuhu!

  # Generate stats files
  echo -e 'vOTU\t'"${file}"'' > "${output_stats}"
  samtools index -b "${sorted_bam}"
  samtools idxstats "${sorted_bam}" | cut -f1,3 | sed '/*/d' >> "${output_stats}"
done


###################################################
# 7. checkV #
###################################################

vim checkv.sh

#!/bin/sh
# -c 16
# --mem=5gb

source ~/.bashrc
mamba activate checkv

in="/path/"
out="/path/"
db="/database/checkv-db-v1.5/"

checkv end_to_end $in/filtered_0.8_virus_3kb.fna $out/ -t 16 -d $db


###################################################
# 8. Host predictions using iPHoP #
###################################################

#!/bin/sh
# -c 8
# --mem=60gb                       # cores

source ~/.bashrc
mamba activate iphop_env

in="/path/"
out="/path/"
db="/iphop_db/Aug_2023_pub_rw"

iphop predict --fa_file $in/filtered_0.8_virus_3kb.fna --db_dir $db --out_dir $out/


###################################################
# 9.vOTU table and normalisation #
###################################################

vim vOTU_table.sh #worked!

#!/bin/bash
# -c 1
# --mem=1gb

# setup path to directories containing files
in="/path/"
out="/path/vOTU_table"

sample_list=$(cat $in/sample_list)

# Temporary directory to store intermediate files
temp_dir="$out/temp"
mkdir -p $temp_dir

# Use vOTUs.txt file for the first column with the header vOTU
cp $in/vOTUs.txt $temp_dir/votus.txt

# Initialize the final vOTU table with the vOTUs (contig names) as the first column
echo -e "vOTU,$(echo $sample_list | tr ' ' ',')" > $out/rpkm_0.1_abundance_table.csv

# Create a file with contig names to start building the table
cp $temp_dir/votus.txt $temp_dir/votu_table.csv

# Loop through each sample and extract the relevant RPKM values
for sample in $sample_list; do
    awk -F'\t' 'NR > 1 {if ($6 >= 0.1) print $11; else print 0}' $in/${sample}_CoverM.tsv > $temp_dir/${sample}_filtered_rpkm.txt
    paste -d',' $temp_dir/votu_table.csv $temp_dir/${sample}_filtered_rpkm.txt > $temp_dir/temp_table.csv
    mv $temp_dir/temp_table.csv $temp_dir/votu_table.csv
done

# Combine the vOTU names and the vOTU table into the final CSV file
awk 'BEGIN {FS=OFS=","} {print $0}' $temp_dir/votu_table.csv >> $out/rpkm_0.1_abundance_table.csv

# Cleanup temporary files
rm -rf $temp_dir

echo "vOTU table created at $out/rpkm_0.1_abundance_table.csv"


### downstream statistics and visualisation done using R
