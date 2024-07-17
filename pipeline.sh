##################
## step 1. QC - trim and filtering ##
##################

#!/bin/sh
#SBATCH -c 8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu 6G
#SBATCH --mail-type END
#SBATCH --mail-user caroline.winther-have@sund.ku.dk

module load bbmap/39.01 seqkit
sample_list=$(cat sample_list.txt)
for i in $sample_list
  do
    cd "$i"
      bbduk.sh in1=${i}_1.fastq.gz in2=${i/_1/_2}_2.fastq.gz out1=/home/zts270/projects/virome_speciale/01_trimmed/$(basename ${i/_1/_1})_trim_1.fq.gz out2=/home/zts270/projects/virome_speciale/01_trimmed/$(basename ${i/_1/_2})_trim_2.fq.gz ref=adapters,phix k=23 hdist=1 tbo tpe;

      bbduk.sh in1=/home/zts270/projects/virome_speciale/01_trimmed/$(basename ${i/_1/_1})_trim_1.fq.gz in2=/home/zts270/projects/virome_speciale/01_trimmed/$(basename ${i/_1/_2})_trim_2.fq.gz out1=/home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_1})_clean_1.fq.gz out2=/home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_2})_c
lean_2.fq.gz qtrim=rl trimq=12 minlength=51 maxns=2 maq=3;

      seqkit rmdup /home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_1})_clean_1.fq.gz -s -o /home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_1})_cleandp_1.fq.gz;
      seqkit rmdup /home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_2})_clean_2.fq.gz -s -o /home/zts270/projects/virome_speciale/02_filtered/$(basename ${i/_1/_2})_cleandp_2.fq.gz;
    cd ..
done


##################
## step 2. Host removal w. kraken2 ##
##################

#!/bin/zsh
#SBATCH --job-name=kraken_viromes.zsh    # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --mem=600gb
#SBATCH --ntasks=20 #cpu/threads jobs are allowed to use (doesnt multiply mem..)
#SBATCH --array=1-84%10
#SBATCH --mail-type END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user caroline.winther-have@sund.ku.dk     # Where to send mail
#SBATCH --time=1-00:00:00              # Time limit days:hrs:min:sec
#SBATCH --output=/home/zts270/projects/virome_speciale/03_noEU_kraken_virome/slurm_out/sample.%A_%a.log

module load kraken2/2.1.2 pigz/2.6.0

i=$(cat names_list | sed -n "${SLURM_ARRAY_TASK_ID}p")
threads=20

/projects/mjolnir1/apps/conda/kraken2-2.1.2/bin/kraken2 \
-db /projects/mjolnir1/data/databases/kraken2/kraken2_RefSeqV205_Complete_500GB \
--paired /home/zts270/projects/virome_speciale/02_filtered/samples_renamed/${i}_v_clean_1.fq.gz /home/zts270/projects/virome_speciale/02_filtered/samples_renamed/${i}_v_clean_2.fq.gz \
--output /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2 \
--report /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2.report \
--use-names --threads $threads \

pigz -p $threads /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2 /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2.report \

/projects/mjolnir1/apps/bin/extract_kraken_readsBB.py -k \
/home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2.gz \
-s1 /home/zts270/projects/virome_speciale/02_filtered/samples_renamed/${i}_v_clean_1.fq.gz \
-s2 /home/zts270/projects/virome_speciale/02_filtered/samples_renamed/${i}_v_clean_2.fq.gz \
-o /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_noEUK_1.fq \
-o2 /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_noEUK_2.fq \
--fastq-output  --exclude --taxid 2759 --include-children  -r /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_clean_12_krk2_500.kraken2.report.gz \

pigz -p $threads /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_noEUK_1.fq.gz /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/${i}_v_noEUK_2.fq.gz

##################
## step 3.1 SPAdes coassembly ##
##################

#script name: cf_spades_coassem_virome.sh

#!/bin/bash
#SBATCH --job-name=cf_spades_coassem_virome.sh    # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --mem=80gb
#SBATCH --ntasks=20 #cpu/threads jobs are allowed to use (doesnt multiply mem..)
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=caroline.winther-have@sund.ku.dk     # Where to send mail
#SBATCH --time=3-00:00:00              # Time limit hrs:min:sec
#SBATCH --output=/home/zts270/projects/virome_speciale/03_noEU_kraken_virome/slurm_coass/sample.coass.cf.%A_%a.log

module purge
module load spades/3.15.5

threads=20

spades.py --careful -k 21,25,33,55,77,95,99,127 \
-1 /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/cf_forward_all_virome.fq.gz \
-2 /home/zts270/projects/virome_speciale/03_noEU_kraken_virome/cf_reverse_all_virome.fq.gz \
-t $threads -m 80 \
-o /home/zts270/projects/virome_speciale/04_new_assemblies_viromes/cf_coassembly_virome_spades_contigs

##################
## step 3.2 SPAdes single sample assembly ##
##################

