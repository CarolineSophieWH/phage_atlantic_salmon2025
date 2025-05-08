#!/bin/sh
## Using deduplicated contigs from both single sample and coassembled SPAdes assemblies:
anvi-script-reformat-fasta contigs.fasta -o contigs-fixed_cf2.fa -l 0 --simplify-names --report_file
mv contigs-fixed_cf2.fa contigs_2cf.fa

# Make db
anvi-gen-contigs-database -f contigs_2cf.fa -o contigs2_cf.db -n CF_virome_mags2_db -T 8

# Run hmms
anvi-run-hmms -c contigs2_cf.db -T 8

# Display contigs stats
anvi-display-contigs-stats contigs2_cf.db -T 8

# Make COGs
anvi-run-ncbi-cogs -c contigs2_cf.db --cog-data-dir /projects/mjolnir1/data/databases/anvio-8.x/ncbi-cogs/20231006/ -T 8

# Profiling script:
#!/bin/sh
module load anvio/8.0

BAM='/path/'
CON='/path/'

cd $BAM
sample_list=$(cat sample_list)
for file in $sample_list
  do
  echo "$file"
  anvi-profile -i $BAM/"$file".bam -c $CON/contigs2_cf.db -o $CON/Profiles/"$file" -T 8 --cluster-contigs --profile-SCVs -M 500 --write-buffer-size 1000
done

# Merge profile databases
anvi-merge */PROFILE.db -o merged_samples -c /contigs2_cf.db -T 8

# Taxonomy 
anvi-run-scg-taxonomy -c contigs2_cf.db --scgs-taxonomy-data-dir /databases/anvio-8.x/scgs-taxonomy/20231006 -T 8

# Binning of clusters - METABAT2
module load metabat2/2.15
anvi-cluster-contigs -p /profile2_cf.db -c contigs2_cf.db -C METABAT2 --driver metabat2 --just-do-it --minCVSum 0 --minClsSize 0 -T 8

# Rename bins
anvi-rename-bins -c contigs2_cf.db -p profile2_cf.db --prefix CF --collection-to-read METABAT2 --collection-to-write MAGs_CF --report-file report_metabat2MAGs_CF.txt --call-MAGs --min-completion-for-MAG 70

# Show collections and bins
anvi-show-collections-and-bins -p profile2_cf.db

# Manual binning
anvi-interactive -p profile2_cf.db -c contigs2_cf.db
anvi-script-add-default-collection -c contigs2_cf.db -p profile2_cf.db -C DEFAULT -b DEFAULT
anvi-refine -c contigs2_cf.db -p profile2_cf.db -C DEFAULT -b DEFAULT
anvi-refine -c contigs2_cf.db -p profile2_cf.db -C DEFAULT -b CF_MAG_001
anvi-refine -c contigs2_cf.db -p profile2_cf.db -C DEFAULT -b CF_MAG_002

# Summarise
anvi-summarize -c contigs2_cf.db -p profile2_cf.db -C somename -o Anvio_summary_CF_no2

# KEGG setup and run
anvi-setup-kegg-data
anvi-run-kegg-kofams -c contigs2_cf.db -T 4

# Refine MAGs
anvi-refine -c contigs2_cf.db -p profile2_cf.db -C Jacob -b Bin_1 --server-only
anvi-rename-bins -c contigs2_cf.db -p profile2_cf.db --prefix CF --collection-to-read Jacob --collection-to-write Refined_MAGs --report-file report_Refined_MAGs_CF.txt --call-MAGs --min-completion-for-MAG 50
anvi-summarize -c contigs2_cf.db -p profile2_cf.db -C Refined_MAGs -o Anvio_summary

# Pangenome
mkdir PANGENOME
cp Anvio_summary/bin_by_bin/CF_MAG_00002/*fa PANGENOME
anvi-gen-contigs-database -f CF_MAG_00002-contigs.fa -o CF_MAG_00002-contigs.db
ls *db > external-genomes.txt
anvi-gen-genomes-storage -e external-genomes.txt -o MYCOPLASMA-GENOMES.db
anvi-pan-genome -g MYCOPLASMA-GENOMES.db --project-name "MYCOPLASMA_PAN" --output-dir MYCOPLASMA_PAN --num-threads 4 --minbit 0.5 --mcl-inflation 10
anvi-compute-genome-similarity --external-genomes external-genomes.txt --program pyANI --output-dir ANI --num-threads 4 --pan-db MYCOPLASMA_PAN/MYCOPLASMA_PAN.db
anvi-display-pan -p MYCOPLASMA_PAN/MYCOPLASMA_PAN.db -g MYCOPLASMA-GENOMES.db

# Metabolism estimation
anvi-estimate-metabolism -C Refined_MAGs -i internal_genomes_CF.tsv --matrix-format
anvi-matrix-to-newick kegg-metabolism-module_pathwise_completeness-MATRIX.txt
anvi-export-state -p profile_metabolism.db -s default -o default.json
anvi-import-state -s default.json -p metabolism_profile.db -n default
anvi-interactive -d kegg-metabolism-module_pathwise_completeness-MATRIX.txt -p MAG_metabolism_profile.db --dry-run --manual-mode
anvi-interactive --manual-mode -d kegg-metabolism-module_pathwise_completeness-MATRIX.txt -t kegg-metabolism-module_pathwise_completeness-MATRIX.txt.newick -p MAG_metabolism_profile.db --title "MAG CF Metabolism Heatmap"
anvi-import-state -s default.json -p MAG_metabolism_profile.db -n default

# Modules info
export ANVIO_MODULES_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/MODULES.db'))"`
echo -e "module\tclass\tcategory\tsubcategory\tname" > modules_info.txt
sqlite3 $ANVIO_MODULES_DB "select module, data_value from modules where data_name='CLASS'" | sed 's/; /|/g' | tr '|' '\t' >> module_class.txt
sqlite3 $ANVIO_MODULES_DB "select module,data_value from modules where data_name='NAME'" | tr '|' '\t' > module_names.txt
paste module_class.txt <(cut -f 2 module_names.txt ) >> modules_info.txt
rm module_names.txt module_class.txt
anvi-import-misc-data modules_info.txt -p MAG_metabolism_profile.db -t items
mv default.json MAG_metabolism.json
anvi-import-state -s MAG_metabolism.json -p MAG_metabolism_profile.db -n default
anvi-interactive --manual-mode -d kegg-metabolism-module_pathwise_completeness-MATRIX.txt -t kegg-metabolism-module_pathwise_completeness-MATRIX.txt.newick -p MAG_metabolism_profile.db --title "MAG CF Metabolism Heatmap"

# E.g. Estimate metabolism Mycoplasma
anvi-estimate-metabolism -c contigs2_cf.db -p profile2_cf.db -C Refined_MAGs -b CF_MAG_00002

# E.g. Anvio phylogeny for Aliivibrio
for file in *_fna
  do
    mv -- "$file" "${file//_fna/.fna}"
done

for i in `ls *.fna | awk 'BEGIN{FS=".fna"}{print $1}'`
do
  echo $i
  anvi-script-reformat-fasta $i.fna -o $i.fixed.fna --simplify-names
  anvi-gen-contigs-database -f $i.fixed.fna -o $i.db -T 4
  anvi-run-hmms -c $i.db -T 4
done

# Make fasta of concatenated hmms
anvi-get-sequences-for-hmm-hits --external-genomes external_genomes_aliivibrio_CF.tsv \
                                -o concatenated-proteins.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate

# Make tree
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa -o phylogenomic-tree.txt

# Visualize with anvi'o
anvi-interactive -p phylogenomic-profile.db \
                 -t phylogenomic-tree.txt \
                 --title "Phylogenomics of Aliivibrio" \
                 --manual
