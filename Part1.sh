#!/bin/bash
#SBATCH -J Fastq_bam_qc
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g

##### Loading modules that are necessary for the analysis
module load bwa/0.7.15 tabix/0.2.6 vcftools/0.1.16 datamash/1.3 bcftools/1.9 plink/1.90 gatk/4.2.2.0 sratoolkit/2.8.2-1 R/4.1.0 samtools trimmomatic/0.36 
module load gsl/2.7.1 zlib/1.2.11 RAiSD/2.8

###### Specifying analysis directories

while getopts i:f:r:b:u:s: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        f) fastq_dir=${OPTARG};;
        r) ref_dir=${OPTARG};;
        b) bam_dir=${OPTARG};;
        u) unpaired_dir=${OPTARG};;
        s) stat_dir=${OPTARG};;
    esac
done
echo "Argument definition:"
echo "-i: input file (sample ID list) in tab separated format"
echo "-f: fastq file directory"
echo "-r: reference genome and bed file directory"
echo "-b: bam file directory"
echo "-v : vcf file (or output) directory"
echo "-u: directory where unpaired fastq files are kept"
echo "-g: gvcf file directory"
echo "-s: directory to save qc output files"

##fastq_dir="/users/kniare/data/shared/scripts"
##ref_dir="/users/kniare/data/kniare/Pf6k_work/Analysis/Pipeline"
##bam_dir="/users/kniare/data/kniare/Pf6k_work/Tools"
##vcf_dir="/users/kniare/data/shared/vcfs/pf3k_release6"
##unpaired_dir="/users/kniare/data/shared/scripts/Unpaired_ir"
##gvcf_dir=""
##stat_dir=""

#####
cd $bam_dir
for i in $(cat $input)
  do
    cd $fastq_dir
    TrimmomaticPE "$i"_R1_001.fastq.gz "$i"_R2_001.fastq.gz "$i"_R1_paired.fq.gz "$i"_R1_unpaired.fq.gz "$i"_R2_paired.fq.gz "$i"_R2_unpaired.fq.gz  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:3  SLIDINGWINDOW:5:20 -threads 12
    mv *unpaired* $unpaired_dir
    cd $bam_dir
    bwa mem  -t 10 -M -R "@RG\tID:"$i"\tLB:"$i"\tPL:illumina\tSM:"$i"\tPU:"$i"" $ref_dir/Pf3D7_human.fa $fastq_dir/"$i"_R1_paired.fq.gz $fastq_dir/"$i"_R2_paired.fq.gz > "$i".sam
    gatk --java-options "-Xmx40g -Xms40g" SamFormatConverter -R $ref_dir/Pf3D7_human.fa -I "$i".sam -O "$i".bam
    gatk --java-options "-Xmx40g -Xms40g" CleanSam -R $ref_dir/Pf3D7_human.fa -I "$i".bam -O "$i".clean.bam
    gatk --java-options "-Xmx40g -Xms40g" SortSam -R $ref_dir/Pf3D7_human.fa -I "$i".clean.bam -O "$i".sorted.bam -SO coordinate --CREATE_INDEX true
    gatk --java-options "-Xmx40g -Xms40g" MarkDuplicatesSpark -R $ref_dir/Pf3D7_human.fa -I "$i".sorted.bam -O "$i".sorted.dup.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/Pf3D7.fasta -L $ref_dir/Pf3D7_core.bed > "$i".sorted.dup.pf.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/genome.fa -L $ref_dir/human.bed > "$i".sorted.dup.hs.bam
    samtools index -bc "$i".sorted.dup.pf.bam
    rm "$i".sam "$i".bam* "$i".clean.bam* "$i".sorted.bam*
  done 
###### Computing the distribution of the read depth using GATK
cd $bam_dir
ls *.sorted.dup.pf.bam > bams_pf.list
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
    do
       gatk --java-options "-Xmx80g -Xms80g" DepthOfCoverage \
       -R $ref_dir/Pf3D7.fasta \
       -O $stat_dir/chr"$i" \
       -L Pf3D7_"$i"_v3 \
       --omit-locus-table true \
       -I bams_pf.list
       awk -F "\t" -v OFS="\t" '{print $0, $NF="chr'$i'"}' $stat_dir/chr"$i".sample_summary > $stat_dir/chr"$i".sample2_summary
       rm $stat_dir/chr"$i".sample_summary
   done
cd $stat_dir
cat *.sample2_summary | awk '!/sample_id/ {print $0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15' > ReadCoverage_final.tsv

######### Insert size calculation 
cd $bam_dir
for i in $(cat $input)
     do
         gatk CollectInsertSizeMetrics -I  "$i".sorted.dup.pf.bam  -O $stat_dir/"$i"_insert.txt -H $stat_dir/"$i"_histo.pdf  -M 0.05
         awk 'FNR>=8 && FNR<=8 {print $1,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$NF="'$i'"}' $stat_dir/"$i"_insert.txt > $stat_dir/"$i"_insert2.txt
         rm $stat_dir/"$i"_insert.txt 
     done
cd $stat_dir
cat *_insert2.txt > $prov/InsertSize_Final.txt
rm *_insert2.txt

#### Generating some additional bam statistics for P. falciparum only
cd $bam_dir

for i in $(cat $input)
  do
    samtools stats "$i".sorted.dup.pf.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat.tsv
      datamash transpose<$stat_dir/"$i"_bamstat.tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_final.tsv
           rm $stat_dir/"$i"_bamstat.tsv
done
cd $stat_dir
cat *_bamstat_final.tsv | sed '1irow_total_reads_pf,	filtered_reads_pf,	sequences_pf,	is_sorted_pf,	1st_fragments_pf,	last_fragments_pf,	reads_mapped_pf,	reads_mapped_and_paired_pf,	reads_unmapped_pf,	reads_properly_paired_pf,	reads_paired_pf,	reads_duplicated_pf,	reads_MQ0_pf,	reads_QC_failed_pf,	non_primary_alignments_pf,	supplementary_alignments_pf,	total_length_pf,	total_first_fragment_length_pf,	total_last_fragment_length_pf,	bases_mapped_pf,	bases_mapped_(cigar)_pf,	bases_trimmed_pf,	bases_duplicated_pf,	mismatches_pf,	error_rate_pf,	average_length_pf,	average_first_fragment_length_pf,	average_last_fragment_length_pf,	maximum_length_pf,	maximum_first_fragment_length_pf,	maximum_last_fragment_length_pf,	average_quality_pf,	insert_size_average_pf,	insert_size_standard_deviation_pf,	inward_oriented pairs_pf,	outward_oriented_pairs_pf,	pairs_with_other_orientation_pf,	pairs_on_different_chromosomes_pf,	percentage_of_properly_paired_reads_(%)_pf,
' > Bam_stats_pf_Final.tsv

rm *_bamstat_final.tsv

######### Generating some additional bam statistics for human genome only
cd $bam_dir

for i in $(cat $input)
  do
    samtools stats "$i".sorted.dup.hs.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat.tsv
      datamash transpose<$stat_dir/"$i"_bamstat.tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_final.tsv
           rm $stat_dir/"$i"_bamstat.tsv
done
cd $stat_dir
cat *_bamstat_final.tsv | sed '1irow_total_reads_hs,	filtered_reads_hs,	sequences_hs,	is_sorted_hs,	1st_fragments_hs,	last_fragments_hs,	reads_mapped_hs,	reads_mapped_and_paired_hs,	reads_unmapped_hs,	reads_properly_paired_hs,	reads_paired_hs,	reads_duplicated_hs,	reads_MQ0_hs,	reads_QC_failed_hs,	non_primary_alignments_hs,	supplementary_alignments_hs,	total_length_hs,	total_first_fragment_length_hs,	total_last_fragment_length_hs,	bases_mapped_hs,	bases_mapped_(cigar)_hs,	bases_trimmed_hs,	bases_duplicated_hs,	mismatches_hs,	error_rate_hs,	average_length_hs,	average_first_fragment_length_hs,	average_last_fragment_length_hs,	maximum_length_hs,	maximum_first_fragment_length_hs,	maximum_last_fragment_length_hs,	average_quality_hs,	insert_size_average_hs,	insert_size_standard_deviation_hs,	inward_oriented pairs_hs,	outward_oriented_pairs_hs,	pairs_with_other_orientation_hs,	pairs_on_different_chromosomes_hs,	percentage_of_properly_paired_reads_(%)_hs,
' > Bam_stats_hs_Final.tsv


rm *_bamstat_final.tsv

##### Generating human/parasite read ratios
pr -m -t -s\ Bam_stats_pf_Final.tsv Bam_stats_hs_Final.tsv | gawk '{print $7,$46}' | awk '!/reads_mapped/ {print $0}' | awk -F "\t" -v OFS="\t" '{print $1, $2, $2/$1}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > Ratios_hs_pf_reads.tsv
