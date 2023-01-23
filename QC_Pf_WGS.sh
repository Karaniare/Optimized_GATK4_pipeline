#!/bin/bash
#SBATCH -J Fastq_bam_qc
#SBATCH -t 24:00:00
#SBATCH -c 24
#SBATCH --mem=5g


###### Setting arguments

usage() { echo -e "Usage: $0 [-i <input list>] [-f <fastq directory>] [-b <bam directory>] [-u <unpaired directory>] [-s <stat directory>] [-k <kit name>]" 2>&1
 echo -e "-i : list of sample IDs in a text file "
 echo -e "-f: full path to the fastq files"
 echo -e "-u: directory to keep unpaired fastqs after trimming "
 echo -e "-s: directory to keep QC stats"
 echo -e "-k: name of the library prep kit (Ex.: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa)"

}

if [[ ! $@ =~ ^\-.+ ]]
then
  usage
fi


while getopts i:f:b:u:s:k: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        f) fastq_dir=${OPTARG};;
        b) bam_dir=${OPTARG};;
        u) unpaired_dir=${OPTARG};;
        s) stat_dir=${OPTARG};;
        k) lib_kit=${OPTARG};;
    esac
done

#ref_dir="/users/kniare/data/shared/scripts/Data_WGS_UNC"
ref_dir=/opt/data

#### Fastq trimming, reads mapping onto the reference genome, bam processing, subsampling human and Plasmodium falciparum specific bams
cd $bam_dir
for i in $(cat $input)
  do
    cd $fastq_dir
    TrimmomaticPE "$i"_1.fastq.gz "$i"_2.fastq.gz "$i"_R1_paired.fq.gz "$i"_R1_unpaired.fq.gz "$i"_R2_paired.fq.gz "$i"_R2_unpaired.fq.gz  ILLUMINACLIP:$lib_kit:2:30:10 LEADING:3 TRAILING:3 MINLEN:3  SLIDINGWINDOW:5:20 -threads 12
    mv *unpaired.fq.gz $unpaired_dir
    cd $bam_dir
    bwa mem -t 20 -M -R "@RG\tID:'$i'\tLB:'$i'\tPL:illumina\tSM:'$i'\tPU:'$i'" $ref_dir/Pf3D7_human.fa $fastq_dir/"$i"_R1_paired.fq.gz $fastq_dir/"$i"_R2_paired.fq.gz > "$i".sam
    gatk --java-options "-Xmx40g -Xms40g" SamFormatConverter -R $ref_dir/Pf3D7_human.fa -I "$i".sam -O "$i".bam
    gatk --java-options "-Xmx40g -Xms40g" CleanSam -R $ref_dir/Pf3D7_human.fa -I "$i".bam -O "$i".clean.bam
    gatk --java-options "-Xmx40g -Xms40g" SortSam -R $ref_dir/Pf3D7_human.fa -I "$i".clean.bam -O "$i".sorted.bam -SO coordinate --CREATE_INDEX true
    gatk --java-options "-Xmx40g -Xms40g" MarkDuplicatesSpark -R $ref_dir/Pf3D7_human.fa -I "$i".sorted.bam -O "$i".sorted.dup.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/Pf3D7.fasta -L $ref_dir/Pf3D7_core.bed > "$i".sorted.dup.pf.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/genome.fa -L $ref_dir/human.bed > "$i".sorted.dup.hs.bam
    samtools index -bc "$i".sorted.dup.pf.bam
    rm "$i".sam "$i".bam* "$i".clean.bam* "$i".sorted.bam*
 done 

######### Insert size calculation 
cd $bam_dir
for i in $(cat $input)
     do
         gatk CollectInsertSizeMetrics -I  "$i".sorted.dup.pf.bam  -O $stat_dir/"$i"_insert.txt -H $stat_dir/"$i"_histo.pdf  -M 0.05
         awk 'FNR>=8 && FNR<=8 {print $1,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$NF="'$i'"}' $stat_dir/"$i"_insert.txt > $stat_dir/"$i"_insert2.txt
         rm $stat_dir/"$i"_insert.txt 
     done
cd $stat_dir
cat *_insert2.txt > InsertSize_Final.tsv
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
cat *_bamstat_final.tsv | sed '1irow_total_reads_pf	filtered_reads_pf	sequences_pf	is_sorted_pf	1st_fragments_pf	last_fragments_pf	reads_mapped_pf	reads_mapped_and_paired_pf	reads_unmapped_pf	reads_properly_paired_pf	reads_paired_pf	reads_duplicated_pf	reads_MQ0_pf	reads_QC_failed_pf	non_primary_alignments_pf	total_length_pf	total_first_fragment_length_pf	total_last_fragment_length_pf	bases_mapped_pf	bases_mapped_(cigar)_pf	bases_trimmed_pf	bases_duplicated_pf	mismatches_pf	error_rate_pf	average_length_pf	average_first_fragment_length_pf	average_last_fragment_length_pf	maximum_length_pf	maximum_first_fragment_length_pf	maximum_last_fragment_length_pf	average_quality_pf	insert_size_average_pf	insert_size_standard_deviation_pf	inward_oriented pairs_pf	outward_oriented_pairs_pf	pairs_with_other_orientation_pf	pairs_on_different_chromosomes_pf	percentage_of_properly_paired_reads_(%)_pf	sample_id' > Bam_stats_pf_Final.tsv

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
cat *_bamstat_final.tsv | sed '1irow_total_reads_pf	filtered_reads_pf	sequences_pf	is_sorted_pf	1st_fragments_pf	last_fragments_pf	reads_mapped_pf	reads_mapped_and_paired_pf	reads_unmapped_pf	reads_properly_paired_pf	reads_paired_pf	reads_duplicated_pf	reads_MQ0_pf	reads_QC_failed_pf	non_primary_alignments_pf	total_length_pf	total_first_fragment_length_pf	total_last_fragment_length_pf	bases_mapped_pf	bases_mapped_(cigar)_pf	bases_trimmed_pf	bases_duplicated_pf	mismatches_pf	error_rate_pf	average_length_pf	average_first_fragment_length_pf	average_last_fragment_length_pf	maximum_length_pf	maximum_first_fragment_length_pf	maximum_last_fragment_length_pf	average_quality_pf	insert_size_average_pf	insert_size_standard_deviation_pf	inward_oriented pairs_pf	outward_oriented_pairs_pf	pairs_with_other_orientation_pf	pairs_on_different_chromosomes_pf	percentage_of_properly_paired_reads_(%)_pf	sample_id' > Bam_stats_hs_Final.tsv


rm *_bamstat_final.tsv

##### Generating human/parasite read ratios

paste Bam_stats_pf_Final.tsv Bam_stats_hs_Final.tsv | awk -v OFS="\t" '!/per/ {print$39,$7,$46}' | awk '{if($2==0) $2=1}1' |  awk '{if($3==0) $3=1}1' | awk -v OFS="\t" '{print$1,$2,$3,($3/$2)}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > Ratios_hs_pf_reads.tsv

############ Computing the distribution of the read depth using GATK

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
