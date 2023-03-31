#!/bin/bash
#SBATCH -J Fastq_bam_qc_"$run_name"
#SBATCH -t 24:00:00
#SBATCH -c 24
#SBATCH --mem-per-cpu=10g


###### Setting arguments

usage() { echo -e "Usage: $0 [-i <input list>] [-f <fastq directory>] [-b <bam directory>] [-u <unpaired directory>] [-s <stat directory>] [-k <kit name>]" 2>&1
 echo -e "-i : list of sample IDs in a text file "
 echo -e "-f: full path to the fastq files"
 echo -e "-u: directory to keep unpaired fastqs after trimming "
 echo -e "-s: directory to keep QC stats"
 echo -e "-k: name of the library prep kit (Ex.: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa)"
 echo -e "-r: run name (Ex: run1), useful if you run the pipeline on separate samplesets simultaneously"
 exit 1
}

if [[ ! $@ =~ ^\-.+ ]]
then
  usage
fi


while getopts i:f:b:u:s:k:r: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        f) fastq_dir=${OPTARG};;
        b) bam_dir=${OPTARG};;
        u) unpaired_dir=${OPTARG};;
        s) stat_dir=${OPTARG};;
        k) lib_kit=${OPTARG};;
        r) run_name=${OPTARG};;
    esac
done

#ref_dir="/users/kniare/data/shared/WGS_works_Karamoko/scripts/Data_WGS_UNC"
ref_dir=/opt/data

#### Fastq trimming, reads mapping onto the reference genome, bam processing, subsampling human and Plasmodium falciparum specific bams
cd $bam_dir
for i in $(cat $input)
  do
  cd $fastq_dir
    TrimmomaticPE "$i"_1.fastq.gz "$i"_2.fastq.gz "$i"_R1_paired.fq.gz "$i"_R1_unpaired.fq.gz "$i"_R2_paired.fq.gz "$i"_R2_unpaired.fq.gz  ILLUMINACLIP:$lib_kit:2:30:10 LEADING:3 TRAILING:3 MINLEN:36  SLIDINGWINDOW:5:20 -threads 12
    mv "$i"_R*_unpaired.fq.gz $unpaired_dir
    cd $bam_dir
    bwa mem -t 20 -M -R "@RG\tID:'$i'\tLB:'$i'\tPL:illumina\tSM:'$i'\tPU:'$i'" $ref_dir/Pf3D7_human.fa $fastq_dir/"$i"_R1_paired.fq.gz $fastq_dir/"$i"_R2_paired.fq.gz > "$i".sam
    gatk --java-options "-Xmx40g -Xms40g" SamFormatConverter -R $ref_dir/Pf3D7_human.fa -I "$i".sam -O "$i".bam
    gatk --java-options "-Xmx40g -Xms40g" CleanSam -R $ref_dir/Pf3D7_human.fa -I "$i".bam -O "$i".clean.bam
    gatk --java-options "-Xmx40g -Xms40g" SortSam -R $ref_dir/Pf3D7_human.fa -I "$i".clean.bam -O "$i".sorted.bam -SO coordinate --CREATE_INDEX true
    gatk --java-options "-Xmx40g -Xms40g" MarkDuplicates -R $ref_dir/Pf3D7_human.fa -I "$i".sorted.bam -O "$i".sorted.dup.bam -M "$i"_dup_metrics.tsv
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/Pf3D7.fasta -L $ref_dir/Pf3D7_core.bed > "$i".sorted.dup.pf.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/genome.fa -L $ref_dir/human.bed > "$i".sorted.dup.hs.bam
    samtools index -bc "$i".sorted.dup.pf.bam
    rm "$i".sam
 done 

######### Insert size calculation 
cd $bam_dir
for i in $(cat $input)
     do
         gatk CollectInsertSizeMetrics -I  "$i".sorted.dup.pf.bam  -O $stat_dir/"$i"_insert_"$run_name".txt -H $stat_dir/"$i"_"$run_name"_histo.pdf  -M 0.05
         awk 'FNR>=8 && FNR<=8 {print $1,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$NF="'$i'"}' $stat_dir/"$i"_insert_"$run_name".txt > $stat_dir/"$i"_insert2_"$run_name".txt
         rm $stat_dir/"$i"_insert_"$run_name".txt 
done
cd $stat_dir
cat *_insert2_"$run_name".txt > InsertSize_Final_"$run_name".tsv
rm *_insert2_"$run_name".txt

#### Generating some additional bam statistics for P. falciparum only
cd $bam_dir

for i in $(cat $input)
  do
######### Gettings bam statistics for total reads

    samtools stats "$i".sorted.dup.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat_total_"$run_name".tsv
     datamash transpose<$stat_dir/"$i"_bamstat_total_"$run_name".tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_total_final_"$run_name".tsv
     rm $stat_dir/"$i"_bamstat_total_"$run_name".tsv

####### Getting bam statistics for Plasmodium falciparum 3D7-specific reads
    
   samtools stats "$i".sorted.dup.pf.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat_"$run_name".tsv
    datamash transpose<$stat_dir/"$i"_bamstat_"$run_name".tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_final_"$run_name".tsv
    rm $stat_dir/"$i"_bamstat_"$run_name".tsv
done

cd $stat_dir
###### Adding header for total reads

cat *_bamstat_total_final_"$run_name".tsv | sed '1irow_total_reads_total filtered_reads_total sequences_total is_sorted_total 1st_fragments_total last_fragments_total reads_mapped_total reads_mapped_and_paired_total reads_unmapped_total reads_properly_paired_total reads_paired_total reads_duplicated_total reads_MQ0_total reads_QC_failed_total non_primary_alignments_total total_length_total total_first_fragment_length_total total_last_fragment_length_total bases_mapped_total bases_mapped_(cigar)_total bases_trimmed_total bases_duplicated_total mismatches_total error_rate_total average_length_total average_first_fragment_length_total average_last_fragment_length_total maximum_length_total maximum_first_fragment_length_total maximum_last_fragment_length_total average_quality_total insert_size_average_total insert_size_standard_deviation_total inward_oriented_pairs_total outward_oriented_pairs_total pairs_with_other_orientation_total pairs_on_different_chromosomes_total percentage_of_properly_paired_reads_(%)_total sample_id' > Bam_stats_total_Final_"$run_name".tsv


rm *_bamstat_final_"$run_name".tsv
########## Adding header for 3D7 reads

cat *_bamstat_final_"$run_name".tsv | sed '1irow_total_reads_pf filtered_reads_pf sequences_pf is_sorted_pf 1st_fragments_pf last_fragments_pf reads_mapped_pf reads_mapped_and_paired_pf reads_unmapped_pf reads_properly_paired_pf reads_paired_pf reads_duplicated_pf reads_MQ0_pf reads_QC_failed_pf non_primary_alignments_pf total_length_pf total_first_fragment_length_pf total_last_fragment_length_pf bases_mapped_pf bases_mapped_(cigar)_pf bases_trimmed_pf bases_duplicated_pf mismatches_pf error_rate_pf average_length_pf average_first_fragment_length_pf average_last_fragment_length_pf maximum_length_pf maximum_first_fragment_length_pf maximum_last_fragment_length_pf average_quality_pf insert_size_average_pf insert_size_standard_deviation_pf inward_oriented_pairs_pf outward_oriented_pairs_pf pairs_with_other_orientation_pf pairs_on_different_chromosomes_pf percentage_of_properly_paired_reads_(%)_pf sample_id' > Bam_stats_pf_Final_"$run_name".tsv



######### Generating some additional bam statistics for human genome only
cd $bam_dir

for i in $(cat $input)
  do
    samtools stats "$i".sorted.dup.hs.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat_hum_"$run_name".tsv
     datamash transpose<$stat_dir/"$i"_bamstat_hum_"$run_name".tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_hum_final_"$run_name".tsv
           rm $stat_dir/"$i"_bamstat_hum_"$run_name".tsv
done
cd $stat_dir
######## Adding header for human reads

cat *_bamstat_hum_final_"$run_name".tsv | sed '1irow_total_reads_hs filtered_reads_hs sequences_hs is_sorted_hs 1st_fragments_hs last_fragments_hs reads_mapped_hs reads_mapped_and_paired_hs reads_unmapped_hs reads_properly_paired_hs reads_paired_hs reads_duplicated_hs reads_MQ0_hs reads_QC_failed_hs non_primary_alignments_hs total_length_hs total_first_fragment_length_hs total_last_fragment_length_hs bases_mapped_hs bases_mapped_(cigar)_hs bases_trimmed_hs bases_duplicated_hs mismatches_hs error_rate_hs average_length_hs average_first_fragment_length_hs average_last_fragment_length_hs maximum_length_hs maximum_first_fragment_length_hs maximum_last_fragment_length_hs average_quality_hs insert_size_average_hs insert_size_standard_deviation_hs inward_oriented_pairs_hs outward_oriented_pairs_hs pairs_with_other_orientation_hs pairs_on_different_chromosomes_hs percentage_of_properly_paired_reads_(%)_hs sample_id' > Bam_stats_hs_Final_"$run_name".tsv


rm *_bamstat_final_"$run_name".tsv

##### Generating human/parasite read ratios

paste Bam_stats_pf_Final_"$run_name".tsv Bam_stats_hs_Final_"$run_name".tsv | awk -v OFS="\t" '!/per/ {print$39,$7,$46}' | awk '{if($2==0) $2=1}1' |  awk '{if($3==0) $3=1}1' | awk -v OFS="\t" '{print$1,$2,$3,($3/$2)}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > Ratios_hs_pf_reads_"$run_name".tsv

############ Computing the distribution of the read depth using GATK
cd $bam_dir
awk -v OFS="\t" '{print$1".sorted.dup.pf.bam"}' $input > bams_pf_"$run_name".list

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
    do
       gatk --java-options "-Xmx80g -Xms80g" DepthOfCoverage \
       -R $ref_dir/Pf3D7.fasta \
       -O $stat_dir/chr"$i"_"$run_name" \
       -L Pf3D7_"$i"_v3 \
       --omit-locus-table true \
       -I bams_pf_"$run_name".list
       awk -F "\t" -v OFS="\t" '{print $0, $NF="chr'$i'"}' $stat_dir/chr"$i"_"$run_name".sample_summary > $stat_dir/chr"$i"_"$run_name".sample2_summary
       rm $stat_dir/chr"$i"_"$run_name".sample_summary
   done
cd $stat_dir
cat *"$run_name".sample2_summary | awk '!/sample_id/ {print $0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15' > ReadCoverage_final_"$run_name".tsv
