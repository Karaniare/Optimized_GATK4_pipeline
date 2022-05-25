#!/bin/bash
#SBATCH -J Fastq_bam_qc
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g

##### Loading modules that are necessary for the analysis
module load bwa/0.7.15 tabix/0.2.6 vcftools/0.1.16 datamash/1.3 bcftools/1.9 plink/1.90 gatk/4.2.2.0 sratoolkit/2.8.2-1 R/4.1.0 samtools trimmomatic/0.36 
module load gsl/2.7.1 zlib/1.2.11 RAiSD/2.8

###### Specifying analysis directories

fastq_dir="/users/kniare/data/shared/scripts"
ref_dir="/users/kniare/data/kniare/Pf6k_work/Analysis/Pipeline"
bam_dir="/users/kniare/data/kniare/Pf6k_work/Tools"
vcf_dir="/users/kniare/data/shared/vcfs/pf3k_release6"
unpaired_dir="/users/kniare/data/shared/scripts/Unpaired_ir"
gvcf_dir=""
stat_dir=""

#####
cd $bam_dir
for i in $(cat ID_list.tsv)
  do
    cd $fastq_dir
    TrimmomaticPE "$i"_R1_001.fastq.gz "$i"_R2_001.fastq.gz "$i"_R1_paired.fq.gz "$i"_R1_unpaired.fq.gz "$i"_R2_paired.fq.gz "$i"_R2_unpaired.fq.gz  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:3  SLIDINGWINDOW:5:20 -threads 12
    mv *unpaired* $unpaired_dir
    cd $bam_dir
    bwa mem  -t 10 -M -R "@RG\tID:"$i"\tLB:"$i"\tPL:illumina\tSM:"$i"\tPU:"$i"" $ref_dir/Pf3D7_human.fa $fastq_dir/"$i"_R1_paired.fq.gz $fastq_dir/"$i"_R2_paired.fq.gz > "$i".sam
    gatk --java-options "-Xmx80g -Xms80g" SamFormatConverter -R $ref/Pf3D7_human.fa -I "$i".sam -O "$i".bam
    gatk --java-options "-Xmx80g -Xms80g" CleanSam $ref/Pf3D7_human.fa -I "$i".bam -O "$i".clean.bam
    gatk --java-options "-Xmx80g -Xms80g" SortSam -R $ref_dir/Pf3D7_human.fa -I "$i".clean.bam -O "$i".sorted.bam -SO coordinate --CREATE_INDEX true
    gatk --java-options "-Xmx80g -Xms80g" MarkDuplicatesSpark -R $ref_dir/Pf3D7_human.fasta -I "$i".sorted.bam -O "$i".sorted.dup.bam
    samtools view -b -h "$i".sorted.dup.bam -T $ref_dir/Pf3D7.fasta -L $ref_dir/Pf3D7_core.bed > "$i".sorted.dup.pf.bam
    rm "$i".sam "$i".bam "$i".clean.bam "$i".sorted.bam
  done 
###### Computing the distribution of the read depth using GATK
cd $bam_dir
ls *.sorted.dup.pf.bam > bam.list
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
    do
       gatk --java-options "-Xmx80g -Xms80g" DepthOfCoverage \
       -R $ref_dir/Pf3D7.fasta \
       -O $stat_dir/chr"$i" \
       --omit-locus-table true \
       -I bams.list
       awk -F "\t" -v OFS="\t" '{print $0, $NF="chr'$i'"}' $stat_dir/"$i".sample_summary > $stat_dir/"$i".sample2_summary
       rm $stat_dir/"$i".sample_summary
   done
cd $stat_dir
cat *.sample2_summary | awk '!/sample_id/ {print $0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15' > ReadCoverage_final.tsv

######### Insert size calculation 
cd $bam_dir
for i in $(cat ID_list.tsv)
     do
         gatk CollectInsertSizeMetrics -I  "$i".sorted.dup.pf.bam  -O $stat_dir/"$i"_insert.txt -H $stat_dir/"$i"_histo.pdf  -M 0.05
         awk 'FNR>=8 && FNR<=8 {print $1,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$NF="'$i'"}' $stat_dir/"$i"_insert.txt > $stat_dir/"$i"_insert2.txt
         rm $stat_dir/"$i"_insert.txt 
     done
cd $stat_dir
cat *_insert2.txt > $prov/InsertSize_Final.txt
rm *_insert2.txt

#### Generating some additional bam statistics
cd $bam_dir

for i in $(cat ID_list.tsv)
  do
    samtools stats "$i".sorted.dup.bam | grep ^SN | cut -f 2- | awk -F"\t" '{print$2}' > $stat_dir/"$i"_bamstat.tsv
      datamash transpose<$stat_dir/"$i"_bamstat.tsv | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "'$i'"; print }' > $stat_dir/"$i"_bamstat_final.tsv
           rm $stat_dir/"$i"_bamstat.tsv
done
cd $stat_dir
cat *_bamstat_final.tsv > Bam_stats_Final.tsv

rm *_bamstat_final.tsv



