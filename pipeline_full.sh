#!/bin/bash
#SBATCH -J Gatk4_testT_1
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g

##### Loading modules that are necessary for the analysis
module load bwa/0.7.15 tabix/0.2.6 vcftools/0.1.16 bcftools/1.9 plink/1.90 gatk/4.2.2.0 sratoolkit/2.8.2-1 R/4.1.0 samtools trimmomatic/0.36 
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


######## Variant Calling starts here
####### Running HaplotypeCaller to generate gVCFs
cd $gvcf
mkdir chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14
cd $bam_dir
for i in $(cat ID_list)
    do 
      for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
       do
          cd $gvcf_dir
          gatk --java-options "-Xmx80G" HaplotypeCaller -R $ref_dir/Pf3D7.fasta \
         -I "$bam_dir/$i".sorted.dup.pf.bam -ERC GVCF -ploidy 6 \
         --native-pair-hmm-threads 16 -O  chr"$j"/"$i".chr"$j".g.vcf --assembly-region-padding 100 \
         --max-num-haplotypes-in-population 128 --kmer-size 10 --kmer-size 25 \
         --min-dangling-branch-length 4 --heterozygosity 0.0029 --indel-heterozygosity 0.0017 \ 
         --min-assembly-region-size 100 -L $ref_dir/core_chr"$j".list \
         -mbq 5 -DF MappingQualityReadFilter  --base-quality-score-threshold 12
       done
done

####### Combining gVCFs into a genomic database 
##### Combining gvcfs 
cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  do
    gatk --java-options "-Xmx40G" GenomicsDBImport \
    --sample-name-map gvcf_chr"$i"_list.tsv \
    --genomicsdb-workspace-path chr"$i"_database \
    --intervals $ref_dir/core_chr"$i".list --batch-size 100 \
    --reader-threads 24 --genomicsdb-segment-size 8048576 \
    --genomicsdb-vcf-buffer-size 160384
done
#### Running joint genotyping by genomic segment (the core genome of each chromosome is split into smaller genomic regions) via multiple slrum jobs
cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  do
   cd $vcf_dir
   for j in $(cat Genomic_region_list.tsv)
     do
     echo -e "#SBATCH -J Genotype_chr"$i"_"$j"" > genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH -t 72:00:00" >> genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH -c 8" >> genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH --mem-per-cpu=10g" >> genotype_chr"$i"_"$j".sh
     echo -e "module load gatk/4.2.2.0" >> genotype_chr"$i"_"$j".sh
     echo -e "vcf_dir="/users/vcf_dir"" >> genotype_chr"$i"_"$j".sh
     echo -e "ref_dir="/users/ref_dir"" >> genotype_chr"$i"_"$j".sh
     echo -e "gatk --java-options "-Xmx80G" GenotypeGVCFs --genomicsdb-use-vcf-codec true  -R $ref_dir/Pf3D7.fasta -V gendb://$vcf_dir/chr"$i"_database --max-genotype-count 1024 -O $vcf_dir/chr"$i"_part"$j".vcf.gz --tmp-dir=$ref_dir -stand-call-conf 30 -L Pf3D7_13_v3:"$j"" >> genotype_chr"$i"_"$j".sh
     chmod +x genotype_chr"$i"_"$j".sh
     sbatch genotype_chr"$i"_"$j".sh
   done
done
##### Combining raw vcfs per chromosome  by genomic region order
cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  do
    cd $vcf_dir
    gatk --java-options "-Xmx80G" GatherVcfs -R $ref_dir/Pf3D7.fasta -I vcf_chr"$i"_list.tsv -O Chr"$i".raw.vcf.gz  --CREATE_INDEX true
    tabix -p vcf Chr"$i".raw.vcf.gz
    rm chr"$i"_part*.vcf.gz
done
###### This script  does raw vcf soft filtering via  Variant Quality Score Recalibration (VQSR) algorithm with our custom training dataset
    ### VariantRecalibrator for Indels
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
 do 
gatk --java-options "-Xmx80g -Xms80g" VariantRecalibrator \
-R Pf3D7.fasta \
-V Chr"$i".raw.vcf.gz \
--trust-all-polymorphic \
-an QD -an DP -an FS -an SOR -an MQ \
-mode INDEL \
--max-gaussians 4 \
-resource:Brown,known=true,training=true,truth=true,prior=15.0 Strains.vcf.gz \
-O  Chr"$i".raw.indel.recal \
--tranches-file  Chr"$i".raw.indel.tranches \
--rscript-file  Chr"$i".raw.indel.plots.R

    ### ApplyVQSR for Indels

gatk --java-options "-Xmx80g -Xms80g" \
   ApplyVQSR \
    -V Chr"$i".raw.vcf.gz \
    --recal-file Chr"$i".raw.indel.recal \
    --tranches-file Chr"$i".raw.indel.tranches \
    --create-output-variant-index true \
    --lod-score-cutoff 0.0 \
--exclude-filtered true \
    -mode INDEL \
    -O Chr"$i".raw.indel.recal.vcf.gz

 ### VariantRecalibrator for SNPs
gatk --java-options "-Xmx80g -Xms80g" VariantRecalibrator \
-R Pf3D7.fasta \
-V Chr"$i".raw.indel.recal.vcf.gz \
--trust-all-polymorphic \
-an QD -an DP -an FS -an SOR -an MQ \
-mode SNP \
--max-gaussians 4 \
-resource:Brown,known=true,training=true,truth=true,prior=15.0 Strains.vcf.gz \
-O Chr"$i".raw.snp.recal \
--tranches-file Chr"$i".raw.snp.tranches \
--rscript-file Chr"$i".raw.snp.plots.R

     #### AppyVQSR for SNPs
gatk --java-options "-Xmx80g -Xms80g" \
   ApplyVQSR \
    -V  Chr"$i".raw.indel.recal.vcf.gz \
    --recal-file  Chr"$i".raw.snp.recal \
    --tranches-file  Chr"$i".raw.snp.tranches \
    --create-output-variant-index true \
    --lod-score-cutoff 0.0 \
--exclude-filtered true \
    -mode SNP \
    -O  Chr"$i".raw.recal.vcf.gz
done

