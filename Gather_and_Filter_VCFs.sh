#!/bin/bash
#SBATCH -J Variant_calling
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g
#
###### Creating arguments
usage() {
   echo "Usage $0 or Gather_and_Filter_VCFs.sh [-ms <Gaussian model for snp>] [-mi <Gaussian model for indel>] [-v <VCF directory>]" 2>&1
   echo "-ms: Gaussian model for snp (ex.: 4)"
   echo "-v : vcf file directory (a.k.a variant call output)"
   echo "-mi: Gaussian model for indel (ex.: 4)"
   exit 1
}

if [[ ! $@ =~ ^\-.+ ]]
then
  usage
fi

snp_model=4
indel_model=4
while getopts i:r:b:v:t:g:s:ms:mi: flag
do
    case "${flag}" in
        ms) snp_model=${OPTARG};;
        mi) indel_model=${OPTARG};;
        v) vcf_dir=${OPTARG};;
    esac
done

##Specifying  data directory

ref_dir="/opt/data"
#

###### Combining raw VCFs per chromosome -- VCFs shoulde be ordered
cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  do
    cd $vcf_dir
    ls chr"$i"_part*.vcf.gz > vcf_chr"$i"_list.tsv 
    gatk --java-options "-Xmx80G" GatherVcfs -R $ref_dir/Pf3D7.fasta -I vcf_chr"$i"_list.tsv -O Chr"$i".raw.vcf.gz  --CREATE_INDEX true
    tabix -p vcf Chr"$i".raw.vcf.gz
    rm chr"$i"_part*.vcf.gz
done
########## Normalizing vcf by splitting multiallelic  snps and indels prior to variant recalibration

cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  do
  cd $vcf_dir
  bcftools norm -m-any Chr"$i".raw.vcf.gz | bcftools norm --check-ref -w -f $ref_dir/Pf3D7.fasta | bcftools annotate \
  -Ob -x 'ID' -I +'%CHROM:%POS:%POS:%REF:%ALT' | bcftools view -i 'AC>0' -Oz -o Chr"$i".raw.norm.vcf.gz 
  tabix -p vcf Chr"$i".raw.norm.vcf.gz
done

##### Normalizing vcf to split  snps and indels prior to variant recalibration

####### This script  does raw vcf soft filtering via  Variant Quality Score Recalibration (VQSR) algorithm with our custom training dataset
#    ### VariantRecalibrator for Indels
    cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
 do 
gatk --java-options "-Xmx80g -Xms80g" VariantRecalibrator \
-R $ref_dir/Pf3D7.fasta \
-V Chr"$i".raw.norm.vcf.gz \
--trust-all-polymorphic \
-an QD -an DP -an FS -an SOR -an MQ \
-mode INDEL \
--max-gaussians $indel_model \
-resource:Brown,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz \
-O  Chr"$i".raw.indel.recal \
--tranches-file  Chr"$i".raw.indel.tranches \
--rscript-file  Chr"$i".raw.indel.plots.R

    ### ApplyVQSR for Indels
#
gatk --java-options "-Xmx80g -Xms80g" \
   ApplyVQSR \
    -V Chr"$i".raw.norm.vcf.gz \
    --recal-file Chr"$i".raw.indel.recal \
    --tranches-file Chr"$i".raw.indel.tranches \
    --create-output-variant-index true \
    --lod-score-cutoff -2.0 \
    --exclude-filtered false \
    -mode INDEL \
    -O Chr"$i".raw.indel.recal.vcf.gz
#
# ### VariantRecalibrator for SNPs
gatk --java-options "-Xmx80g -Xms80g" VariantRecalibrator \
-R $ref_dir/Pf3D7.fasta \
-V Chr"$i".raw.indel.recal.vcf.gz \
--trust-all-polymorphic \
-an QD -an DP -an FS -an SOR -an MQ \
-mode SNP \
--max-gaussians $snp_indel \
-resource:Brown,known=true,training=true,truth=true,prior=15.0 Strains.2kb.vcf.gz \
-O Chr"$i".raw.snp.recal \
--tranches-file Chr"$i".raw.snp.tranches \
--rscript-file Chr"$i".raw.snp.plots.R
#
#     #### AppyVQSR for SNPs
gatk --java-options "-Xmx80g -Xms80g" \
   ApplyVQSR \
    -V  Chr"$i".raw.indel.recal.vcf.gz \
    --recal-file  Chr"$i".raw.snp.recal \
    --tranches-file  Chr"$i".raw.snp.tranches \
    --create-output-variant-index true \
    --lod-score-cutoff 0.0 \
    --exclude-filtered false \
    -mode SNP \
    -O  Chr"$i".raw.recal.vcf.gz
done

################## Merging multi-allelic sites after variant recalibration

cd $vcf_dir
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
 do
  gatk SelectVariants \
  -R $ref_dir/Pf3D7.fasta \
  -V Chr"$i".raw.recal.vcf.gz \
  -O Chr"$i".raw.recal.pass.vcf.gz \
  --exclude-filtered true
  bcftools norm -m+any Chr"$i".raw.recal.pass.vcf.gz --check-ref -f $ref_dir/Pf3D7.fasta  -Oz -o Chr"$i".pass.merged.vcf.gz --threads 12
  rm Chr"$i".raw.recal.pass.vcf.gz*
  tabix -p vcf Chr"$i".pass.merged.vcf.gz

done

