#!/bin/bash
#SBATCH -J Variant_calling
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g
#
###### Creating arguments
usage() {
  echo "Usage: $0 or bash /opt/resources/Making_gVCFs.sh [-i <input>] [-c <chromosome number>] [-b <bam directory>] [-g <gVCF directory>]" 1>&2
  echo "-i: input file (sample ID list) in tab separated value format"
  echo "-c: chromosome number (ex.: 1)" 
  echo "-b: bam file directory"
  echo "-g: directory to save gvcf files"
  exit 1
}

if [[ ! $@ =~ ^\-.+ ]]
then
  usage
fi


while getopts i:c:b:g: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        c) chrom=${OPTARG};;
        b) bam_dir=${OPTARG};;
        g) gvcf_dir=${OPTARG};;
    esac
done
echo "Argument definition:"
echo "-i: input file (sample ID list) in tab separated format"
echo "-t: genomic region file"
echo "-c: chromosome number"
echo "-b: bam file directory"
echo "-v : vcf file (or output) directory"
echo "-g: gvcf file directory"
region="/opt/data"
#ref_dir="/opt/data"
ref_dir="/users/kniare/data/shared/scripts/Data_WGS_UNC"

######### Variant Calling starts here
######## Running HaplotypeCaller to generate gVCFs
cd $gvcf_dir
mkdir chr"$chrom"
cd $bam_dir
for i in $(cat $input)
    do 
      for j in $chrom
       do
          cd $gvcf_dir
          gatk --java-options "-Xmx80G" HaplotypeCaller -R $ref_dir/Pf3D7.fasta -I "$bam_dir/$i".sorted.dup.pf.bam -ERC GVCF -ploidy 2 \
         --native-pair-hmm-threads 16 -O  chr"$j"/"$i".chr"$j".g.vcf --assembly-region-padding 100 \
         --max-num-haplotypes-in-population 128 --kmer-size 10 --kmer-size 25 \
         --min-dangling-branch-length 4 --heterozygosity 0.0029 --indel-heterozygosity 0.0017 --min-assembly-region-size 100 -L $ref_dir/core_chr"$j".list -mbq 5 -DF MappingQualityReadFilter --base-quality-score-threshold 12
       done
done

