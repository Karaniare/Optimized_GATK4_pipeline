#!/bin/bash
#SBATCH -J Variant_calling
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g
#
###### Creating arguments
usage() {
  echo "Usage: $0 [-i <input>] [-v <VCF directory>]" 1>&2
  echo "-i: list of VCFs file name to annotate (only names without .vcf.gz) in tab separated value format"
  echo "-v : VCF file directory (directory of variant call outputs)"
  exit 1
}

if [[ ! $@ =~ ^\-.+ ]]
then
  usage
fi


while getopts i:b:v:g:c:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        v) vcf_dir=${OPTARG};;
       
 esac
done


##### Starts annotating VCF
cd $vcf_dir
for i in $(cat $input)
 do
 cd /opt/resources/snpEff/data

  java -Xmx80g -jar snpEff.jar Pfalciparum.noseq $vcf_dir/"$i".vcf.gz -no-upstream -no-downstream | bgzip -c > $vcf_dir/"$i".annotated.vcf.gz
  tabix -p vcf $vcf_dir/"$i".annotated.vcf.gz
done
