#!/bin/bash
#SBATCH -J Variant_calling
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=60g
#
###### Creating arguments
usage() {
  "Usage:$0 or bash /opt/resources/Making_VCFs.sh [-i <input>] [-b <bam directory>] [-v <VCF directory>] [-g <gVCF directory>] [-c <chromosome number>] [-t <region>]" 1>&2
  echo "-i: input file (sample ID list) in tab separated value format"
  echo "-c: chromosome number (ex.: 1)"
  echo "-b: bam file directory"
  echo "-v: directory to save VCF files"
  echo "-g: gVCF file directory"
  echo "-t: file containing a list of genomic regions to call variants at, one region per line (ex.: 1-20,000  20,001-40,000 40,001-60,000 60,001-80,000)"
  exit 1
}


while getopts i:b:v:g:c:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        b) bam_dir=${OPTARG};;
        v) vcf_dir=${OPTARG};;
        g) gvcf_dir=${OPTARG};;
	c) chrom=${OPTARG};;
	t) region=${OPTARG};;
    esac
done
echo "##################This package works on a server with slurm job submissions#################"
echo "Argument definition:"
echo "-i: input file (sample ID list) in tab separated format"
echo "-t: genomic region file"
echo "-b: bam file directory"
echo "-v : vcf file (or output) directory"
echo "-g: gvcf file directory"
echo "-t: specify a file containing a list of genomic regions to call variants such as 1-20,000  20,001-40,000 40,001-60,000 60,001-80,000"
###
ref_dir="/opt/data"

######## Combining gVCFs into a genomic database 
###### Combining gvcfs 
cd $vcf_dir
for i in $chrom
  do
    gatk --java-options "-Xmx40G" GenomicsDBImport \
    --sample-name-map gvcf_chr"$i"_list.tsv \
    --genomicsdb-workspace-path chr"$i"_database \
    --intervals $ref_dir/core_chr"$i".list --batch-size 100 \
    --reader-threads 24 --genomicsdb-segment-size 8048576 \
    --genomicsdb-vcf-buffer-size 160384
done
##### Running joint genotyping by genomic segment (the core genome of each chromosome is split into smaller genomic regions) via multiple slrum jobs
cd $vcf_dir
for i in $chrom
  do
   cd $ref_dir
   for j in $(cat $region)
     do
     cd $vcf_dir
     echo -e "#!/bin/bash" > genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH -J Genotype_chr"$i"_"$j"" >> genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH -t 120:00:00" >> genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH -c 8" >> genotype_chr"$i"_"$j".sh
     echo -e "#SBATCH --mem-per-cpu=10g" >> genotype_chr"$i"_"$j".sh
     echo -e "module load gatk/4.1.6.0 samtools/ bcftools/1.9 bcftools/1.9" >> genotype_chr"$i"_"$j".sh
     echo -e "ulimit -c unlimited" >> genotype_chr"$i"_"$j".sh
     echo -e "module load java/jdk-17.0.2" >> genotype_chr"$i"_"$j".sh
     echo -e "gatk --java-options "-Xmx80g -Xms80g"  GenotypeGVCFs --genomicsdb-use-bcf-codec true  -R $ref_dir/Pf3D7.fasta -V gendb://$vcf_dir/chr"$i"_database --max-genotype-count 1024 -O $vcf_dir/chr"$i"_part"$j".vcf.gz --tmp-dir $ref_dir -stand-call-conf 30 -L "$j"" >> genotype_chr"$i"_"$j".sh
     chmod +x genotype_chr"$i"_"$j".sh
     sed -i 's/-Xmx80g -Xms80g/"-Xmx80g -Xms80g"/g' genotype_chr"$i"_"$j".sh
     sbatch genotype_chr"$i"_"$j".sh
   done
done

