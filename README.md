# Optimized_GATK4_pipeline
**QC_Pf_WGS.sh** performes fastq and bam processing and quality check. 
The scripts **Making_gVCFs.sh**, **Making_VCFs.sh**, **Gather_and_Filter_VCFs.sh** and **Annotating_VCFs.sh** perform variant calling in the core genome of _Plasmodium falciparum_ and filter and annotate the VCFs.
Each script of the current version of the pipeline should be submitted indepently as slurm jobs. 
It is important to ensure that enough cores and memories are accessible on the server based on the number of samples and the size of the genome of interest before you run the packages. Clone the files and follow directions below.
A singularity container containing all these packages is also available (highly recommended). Just run wget https://baileylab.brown.edu/opti_gatk4/opti_gatk4_230911.sif to download the container. A manual on how to run the analysis in the container can be found in this link https://docs.google.com/document/d/1SYYFo9kWHRwo20rhY1VU84RPaQKg8zeXyM6Cy7DJ4Y0/edit

## Dependencies

* The pipeline requires GATK4, bwa, trimmomatic, samtools and datamash. If they are already installed on your server as modules, you do not have to load them as the pipeline does it automatically. Some additional resources are needed including Pf3D7 reference genome and index files, homo sapiens reference genome and index files, concatenaed Pf3D7 and human (hibrid) reference genome and index files, Training dataset (VCF), bed files specifying core regions of chromosomes in gatk format (.list), bed files specifying Pf3D7 and human genome chromosome locationsin samtools format (tab-separated value format) and snpeff database with the Pf3D7_v3 assembly. All the tools and resources are available in the the singularity container.

## Running optimized GATK4 pipeline
* Please use the packages with arguments as specified below. 

**QC_Pf_WGS.sh** (does quality control)
* sh QC_Pf_WGS.sh [-i <input list>] [-f <fastq directory>] [-b <bam directory>] [-u <unpaired directory>] [-s <stat directory>] [-k <kit name>] [-r <run name>].
 
  -  -i : list of sample IDs in a text file.
  -  -f: full path to the fastq files.
  -  -u: directory to keep unpaired fastqs after trimming.
  -  -s: directory to keep QC stats.
  -  -b: bam file directory
  -  -k: name of the library prep kit (Ex.: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa).
  -  -r: run name (Ex: run1), useful if you run the pipeline on separate samplesets simultaneously

 
**Making_gVCFs.sh** (generates gVCFs based on HaplotypeCaller).
* sh Making_gVCFs.sh [-i <input>] [-c <chromosome number>] [-b <bam directory>] [-g <gVCF directory>].
  - -i: input file (sample ID list) in tab separated value format.
  - -c: chromosome number (ex.: 1).
  - -b: bam file directory.
  - -g: directory to save gVCF files.


**Making_VCFs.sh** (generates VCFs from gVCFs)
* sh Making_VCFs.sh [-v <VCF directory>] [-c <chromosome number>] [-t <region>].
   - -c: chromosome number (ex.: 1).
   - -v: directory to save VCF files.
   - -t: file containing a list of genomic regions to call variants at, one region per line (ex.: Pf3D7_01_v3:1-20000  Pf3D7_01_v3:20001-40000 Pf3D7_01_v3:40001-60000 Pf3D7_01_v3:60001-80000).
 
 * **NB:** files specifying the gVCF IDs (first column) and the full paths to them (second column) should be created in the VCF directory by the user. They should be named as follows: gvcf_chr1_list.tsv gvcf_chr2_list.tsv gvcf_chr3_list.tsv gvcf_chr4_list.tsv gvcf_chr5_list.tsv gvcf_chr6_list.tsv gvcf_chr7_list.tsv gvcf_chr8_list.tsv gvcf_chr9_list.tsv gvcf_chr10_list.tsv gvcf_chr11_list.tsv gvcf_chr12_list.tsv gvcf_chr13_list.tsv and gvcf_chr14_list.tsv
 
 **Gather_and_Filter_VCFs**.sh (joins small VCFs and recalibrates variants)
 * sh Gather_and_Filter_VCFs.sh [-ms <Gaussian model for snp>] [-mi <Gaussian model for indel>] [-v <VCF directory>].
    - -ms: Gaussian model for snp (ex.: 4).
    - -v : vcf file directory (a.k.a variant calling output).
    - -mi: Gaussian model for indel (ex.: 4).
 
 
**Annotating_VCFs.sh** (performs functional annotation of VCFs)
 * sh  Annotating_VCFs.sh [-i <input>] [-v <VCF directory>].
   - -i: list of VCFs file name to annotate (only names without .vcf.gz) in tab separated value format and should be in the VCF directory.
   - -v : VCF file directory (directory of variant call outputs).
