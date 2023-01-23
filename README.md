# Optimized_GATK4_pipeline
**QC_Pf_WGS.sh** performes fastq and bam processing and quality check. 
The scripts **Making_gVCFs.sh**, **Making_VCFs.sh**, **Gather_and_Filter_VCFs.sh** and **Annotating_VCFs.sh** perform variant calling in the core genome of _Plasmodium falciparum_ and filter and annotate the VCFs
Each script of the current version of the pipeline should be submitted as slurm jobs. 
It is important to ensure that enough cores and memories are accessible on the server based on the number of samples and the size of the genome of interest before you run the packages. Clone the files and follow directions below.
A singularity container containing all these packages is also available.

## Dependencies

* The pipeline requires GATK4, bwa, trimmomatic, samtools and datamash. If they are already installed on your server as modules, you do not have to load them as the pipeline does it automatically. Some additional resources are needed including Pf3D7 reference genome and index files, homo sapiens reference genome and index files, concatenaed Pf3D7 and human (hibrid) reference genome and index files, Training dataset (VCF), bed files specifying core regions of chromosomes in gatk format (.list), bed files specifying Pf3D7 and human genome chromosome locationsin samtools format (tab-separated format) and snpeff database with the Pf3D7_v3 assembly. All the tools and resources are available in the the singularity container.

## Running optimized GATK4 pipeline
* If want to use the non-container version, download each package packages and execute them with arguments as specified below. 
* 
* usage() { echo -e "Usage: $0 [-i <input list>] [-f <fastq directory>] [-b <bam directory>] [-u <unpaired directory>] [-s <stat directory>] [-k <kit name>]" 2>&1
 echo -e "-i : list of sample IDs in a text file "
 echo -e "-f: full path to the fastq files"
 echo -e "-u: directory to keep unpaired fastqs after trimming "
 echo -e "-s: directory to keep QC stats"
 echo -e "-k: name of the library prep kit (Ex.: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa)"

}
* sh Part2.sh -i [OPTION] -t [OPTION] -r [OPTION] -b [OPTION] -v [OPTION] -g [OPTION]
* Argument definitions:
   - -i: input file (sample ID list) in tab separated format
   - -t: genomic region file
   - -r: reference genome and bed file directory
   - -b: bam file directory
   - -v : vcf file (or output) directory
   - -g: gvcf file directory
* For Part2.sh, you need to create a gvcf list for each chromosome in the vcf (output) directory. 
   - cd vcf_directory
   - awk -F "\t" -v OFS="\t" '{print $1,"/full/path/to/the/chr1/gvcf/files/"$1}' sample_ID_list.tsv > gvcf_chr1_list.tsv
