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

**QC_Pf_WGS.sh**
* sh QC_Pf_WGS.sh [-i <input list>] [-f <fastq directory>] [-b <bam directory>] [-u <unpaired directory>] [-s <stat directory>] [-k <kit name>].
 
-  -i : list of sample IDs in a text file.
-  -f: full path to the fastq files.
-  -u: directory to keep unpaired fastqs after trimming.
-  -s: directory to keep QC stats.
-  -k: name of the library prep kit (Ex.: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa).

Making_gVCFs.sh
* sh Making_gVCFs.sh [-i <input>] [-c <chromosome number>] [-b <bam directory>] [-g <gVCF directory>].
- -i: input file (sample ID list) in tab separated value format.
- -c: chromosome number (ex.: 1).
- -b: bam file directory.
- -g: directory to save gVCF files.

