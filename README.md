# Optimized_GATK4_pipeline
QC_Pf_WGS.sh performes fastq and bam processing and quality check. 
The scripts Making_gVCFs.sh, Making_VCFs.sh, Gather_and_Filter_VCFs and Annotating_VCFs.sh perform variant calling in the core genome of _Plasmodium falciparum_ and filter and annotate the VCFs
Each script of the current version of the pipeline should be submitted as slurm jobs. 
It is important to ensure that enough cores and memories are accessible on the server based on the number of samples and the size of the genome of interest before you run the packages. Clone the files and follow directions below.

## Dependencies

* The pipeline requires GATK4, bwa, trimmomatic, samtools and datamash and need to be installed. If they are already installed on your server as modules, you do not have to load them as the pipeline does it automatically.

## Running optimized GATK4 pipeline
* Download Part1.sh and Part2.sh packages and execute them with arguments as specified below. For Part2.sh, you need to downsload SnpEff folder and keep it in the vcf (output) directory.

*sh Part1.sh -i [input list] -f [fastq directory] -b [bam directory] -u [unpaired directory] -s [stat directory] -k [kit name]
   - -i : list of sample IDs in a text file
   - -f: full path to the fastq file
   - -u: directory to keep unpaired fastqs after trimming
   - -s: directory to keep QC stats
   - -k: name of the library prep kit (Examples: TruSeq3-PE.fa, /opt/data/Nextera-PE.fa)

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
