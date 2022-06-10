# Optimized_GATK4_pipeline
Part1.sh performes fastq and bam processing and quality check. 
Part2.sh performes variant calling in the core genome of _Plasmodium falciparum_ and does vcf filtering and annotation
The current version of these scripts should be submitted as slurm jobs. 
You need to ensure that you have access to enough cores and memories on your server based on the number of samples and the size of your genome of interest before you run the packages. Downsload the files and follow directions below.

## Dependencies

* The pipeline requires GATK4, bwa, trimmomatic, samtools and datamash and need to be installed. If they are already installed on your server as modules, you do not have to load them as the pipeline does it automatically.

## Running optimized GATK4 pipeline
* Download Part1.sh and Part2.sh packages and execute them with arguments as specified below. For Part2.sh, you need to downsload SnpEff folder and keep it in the vcf (output) directory.
* sh Part1.sh -i [OPTION] -f [OPTION] -r [OPTION] -b [OPTION] -u [OPTION] -s [OPTION]
* Argument definitions:

   - -i: input file (sample ID list) in tab separated format
   - -f: fastq file directory
   - -r: reference genome and bed file directory
   - -b: bam file directory
   - -v : vcf file (or output) directory
   - -u: directory where unpaired fastq files will be kept
   - -g: gvcf file directory
   - -s: directory to save qc output files

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
