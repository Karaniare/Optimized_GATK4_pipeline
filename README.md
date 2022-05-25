# Optimized_GATK4_pipeline
Part1.sh performes fastq and bam processing and quality check. 
Part2.sh performes variant calling, filtering and annotation
The current version of these scripts should be submitted as slurm jobs. 
You need to ensure that you have access to enough cores and memories on your server based on the number of samples and the size of your genome of interest before you run the packages. Downsload the files and follow directions below.

## Specifying necessary directories
This step is required for both packages prior to submitting jobs. Open the package you want to run and edit the directory part accrodingly:

ref_dir="/full/path/to/reference/genome/directory"

fastq_dir="/full/path/to/fastq/directory"

bam_dir="/full/path/to/bam/file/directory"

vcf_dir="/full/path/to/vcfs/directory"

unpaired_dir="/full/path/to/keep/unpaired/fastq/directory"

gvcf_dir="/full/path/to/gvcf/directory"

stat_dir="/full/path/to/qc/directory"

## Creating sample ID list 
This step is required for both packages prior to submitting jobs.

Create a sample ID file and call it ID_list.tsv in your bam file directory.

## Submitting slurm jobs
sbatch Part1.sh

sbatch Part2.sh
