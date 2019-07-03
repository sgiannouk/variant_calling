###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
import shutil, time, glob, sys, os, re


inputFolder = "/shared/projects/jill_dnapanel_pangaea"
fastQscreen_config = "~/playground/progs/16S_subsidiary_files/fastq_screen.conf"
refGenome_GRCh38 = "~/playground/progs/reference_files/reference_genome/GRCh38_primAssembly/GRCh38_primary_assembly_genome.fa"
reference_annotation_gtf = "~/playground/progs/reference_files/gene_annotation/gencode.v29.primary_assembly.annotation.gtf"
reference_annotation_bed = "~/playground/progs/reference_files/gene_annotation/hg38_Gencode_V28.bed"

usage = "variant_calling [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Create required section in help
requiredArgs = parser.add_argument_group('required arguments')
# Input folder option
requiredArgs.add_argument('-i', '--project_dir', required=False, metavar='', 
					   	help="Path of the input directory that contains the raw data.")
# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=str(30), metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()



args.project_dir = inputFolder


### Main directories

analysis_dir = os.path.join(os.getcwd(), "analysis")
temp = os.path.join(analysis_dir, "temp")
reports_dir = os.path.join(analysis_dir, "preprocessing_reports")
preprocesed_dir = os.path.join(analysis_dir, "preprocesed_data")
alignment_dir = os.path.join(analysis_dir, "alignments")
varcall_dir = os.path.join(analysis_dir, "variant_calling")
postalign_dir = os.path.join(analysis_dir, "post-alignment_reports")


def quality_control():

	if not os.path.exists(reports_dir): os.makedirs(reports_dir)
	if not os.path.exists(temp): os.makedirs(temp)

	for path, subdir, folder in os.walk(args.project_dir):
		for i, dirs in enumerate(subdir, 0):
			
			print("{0}/{1} | Analysing {2}".format(i, len(subdir), dirs))
			mfiltered_data = ' '.join([f for f in glob.glob(os.path.join(os.path.join(path, dirs), "*.fastq.gz"))])

			print("fastQscreen - Checking random reads for possible contamination: in progress ..")
			fastQscreen = ' '.join([
			"fastq_screen",  # Call fastQ screen to check contamination in the processed data
			"--threads", args.threads,  # Number of threads to use
			"--outdir",  temp,  # Directory in which the output files will be saved
			"--quiet",  # Suppress all progress reports on stderr and only report errors
			"--conf", fastQscreen_config,  # Location of the required configuration file
			mfiltered_data,
			"2>>", os.path.join(temp, "{0}_fastQscreen_init_report.txt".format(dirs))])  # Output fastQ screen report
			subprocess.run(fastQscreen, shell=True)
		
			print("fastQC - Quality Control reports are being generated: in progress ..")
			fastQC = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--quiet",  # Print only log warnings
			"--outdir", temp,  # Create all output files in this specified output directory
			mfiltered_data,  # String containing all samples that are about to be checked
			"2>>", os.path.join(temp, "{0}_fastQC_init_report.txt".format(dirs))])  # Output fastQC report
			subprocess.run(fastQC, shell=True)

			print("fastP - Quality Control of all reads: in progress ..")
			for files in glob.glob(os.path.join(dirs, "*.fastq.gz")):
				fastP = ' '.join([
				"fastp",  # Call fastQC to quality control all processed data
				"--thread", args.threads,  # Number of threads to use
				"--in1", files,  # Input read1 file
				"--disable_adapter_trimming",  # Adapter trimming is disabled
				"--disable_trim_poly_g",  # Disable polyG tail trimming
				"--disable_quality_filtering",  # Quality filtering is disabled
				"--disable_length_filtering",  # Length filtering is disabled
				"--overrepresentation_analysis",  # Enable overrepresented sequence analysis
				"--html", os.path.join(temp, "{0}_fastp.html".format(os.path.basename(files)[:-9])),  # Create ftml file in this specified output directory
				"--json", os.path.join(temp, "{0}_fastp.json".format(os.path.basename(files)[:-9])),  # Create json output file in this specified output directory
				"2>>", os.path.join(temp, "{0}_fastP_init_report.txt".format(dirs))])  # Output fastP report
				subprocess.run(fastP, shell=True) 

			print("multiQC - Summing all QC reports: in progress ..\n")
			multiQC = " ".join([
			"/home/stavros/anaconda3/bin/multiqc",  # Call MultiQC
			"--quiet",  # Print only log warnings
			"--outdir", temp,  # Create report in the FastQC reports directory
			"--filename", "{0}_init_summarised_report".format(dirs),  # Name of the output report 
			temp,  # Directory where all FastQC and Cutadapt reports reside
			"2>>", os.path.join(temp, "{0}_multiQC_init_report.txt".format(dirs))])  # Output multiQC report
			subprocess.run(multiQC, shell=True)

			os.system('mv {0}/* {1}'.format(temp, reports_dir))
	
	for path, subdir, folder in os.walk(reports_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or (name.endswith("multiQC_report.txt") and os.stat(file).st_size == 583):  
				os.remove(file)
	
	os.system('rm {0}/*fastqc.zip'.format(reports_dir))
	os.system('rm {0}/*screen.png'.format(reports_dir))
	os.system('rm {0}/*screen.html'.format(reports_dir))
	os.system('rm {0}/*fastQC_init_report.txt'.format(reports_dir))
	os.system('rm -r {0}/*init_summarised_report_data'.format(reports_dir))
	return

def preprocessing_rawdata():
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whos quality is below Q20). That way we remove obvious trash and 
	trying to improve mapping. """
	if not os.path.exists(preprocesed_dir): os.makedirs(preprocesed_dir)
	raw_data = glob.glob("{0}/*/*.fastq.gz".format(args.project_dir))

	print("BBDUK - Quality filtering of the raw data: in progress...")
	for file in raw_data:
		# Obtaining the sample name
		sample_name = os.path.basename(file).split(".")[0]
		bbduk = ' '.join([
		"bbduk.sh",  # Call BBDuck (BBTools) to preprocess the raw data
		"threads={0}".format(args.threads),  # Set number of threads to use
		"in={0}".format(file),  # Input of the forward file
		"out={0}".format(os.path.join(preprocesed_dir, "{0}.qt.fastq.gz".format(sample_name))),  # Export edited forward read to file
		"trimq=20",  # Regions with average quality BELOW this will be trimmed (20)
		"qtrim=r",  # Trim read ends to remove bases with Q<15
		"ordered=t",  # Keeps the reads in the same order as we gave them to the software
		"maq={0}".format(20),  # Reads with average quality (after trimming) below this will be discarded (15)
		"2>>", os.path.join(reports_dir, "bbduk_qtrim_report.txt")])  # Output trimming report
		subprocess.run(bbduk, shell=True)
	return

def aligning():
	if not os.path.exists(alignment_dir): os.makedirs(alignment_dir)
	
	preprocesed_data = glob.glob("{0}/*.qt.fastq.gz".format(preprocesed_dir))
	
	print("minimap2 - aligning in total {0} samples against the reference genome...".format(len(preprocesed_data)))
	for i, file in enumerate(preprocesed_data):
		file_name = os.path.basename(file).split(".")[0]
		minimap2_genome = " ".join([
		"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
		"-t", args.threads,  # Number of threads to use
		"-ax sr",   # Genomic short-read mapping mode and output in SAM format (-a)
		# "-o", os.path.join(alignment_dir, "{0}.genome.paf".format(file_name)),
		refGenome_GRCh38,  # Inputting the reference genome
		file,  # Input .fastq.gz file
		"|", "samtools view",  # Calling 'samtools view' to compress the Bowtie's output file to BAM
  		"--threads", args.threads,  # Number of threads to be used by Samtools in the conversion of the SAM files to BAM
  		"-S -u1",  # Input format is auto-detected, the output file should be in BAM format, and use fast compression
  		"-",  # Piping the input file
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
  		"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
  		"-o", os.path.join(alignment_dir, "{0}.qt.genome.bam".format(file_name)), "-",  # Sorted output  BAM file
		"2>>", os.path.join(reports_dir, "minimap2_genome_report.txt")])  # minimap2 report
		subprocess.run(minimap2_genome, shell=True)

	
	# Post Alignment QC
	if not os.path.exists(postalign_dir): os.makedirs(postalign_dir)
	aligned_data = glob.glob("{0}/*.qt.genome.bam".format(alignment_dir))

	print("fastQC - Alignment Quality Control reports are being generated: in progress ..")
	fastQC_align = ' '.join([
	"fastqc",  # Call fastQC to quality control all processed data
	"--threads", args.threads,  # Number of threads to use
	"--quiet",  # Print only log warnings
	"--outdir", temp,  # Create all output files in this specified output directory
	' '.join(aligned_data),  # String containing all samples that are about to be checked
	"2>>", os.path.join(temp, "fastQC_alignment_report.txt")])  # Output fastQC report
	subprocess.run(fastQC_align, shell=True)

	print("RSeQC & Picard - post alignment stats: in progress ..")
	for i, file in enumerate(aligned_data):
		file_name = os.path.basename(file).split(".")[0]
		
		bam_stat = ' '.join([
		"bam_stat.py",  # Call samtools flagstat
		"-i", file,  # Input BAM file
		"> {0}/{1}_bam_stat.txt".format(temp, file_name),  # Output file
		"2>>", os.path.join(temp, "{0}_bamstat_report.txt".format(file_name))])
		# subprocess.run(bam_stat, shell=True)

		read_distribution = ' '.join([
		"read_distribution.py",  # Call samtools flagstat
		"-i", file,  # Input BAM file
		"-r", reference_annotation_bed,
		"> {0}/{1}.fragSize".format(temp, file_name),  # Output file
		"2>>", os.path.join(temp, "{0}_read_distribution_report.txt".format(file_name))])
		# subprocess.run(read_distribution, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}_alignment_metrics.txt".format(temp, file_name),  # Output
		"REFERENCE_SEQUENCE= {0}".format(refGenome_GRCh38),  # Reference sequence file
		"2>>", os.path.join(temp, "{0}_CollectAlignmentSummaryMetrics_report.txt".format(file_name))])
		# subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 

	print("multiQC - Summing all QC reports: in progress ..\n")
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", temp,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	temp,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(temp, "post-alignment_multiQC_report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)
	

	for path, subdir, folder in os.walk(temp):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or (name.endswith("multiQC_report.txt") and os.stat(file).st_size == 583):  
				os.remove(file)
	
	os.system('rm {0}/*fastqc.zip'.format(reports_dir))
	os.system('rm -r {0}/*summarised_report_data'.format(temp))
	os.system('mv {0}/* {1}'.format(temp, postalign_dir))
	return

def variant_calling():

	if not os.path.exists(varcall_dir): os.makedirs(varcall_dir)
	aligned_data = glob.glob("{0}/*.qt.genome.bam".format(alignment_dir))
	
	
	for i, file in enumerate(aligned_data):

		print("{0}/{1} | Variant Calling {2} using GATK4...".format(i, len(aligned_data), file))
		file_name = os.path.basename(file).split(".")[0]
		# Marking duplicate reads
		markDup = " ".join([
		"gatk MarkDuplicates",  # Call gatk MarkDuplicates(v4.1.2.0)
		"--CREATE_INDEX:true",  # Create a BAM index when writing a coordinate-sorted BAM file
		"--INPUT: {0}".format(file),  # BAM files to analyze
		"--VERBOSITY: ERROR",  # Control verbosity of logging
		"--VALIDATION_STRINGENCY: STRICT",  # Validation stringency for all BAM files read by this program
		"--METRICS_FILE: {0}/markDup_{1}_metrics.txt".format(reports_dir, file_name),
		"--OUTPUT: {0}/{1}.mdqt.genome.bam".format(alignment_dir, file_name),
		"2>>", os.path.join(reports_dir, "markDup_genome_report.txt")])
		subprocess.run(markDup, shell=True)



	return

def main():
	

	# quality_control()

	# preprocessing_rawdata()
	
	# aligning()

	variant_calling()

if __name__ == "__main__": main()