###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
import shutil, time, glob, sys, os, re


inputFolder = "/shared/projects/jill_panel_pangaea"
fastQscreen_config = "/home/stavros/playground/16S_metagenomics/subsidiary_files/fastq_screen.conf"
refGenome_GRCh38 = "/home/stavros/playground/progs/reference_files/reference_genome/GRCh38_primAssembly"
reference_annotation = "/home/stavros/playground/progs/reference_files/gene_annotation/gencode.v29.primary_assembly.annotation.gtf"

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
parser.add_argument('-th', '--threads', dest='threads', default=20, metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()



args.project_dir = inputFolder


### Main directories
analysisDir = os.path.join(os.getcwd(), "analysis")
qcDir = os.path.join(analysisDir, "qc_analysis")
alignment_dir = os.path.join(analysisDir, "alignments")
temp = os.path.join(analysisDir, "temp")


if not os.path.exists(qcDir): os.makedirs(qcDir)

def quality_control():

	if not os.path.exists(qcDir): os.makedirs(qcDir)

	for path, subdir, folder in os.walk(args.project_dir):
		for dirs in subdir:
			if dirs == "Tanda_106":
				mfiltered_data = ' '.join([f for f in glob.glob(os.path.join(os.path.join(path, dirs), "*.fastq.gz"))])

				print("Checking random reads for possible contamination: in progress ..")
				fastQscreen = ' '.join([
				"fastq_screen",  # Call fastQ screen to check contamination in the processed data
				"--threads", args.threads,  # Number of threads to use
				"--outdir",  qcDir,  # Directory in which the output files will be saved
				"--quiet",  # Suppress all progress reports on stderr and only report errors
				"--conf", fastQscreen_config,  # Location of the required configuration file
				mfiltered_data, # mfiltered_data.replace("_R1_", "_R2_"),  # Input PE files
				"2>>", os.path.join(qcDir, "{0}_fastQscreen_report.txt".format(dirs))])  # Output fastQ screen report
				# subprocess.run(fastQscreen, shell=True)
			
				print("Quality Control reports are being generated: in progress ..")
				fastQC = ' '.join([
				"fastqc",  # Call fastQC to quality contol all processed data
				"--threads", args.threads,  # Number of threads to use
				"--quiet",  # Print only log warnings
				"--outdir", qcDir,  # Create all output files in this specified output directory
				mfiltered_data,  # String containing all samples that are about to be checked
				"2>>", os.path.join(qcDir, "{0}_fastQC_report.txt".format(dirs))])  # Output fastQC report
				subprocess.run(fastQC, shell=True)

				print("Quality Control of all PE reads: in progress ..")
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
					"--html", os.path.join(qcDir, "{0}_fastp.html".format(os.path.basename(files)[:-9])),  # Create ftml file in this specified output directory
					"--json", os.path.join(qcDir, "{0}_fastp.json".format(os.path.basename(files)[:-9])),  # Create json output file in this specified output directory
					"2>>", os.path.join(qcDir, "{0}_fastP_report.txt".format(dirs))])  # Output fastP report
					subprocess.run(fastP, shell=True) 

				multiQC = " ".join([
				"/opt/anaconda3/bin/multiqc",  # Call MultiQC
				"--quiet",  # Print only log warnings
				"--outdir", qcDir,  # Create report in the FastQC reports directory
				"--filename", "summarised_initqc_report",  # Name of the output report 
				qcDir,  # Directory where all FastQC and Cutadapt reports reside
				"2>>", os.path.join(qcDir, "{0}_multiQC_report.txt".format(dirs))])  # Output multiQC report
				subprocess.run(multiQC, shell=True)


	return



def main():
	
	quality_control()
	
	# Preprocessing of the data if it isn't already being done
	# if len(glob.glob("{0}/*_trimmed.fq.gz".format(preprocessedFiles))) == 0:
	# 	preprocessing_samples()  # Preprocessing the raw reads


if __name__ == "__main__": main()