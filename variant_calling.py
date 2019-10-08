###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.1"

import argparse, subprocess
from datetime import datetime
import shutil, time, glob, sys, os, re


inputFolder = "/shared/projects/jill_dnapanel_pangaea"
fastQscreen_config = "/home/stavros/references/fastQscreen_references/fastq_screen.conf"
### REFERENCE FILES
refGenome_GRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
reference_annotation_gtf = "/home/stavros/references/reference_annotation/GRCh38_gencode.v31.primAssembly_psudo_trna.annotation.gtf.gz"
reference_annotation_bed = "/home/stavros/references/reference_annotation/hg38_gencode.v31.allComprehensive_pseudo.annotation.bed"
refGenome_GRCh38sdf = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome"
### GATK FILES
gatk_functional_annot = "/home/stavros/references/gatk_subfiles/functional_annot/funcotator_dataSources.v1.6.20190124s"
gatk_wg_intervals = "/home/stavros/references/gatk_subfiles/intervals/wgs_calling_regions.hg38.interval_list"
gatk_liftover_GRCh38 = "/home/stavros/references/gatk_subfiles/liftover_GRCh38/af-only-gnomad.hg38.vcf.gz"
gatk3_8 = "java -jar /home/stavros/playground/progs/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar"
gatk_pon = "/home/stavros/references/gatk_subfiles/panel_of_normals/1000g_pon.hg38.vcf.gz"
gatk_known_sites = "/home/stavros/references/gatk_subfiles/known_sites"
starfish = "python3 /home/stavros/playground/progs/starfish/starfish.py"



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
parser.add_argument('-th', '--threads', dest='threads', default=str(10), metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()
args.project_dir = inputFolder


### Directories
analysis_dir = os.path.join(os.getcwd(), "vc_analysis")
preprocesed_dir = os.path.join(analysis_dir, "preprocesed_data")
prereports_dir = os.path.join(analysis_dir, "preprocessing_reports")
alignment_dir = os.path.join(analysis_dir, "alignments")
postprocesed_dir = os.path.join(analysis_dir, "postprocessed_data")
postreports_dir = os.path.join(analysis_dir, "postprocessing_reports")
varcall_dir = os.path.join(analysis_dir, "variant_calling")
prefilter_plots = os.path.join(varcall_dir, "prefilter_plots")
filtered_plots = os.path.join(varcall_dir, "filtered_plots")
intersection_analysis = os.path.join(varcall_dir, "intersection_analysis")
temp = os.path.join(analysis_dir, "temp")


def quality_control():
	print("{0} PREPROCESSING & QUALITY CONTROL".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	if not os.path.exists(prereports_dir): os.makedirs(prereports_dir)
	if not os.path.exists(preprocesed_dir): os.makedirs(preprocesed_dir)
	if not os.path.exists(temp): os.makedirs(temp)
	
	
	for path, subdir, folder in os.walk(args.project_dir):
		for i, dirs in enumerate(subdir, 0):
			if dirs == "Tanda_113":
				print("{0} {1}/{2} | Analysing {3}".format(datetime.now().strftime("%d.%m %H:%M"),i, len(subdir), dirs))
				mfiltered_data = [f for f in glob.glob(os.path.join(os.path.join(path, dirs), "*.fastq.gz"))]

				print(">> {0} 1/5 | fastQscreen - Checking random reads for possible contamination: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
				fastQscreen = ' '.join([
				"fastq_screen",  # Call fastQ screen to check contamination in the processed data
				"--threads", args.threads,  # Number of threads to use
				"--outdir",  temp,  # Directory in which the output files will be saved
				"--quiet",  # Suppress all progress reports on stderr and only report errors
				"--conf", fastQscreen_config,  # Location of the required configuration file
				' '.join(mfiltered_data),
				"2>>", os.path.join(temp, "{0}_fastQscreen_report.txt".format(dirs))])  # Output fastQ screen report 
				subprocess.run(fastQscreen, shell=True)
			
				print(">> {0} 2/5 | fastQC - Quality Control reports are being generated: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
				fastQC = ' '.join([
				"fastqc",  # Call fastQC to quality control all processed data
				"--threads", args.threads,  # Number of threads to use
				"--quiet",  # Print only log warnings
				"--outdir", temp,  # Create all output files in this specified output directory
				' '.join(mfiltered_data),  # String containing all samples that are about to be checked
				"2>>", os.path.join(temp, "{0}_fastQC_init_report.txt".format(dirs))])  # Output fastQC report
				subprocess.run(fastQC, shell=True)

				print(">> {0} 3/5 | fastP - Quality Control of all reads: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
				for file in mfiltered_data:
					sample_name = os.path.basename(file).split(".")[0]
					fastP = ' '.join([
					"fastp",  # Call fastQC to quality control all processed data
					"--thread", args.threads,  # Number of threads to use
					"--in1", file,  # Input read1 file
					"--disable_adapter_trimming",  # Adapter trimming is disabled
					"--disable_trim_poly_g",  # Disable polyG tail trimming
					"--disable_quality_filtering",  # Quality filtering is disabled
					"--disable_length_filtering",  # Length filtering is disabled
					"--overrepresentation_analysis",  # Enable overrepresented sequence analysis
					"--html", os.path.join(temp, "{0}_fastp.html".format(sample_name)),  # Create ftml file in this specified output directory
					"--json", os.path.join(temp, "{0}_fastp.json".format(sample_name)),  # Create json output file in this specified output directory
					"2>>", os.path.join(temp, "{0}_fastP_report.txt".format(sample_name))])  # Output fastP report
					subprocess.run(fastP, shell=True) 

					# Calling BBDuk to quality trim the data
					preprocessing_rawdata(file, sample_name)

				print(">> {0} 4/5 | bbduk - Quality trimming and filtering of all reads: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
				print(">> {0} 5/5 | multiQC - Summing all QC reports: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
				multiQC = " ".join([
				"/home/stavros/anaconda3/bin/multiqc",  # Call MultiQC
				"--quiet",  # Print only log warnings
				"--outdir", temp,  # Create report in the FastQC reports directory
				"--filename", "{0}_init_summarised_report".format(dirs),  # Name of the output report 
				temp,  # Directory where all FastQC and Cutadapt reports reside
				"2>>", os.path.join(temp, "{0}_multiQC_init_report.txt".format(dirs))])  # Output multiQC report
				subprocess.run(multiQC, shell=True)

	os.system('mv {0}/* {1}'.format(temp, prereports_dir))
	for path, subdir, folder in os.walk(prereports_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or (name.endswith("multiQC_init_report.txt") and os.stat(file).st_size == 583):  
				os.remove(file)
	
	os.system('rm {0}/*fastp.json'.format(prereports_dir))
	os.system('rm {0}/*fastqc.zip'.format(prereports_dir))
	os.system('rm {0}/*screen.png'.format(prereports_dir))
	os.system('rm {0}/*screen.html'.format(prereports_dir))
	os.system('rm {0}/*fastQC_init_report.txt'.format(prereports_dir))
	os.system('rm -r {0}/*init_summarised_report_data'.format(prereports_dir))
	return

def preprocessing_rawdata(file, sample_name):
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whos quality is below Q20). That way we remove obvious trash and 
	trying to improve mapping. """
	bbduk = ' '.join([
	"bbduk.sh",  # Call BBDuck (BBTools) to preprocess the raw data
	"threads={0}".format(args.threads),  # Set number of threads to use
	"in={0}".format(file),  # Input of the forward file
	"out={0}".format(os.path.join(preprocesed_dir, "{0}.qt.fastq.gz".format(sample_name))),  # Export edited forward read to file
	"trimq=20",  # Regions with average quality BELOW this will be trimmed (20)
	"qtrim=r",  # Trim read ends to remove bases with Q<15
	"ordered=t",  # Keeps the reads in the same order as we gave them to the software
	"maq={0}".format(20),  # Reads with average quality (after trimming) below this will be discarded (15)
	"2>>", os.path.join(temp, "{0}_bbduk_report.txt".format(sample_name))])  # Output trimming report
	subprocess.run(bbduk, shell=True)
	return

def aligning():
	print("\n{0}ALIGNING AGAINST THE REFERENCE GENOME (GenCode.v31)".format(datetime.now().strftime("%d.%m %H:%M")))
	if not os.path.exists(alignment_dir): os.makedirs(alignment_dir)
	preprocesed_data = glob.glob("{0}/*.qt.fastq.gz".format(preprocesed_dir))
	

	print(">> {0} 1/1 | minimap2 - aligning in total {0} samples against the reference genome: in progress ..".format(datetime.now().strftime("%d.%m %H:%M"), len(preprocesed_data)))
	for i, file in enumerate(preprocesed_data, 1):
		file_name = os.path.basename(file).split(".")[0]
		minimap2_genome = " ".join([
		"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
		"-t", args.threads,  # Number of threads to use
		"-ax sr",   # Genomic short-read mapping mode and output in SAM format (-a)
		"-R", "\'@RG\\tID:{0}\\tSM:{0}\\tPL:ILLUMINA\'".format(file_name),  # SAM read group line in a format like '@RG\\tID:foo\\tSM:bar'
		# "-o", os.path.join(alignment_dir, "{0}.qt.genome.sam".format(file_name)),
		refGenome_GRCh38,  # Inputting the reference genome
		file,  # Input .fastq.gz file
		"|", "samtools view",  # Calling 'samtools view' to compress the Bowtie's output file to BAM
  		"--threads", args.threads,  # Number of threads to be used by Samtools in the conversion of the SAM files to BAM
  		"-Sb1h",  # The output file should be in BAM format, and use fast compression
  		"-",  # Piping the input file
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
  		"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
  		"-o", os.path.join(alignment_dir, "{0}.qt.genome.bam".format(file_name)), "-",  # Sorted output  BAM file
		"2>>", os.path.join(prereports_dir, "minimap2_genome_report.txt")])  # minimap2 report
		subprocess.run(minimap2_genome, shell=True)
	mapping_qc()
	return 

def mapping_qc():
	print("\n{0} MAPPING QUALITY CONTROL AND METRICS".format(datetime.now().strftime("%d.%m %H:%M")))
	if not os.path.exists(postreports_dir): os.makedirs(postreports_dir)
	aligned_data = glob.glob("{0}/*.qt.genome.bam".format(alignment_dir))


	print(">> {0} 1/3 | fastQC - Alignment Quality Control reports are being generated: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	fastQC_align = ' '.join([
	"fastqc",  # Call fastQC to quality control all processed data
	"--threads", args.threads,  # Number of threads to use
	"--quiet",  # Print only log warnings
	"--outdir", temp,  # Create all output files in this specified output directory
	' '.join(aligned_data),  # String containing all samples that are about to be checked
	"2>>", os.path.join(temp, "fastQC_alignment_report.txt")])  # Output fastQC report
	subprocess.run(fastQC_align, shell=True)

	print(">> {0} 2/3 | Picard & RSeQC - removing duplicates and generating post-alignment stats: in progress ..".format(datetime.now().strftime("%d.%m %H:%M"),))
	for i, file in enumerate(aligned_data, 1):
		file_name = os.path.basename(file).split(".")[0]
		# Removing duplicate reads
		rmDuplicates = " ".join([
		"samtools markdup",  # Call samtools markdup
		"--threads", args.threads,  # Number of threads to be used
		"--output-fmt BAM",  # Output in BAM format
		"-r",  # Remove duplicate reads
		"-s",  # Report stats
		file,  # Input BAM files to analyze
		"{0}/{1}.rdqt.genome.bam".format(alignment_dir, file_name),
		"2>>", os.path.join(temp, "{0}_rmDuplicates_report.txt".format(file_name))])
		subprocess.run(rmDuplicates, shell=True)
		
		# Indexing the output sorted deduplicated bam file
		samtools_index = " ".join([
		"samtools index",  # Indexing the concat_samples.bam file
		"-@", args.threads,  # Number of threads to be used
		"{0}/{1}.rdqt.genome.bam".format(alignment_dir, file_name),  # Input BAM file
		"2>>", os.path.join(temp, "samtools_index_report.txt")])  # Output samtools index report
		subprocess.run(samtools_index, shell=True)

		# BAM stats
		bam_stat = ' '.join([
		"bam_stat.py",
		"-i", file,  # Input BAM file
		"> {0}/{1}_bam_stat.txt".format(temp, file_name),  # Output file
		"2>>", os.path.join(temp, "{0}_bamstat_report.txt".format(file_name))])
		subprocess.run(bam_stat, shell=True)

		# BAM read distribution
		read_distribution = ' '.join([
		"read_distribution.py",  # Call samtools flagstat
		"-i", file,  # Input BAM file
		"-r", reference_annotation_bed,
		"> {0}/{1}.fragSize".format(temp, file_name),  # Output file
		"2>>", os.path.join(temp, "{0}_read_distribution_report.txt".format(file_name))])
		subprocess.run(read_distribution, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}_alignment_metrics.txt".format(temp, file_name),  # Output
		"REFERENCE_SEQUENCE= {0}".format(refGenome_GRCh38),  # Reference sequence file
		"2>>", os.path.join(temp, "{0}_CollectAlignmentSummaryMetrics_report.txt".format(file_name))])
		subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 
		os.system('rm {0}'.format(file))
	
	print(">> {0} 3/3 | multiQC - Summing all QC reports: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
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
	
	os.system('rm {0}/*fastqc.zip'.format(temp))
	os.system('rm {0}/*alignment_report.txt'.format(temp))
	os.system('rm -r {0}/*summarised_report_data'.format(temp))
	os.system('mv {0}/* {1}'.format(temp, postreports_dir))
	return

def post_processing():
	# time.sleep(10)
	print("\n{0} POST-PROCESSING AND PREPARING DATA FOR VARIANT CALLING".format(datetime.now().strftime("%d.%m %H:%M")))
	if not os.path.exists(postprocesed_dir): os.makedirs(postprocesed_dir)
	data = glob.glob("{0}/*.rdqt.genome.bam".format(alignment_dir))
	list_of_aligned = ' '.join(["--input_file " + files for files in data])
	known_sites = [f for f in glob.glob("{0}/*hg38.vcf".format(gatk_known_sites))]

	if not os.path.exists("{0}.dict".format(refGenome_GRCh38[:-3])):
		print("gatk - Building the reference dictionary for GATK: in progress ..")
		ref_dict = " ".join([
		"gatk CreateSequenceDictionary",  # Calling gatk CreateSequenceDictionary
		"--REFERENCE", refGenome_GRCh38,  # Input reference genome
		"--OUTPUT", "{0}.dict".format(refGenome_GRCh38[:-3]),  # Output file
		"2>>", os.path.join(temp, "ref_dict_report.txt")])  # Output report
		subprocess.run(ref_dict, shell=True)
	
	if not os.path.exists("{0}.fai".format(refGenome_GRCh38)):
		print("samtools faidx - Indexing the reference genome: in progress ..")
		ref_idx = " ".join([
		"samtools faidx",  # Calling samtools faidx
		refGenome_GRCh38,  # Input reference genome
		"2>>", os.path.join(temp, "samtools_idx_report.txt")])  # Output report
		subprocess.run(ref_idx, shell=True)
	
	print(">> {0} 1/6 | gatk RealignerTargetCreator - Define intervals to target for local realignment: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	known_indels = [f for f in glob.glob("{0}/*indels*.vcf".format(gatk_known_sites))]
	## IDefine intervals to target for local realignment
	realignerTargetCreator = " ".join([
	gatk3_8, "--analysis_type RealignerTargetCreator",  # Call gatk RealignerTargetCreator (v3.8.1.0)
	list_of_aligned,  # List of all aligned BAM files
	"--reference_sequence", refGenome_GRCh38,
	"--known", known_indels[0],  # databases of known InDels
	"--known", known_indels[1],  # databases of known InDels
	"--num_threads", args.threads,  # Number of data threads to allocate to this analysis
	"--out", "{0}/indelRealigner.intervals".format(postprocesed_dir),
	"2>>", os.path.join(postreports_dir, "realignerTargetCreator_report.txt")])
	subprocess.run(realignerTargetCreator, shell=True)

	
	print(">> {0} 2/6 | gatk IndelRealigner - Base Quality Score Recalibration: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	# InDel realignemnt
	indelRealigner = " ".join([
	gatk3_8, "--analysis_type IndelRealigner",  # 
	list_of_aligned,  # List of all aligned BAM files
	"--reference_sequence", refGenome_GRCh38,
	"--targetIntervals", "{0}/indelRealigner.intervals".format(postprocesed_dir),
	"--knownAlleles", known_sites[0],  # databases of known InDels
	"--knownAlleles", known_sites[1],  # databases of known InDels
	"--knownAlleles", known_sites[2],  # databases of known InDels
	"--nWayOut", "\'.realigned.bam\'",  # list of output files
	"--noOriginalAlignmentTags",  # Don't output the original cigar for each realigned read
	"2>>", os.path.join(postreports_dir, "indelRealigner_report.txt")])
	subprocess.run(indelRealigner, shell=True, cwd=postprocesed_dir)


	real_data = glob.glob("{0}/*rdqt.genome.realigned.bam".format(postprocesed_dir))
	for i, file in enumerate(real_data, 1):
		file_name = os.path.basename(file).split(".")[0]
		# Phred score recalibration
		if i==1:print(">> {0} 3/6 | gatk BaseRecalibrator 1 - Base Quality Score Recalibration: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
		preBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk BaseRecalibrator (v4.1.2.0)
		"--input", file,
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--reference", refGenome_GRCh38,
		"--intervals", gatk_wg_intervals,
		"--output", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),
		"2>>", os.path.join(postreports_dir, "preBaseRecalibrator_report.txt")])
		subprocess.run(preBaseRecalibrator, shell=True)

		# Generation of the recalibrated reads
		if i==1:print(">> {0} 4/6 | gatk applyBQSR - Generation of the recalibrated reads: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
		applyBQSR = " ".join([
		"gatk ApplyBQSR",  # Call gatk ApplyBQSR (v4.1.2.0)
		"--intervals", gatk_wg_intervals,
		"--input", file,  # Input file containing sequence data
		"--bqsr-recal-file", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),  # File with base quality score recalibration
		"--output", "{0}/{1}.real.recal.bam".format(postprocesed_dir, file_name),  # Output file
		"2>>", os.path.join(postreports_dir, "applyBQSR_report.txt")])
		subprocess.run(applyBQSR, shell=True)

		# Phred score recalibration
		if i==1:print(">> {0} 5/6 | gatk BaseRecalibrator 2 - Base Quality Score Recalibration: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
		afterBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk BaseRecalibrator (v4.1.2.0)
		"--reference", refGenome_GRCh38,
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--intervals", gatk_wg_intervals,
		"--input", "{0}/{1}.real.recal.bam".format(postprocesed_dir, file_name),
		"--output", "{0}/{1}_post_recal_data.table".format(postprocesed_dir, file_name),
		"2>>", os.path.join(postreports_dir, "afterBaseRecalibrator_report.txt")])
		subprocess.run(afterBaseRecalibrator, shell=True)

		# Evaluate and compare base quality score recalibration (BQSR) tables
		if i==1:print(">> {0} 6/6 | gatk AnalyzeCovariates - Evaluate and compare base quality score recalibration: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
		analyzeCovariates = " ".join([
		"gatk AnalyzeCovariates",  # Call gatk ApplyBQSR (v4.1.2.0)
		"--before-report-file", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),
		"--after-report-file", "{0}/{1}_post_recal_data.table".format(postprocesed_dir, file_name),
		"--intermediate-csv-file", "{0}/{1}.BQSR.csv".format(postreports_dir, file_name),	
		"--plots-report-file", "{0}/{1}BQSR.pdf".format(postreports_dir, file_name),	
		"2>>", os.path.join(postreports_dir, "analyzeCovariates_report.txt")])
		subprocess.run(analyzeCovariates, shell=True)
	os.system('rm {0}/*.rdqt.genome.realigned.*'.format(postprocesed_dir))
	return

def calling_variants():
	print("\n{0} CALLING VARIANTS".format(datetime.now().strftime("%d.%m %H:%M")))

	if not os.path.exists(varcall_dir): os.makedirs(varcall_dir)
	calibrated_data = glob.glob("{0}/*.real.recal.bam".format(postprocesed_dir))
	sample_list = os.path.join(varcall_dir, "sample_list.txt")
	
	# Creating sample list that is needed from VarScan2
	# if not os.path.exists(sample_list):
	with open(sample_list, "w") as fout: 
		for items in calibrated_data:
			fout.write("{0}\n".format(os.path.basename(items).split(".")[0]))


	### CALL CANDIDATE VARIANTS
	## 1. Run freebayes to call variants
	print(">> {0} 1/4 | freebayes - Call somatic SNVs and indels: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	freebayes = " ".join([
	"freebayes",
	"--bam", ' '.join(calibrated_data),
	"--fasta-reference", refGenome_GRCh38,
	"|", "bgzip",
	"--threads", args.threads, 
	"> {0}/variants_freebayes.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "freebayes_report.txt")])
	subprocess.run(freebayes, shell=True)


	## 2. Run VarScan2 to call variants
	print(">> {0} 2/4 | VarScan2 - Call somatic SNVs and indels: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	varscan2 = " ".join([
	"samtools mpileup",  # Converting BAM input to vcf
	"--redo-BAQ",  # Recalculate BAQ on the fly
	"--fasta-ref", refGenome_GRCh38,  # Indexed reference sequence file
	' '.join(calibrated_data),
	"| varscan mpileup2cns",  # Call consensus & variants from the pileup file
	"--vcf-sample-list", sample_list,
	"--p-value 0.01", # Using the default p-value
	"--output-vcf",  # Outputs in VCF format
	"--variants",  # Report only variant positions
	"|", "bgzip",  # GZ the output
	"--threads", args.threads, 
	"> {0}/variants_varscan2.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "varscan2_report.txt")])
	subprocess.run(varscan2, shell=True)


	## 3. Run Mutect2 to call variants
	print(">> {0} 3/4 | Mutect2 - Call somatic SNVs and indels via local assembly of haplotypes: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	list_of_aligned = ' '.join(["--input " + files for files in calibrated_data])
	mutect2 = " ".join([
	"gatk Mutect2",
	list_of_aligned,
	"--panel-of-normals", gatk_pon,
	"--reference", refGenome_GRCh38,
	"--intervals", gatk_wg_intervals,
	"--germline-resource", gatk_liftover_GRCh38,
	"--native-pair-hmm-threads", args.threads,
	"--output", "{0}/variants_mutect2.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "mutect2_report.txt")])
	subprocess.run(mutect2, shell=True)

	called_variants = glob.glob("{0}/*.vcf.gz".format(varcall_dir))
	for file in called_variants:
		if not os.path.exists("{0}.tbi".format(file)):
			subprocess.run("tabix -p vcf {0}".format(file), shell=True)


	## 4. Create intersections, unions and complements of the VCF files
	if not os.path.exists(intersection_analysis): os.makedirs(intersection_analysis)
	print(">> {0} 4/4 | StarFish - Generating the intersections of the VCF files: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	intersecvcfs = " ".cjoin([
	starfish,  # Calling starfish
	"--rtg rtg",
	"--all_records",  # Intersect all records
	"--threads", args.threads,
	"--sdf", refGenome_GRCh38sdf,
	"--output", intersection_analysis,
	"--variants", ' '.join(called_variants),
	"--vennout", "{0}/venn_diagram.pdf".format(intersection_analysis),
	"--sample", os.path.basename(calibrated_data[0]).split(".")[0],
	"--names", ' '.join([os.path.basename(f).split("_")[1].split(".")[0] for f in called_variants]),
	"2>>", os.path.join(postreports_dir, "starfisg_vcfs-intersection-report.txt")])
	subprocess.run(intersecvcfs, shell=True)
	return

def stats_n_annotation():
	print("\n{0} EXPORTING STATS, FILTERING & ANNOTATING THE VARIANTS".format(datetime.now().strftime("%d.%m %H:%M")))
	if not os.path.exists(prefilter_plots): os.makedirs(prefilter_plots)
	if not os.path.exists(filtered_plots): os.makedirs(filtered_plots)
	
	### STATS AND FILTERING STEPS
	## 1. Exporting stats of the called variants
	print(">> {0} 1/8 | rtg stats - Exporting stats of the called variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	rtg_stats = " ".join([
	"rtg vcfstats",
	"{0}/ABC.vcf.gz".format(intersection_analysis),
	">", "{0}/rtg_stats.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "rtg_stats_report.txt")])
	subprocess.run(rtg_stats, shell=True)

	print(">> {0} 2/8 | Bcftools stats - Exporting stats of the called variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	bcftools_stats = " ".join([
	"bcftools stats",
	"--samples -",  # list of samples for sample stats (include all)
	"--fasta-ref", refGenome_GRCh38,  # reference to determine INDEL context
	"{0}/ABC.vcf.gz".format(intersection_analysis),  # input file
	">", "{0}/bcftools_stats.txt".format(varcall_dir),  # output file
	"2>>", os.path.join(postreports_dir, "bcftools_stats_report.txt")])
	subprocess.run(bcftools_stats, shell=True)

	print(">> {0} 3/8 | Plot vcfstats - Exporting plots based on the stats of the called variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	plot_vcfstats = " ".join([
	"plot-vcfstats",
	"--prefix", prefilter_plots,
	"{0}/bcftools_stats.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "plot_vcfstats_report.txt")])
	subprocess.run(plot_vcfstats, shell=True)
	
	## 2. Filtering variants
	print(">> {0} 4/8 | Bcftools filter - Applying several filters to the called variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	bcftools_filter = " ".join([
	"bcftools filter",
	"--output-type z",
	"--threads", args.threads,
	"-e", "\'%QUAL<=30\'",  # %DP<=10 || %QUAL/INFO/AO<=2 || ||  SAF<=2 || SAR<=2
	"{0}/ABC.vcf.gz".format(intersection_analysis),
	">", "{0}/filtered_variants.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "bcftools_stats_report.txt")])
	subprocess.run(bcftools_filter, shell=True)

	# 3. Indexing the filtered vcf file
	subprocess.run("tabix -p vcf {0}/filtered_variants.vcf.gz".format(varcall_dir), shell=True) 

	## 4. Exporting stats of the filtered variants
	print(">> {0} 5/8 | rtg stats - Exporting stats of the filtered variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	filt_rtg_stats = " ".join([
	"rtg vcfstats",
	"{0}/filtered_variants.vcf.gz".format(varcall_dir),
	">", "{0}/rtg_stats.filtered.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "rtg_stats_filtered_report.txt")])
	subprocess.run(filt_rtg_stats, shell=True)

	print(">> {0} 6/8 | Bcftools stats - Exporting stats of the filtered variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	filt_bcftools_stats = " ".join([
	"bcftools stats",
	"--samples -",  # list of samples for sample stats (include all)
	"--fasta-ref", refGenome_GRCh38,  # reference to determine INDEL context
	"{0}/filtered_variants.vcf.gz".format(varcall_dir),  # input file
	">", "{0}/bcftools_stats.filtered.txt".format(varcall_dir),  # output file
	"2>>", os.path.join(postreports_dir, "bcftools_stats_filtered_report.txt")])
	subprocess.run(filt_bcftools_stats, shell=True)

	print(">> {0} 7/8 | Plot vcfstats - Exporting plots based on the stats of the filtered variants: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	filt_plot_vcfstats = " ".join([
	"plot-vcfstats",
	"--prefix", filtered_plots,
	"{0}/bcftools_stats.filtered.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "plot_vcfstats_filtered_report.txt")])
	subprocess.run(filt_plot_vcfstats, shell=True)

	### ANNOTATE VARIANTS
	## 5. Adding information to the discovered variants
	## At this step we run tools to add information to the discovered variants in our dataset.
	## One of those tools, Funcotator, can be used to add gene-level information to each variant.
	print(">> {0} 8/8 | Funcotator - Adding information to the discovered variants in our dataset: in progress ..".format(datetime.now().strftime("%d.%m %H:%M")))
	funcotator = " ".join([
	"gatk Funcotator",
	"--ref-version hg38",
	"--reference", refGenome_GRCh38,
	"--variant", "{0}/filtered_variants.vcf.gz".format(varcall_dir),
	"--output-file-format VCF",  # The output file format
	"--output", "{0}/filtered_variants.annot.vcf".format(varcall_dir),
	"--data-sources-path", gatk_functional_annot,
	"2>>", os.path.join(postreports_dir, "functional_annotator_report.txt")])
	subprocess.run(funcotator, shell=True)
	return

def main():
	

	# quality_control()
	
	# aligning()

	# post_processing()

	calling_variants()

	stats_n_annotation()


	## CHECK FOR "WARN" IN EVERY REPORT file! and print

if __name__ == "__main__": main()