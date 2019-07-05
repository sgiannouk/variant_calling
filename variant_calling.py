###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
import shutil, time, glob, sys, os, re


inputFolder = "/shared/projects/jill_dnapanel_pangaea"
fastQscreen_config = "/home/stavros/playground/progs/16S_subsidiary_files/fastq_screen.conf"
refGenome_GRCh38 = "/home/stavros/playground/progs/reference_files/reference_genome/GRCh38_primAssembly/GRCh38_primary_assembly_genome.fa"
reference_annotation_gtf = "/home/stavros/playground/progs/reference_files/gene_annotation/gencode.v29.primary_assembly.annotation.gtf"
reference_annotation_bed = "/home/stavros/playground/progs/reference_files/gene_annotation/hg38_Gencode_V28.bed"
gatk_knownsites = "/home/stavros/playground/progs/reference_files/gatk_knownsites"


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

preprocesed_dir = os.path.join(analysis_dir, "preprocesed_data")
prereports_dir = os.path.join(analysis_dir, "preprocessing_reports")

alignment_dir = os.path.join(analysis_dir, "alignments")

postprocesed_dir = os.path.join(analysis_dir, "postprocessed_data")
postreports_dir = os.path.join(analysis_dir, "postprocessing_reports")

varcall_dir = os.path.join(analysis_dir, "variant_calling")



def quality_control():

	if not os.path.exists(prereports_dir): os.makedirs(prereports_dir)
	if not os.path.exists(preprocesed_dir): os.makedirs(preprocesed_dir)
	if not os.path.exists(temp): os.makedirs(temp)
	print("PREPROCESSING & QUALITY CONTROL")
	for path, subdir, folder in os.walk(args.project_dir):
		for i, dirs in enumerate(subdir, 0):
			if dirs == "Tanda_106":
				print("{0}/{1} | Analysing {2}".format(i, len(subdir), dirs))
				mfiltered_data = [f for f in glob.glob(os.path.join(os.path.join(path, dirs), "*.fastq.gz"))]

				print("fastQscreen - Checking random reads for possible contamination: in progress ..")
				fastQscreen = ' '.join([
				"fastq_screen",  # Call fastQ screen to check contamination in the processed data
				"--threads", args.threads,  # Number of threads to use
				"--outdir",  temp,  # Directory in which the output files will be saved
				"--quiet",  # Suppress all progress reports on stderr and only report errors
				"--conf", fastQscreen_config,  # Location of the required configuration file
				' '.join(mfiltered_data),
				"2>>", os.path.join(temp, "{0}_fastQscreen_report.txt".format(dirs))])  # Output fastQ screen report
				subprocess.run(fastQscreen, shell=True)
			
				print("fastQC - Quality Control reports are being generated: in progress ..")
				fastQC = ' '.join([
				"fastqc",  # Call fastQC to quality control all processed data
				"--threads", args.threads,  # Number of threads to use
				"--quiet",  # Print only log warnings
				"--outdir", temp,  # Create all output files in this specified output directory
				' '.join(mfiltered_data),  # String containing all samples that are about to be checked
				"2>>", os.path.join(temp, "{0}_fastQC_init_report.txt".format(dirs))])  # Output fastQC report
				subprocess.run(fastQC, shell=True)

				print("fastP - Quality Control of all reads: in progress ..")
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

				print("multiQC - Summing all QC reports: in progress ..")
				print("bbduk - Quality trimming and filtering of all reads: in progress ..")
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
	if not os.path.exists(alignment_dir): os.makedirs(alignment_dir)
	if not os.path.exists(postprocesed_dir): os.makedirs(postprocesed_dir)
	preprocesed_data = glob.glob("{0}/*.qt.fastq.gz".format(preprocesed_dir))
	print("\nALIGNING AND POST-ALIGNMENT QUALITY CONTROL")
	print("minimap2 - aligning in total {0} samples against the reference genome: in progress ..".format(len(preprocesed_data)))
	for i, file in enumerate(preprocesed_data):
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

	
	# Post Alignment QC
	if not os.path.exists(postreports_dir): os.makedirs(postreports_dir)
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

	print("Picard & RSeQC - removing duplicates and generating post-alignment stats: in progress ..")
	for i, file in enumerate(aligned_data):
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
		# subprocess.run(read_distribution, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}_alignment_metrics.txt".format(temp, file_name),  # Output
		"REFERENCE_SEQUENCE= {0}".format(refGenome_GRCh38),  # Reference sequence file
		"2>>", os.path.join(temp, "{0}_CollectAlignmentSummaryMetrics_report.txt".format(file_name))])
		subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 
		# os.system('rm {0}'.format(file))
	
	print("multiQC - Summing all QC reports: in progress ..")
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", temp,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	temp,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(temp, "post-alignment_multiQC_report.txt")])  # Output multiQC report
	# subprocess.run(multiQC, shell=True)
	

	# for path, subdir, folder in os.walk(temp):
	# 	for name in folder:
	# 		file = os.path.join(path, name)
	# 		if os.stat(file).st_size == 0 or (name.endswith("multiQC_report.txt") and os.stat(file).st_size == 583):  
	# 			os.remove(file)
	
	# os.system('rm {0}/*fastqc.zip'.format(temp))
	# os.system('rm {0}/*alignment_report.txt'.format(temp))
	# os.system('rm -r {0}/*summarised_report_data'.format(temp))
	os.system('mv {0}/* {1}'.format(temp, postreports_dir))
	return

def post_processing():
	# time.sleep(10)
	print("\nPOST-PROCESSING AND PREPARING DATA FOR VARIANT CALLING")
	data = glob.glob("{0}/*.rdqt.genome.bam".format(alignment_dir))
	
	# Merging all bam files into one (keeping group info)
	print("Samtools - Merge all BAM files to one: in progress ..")
	samtools_merge = " ".join([
	"samtools merge",  # Call samtools merge
	"--threads", args.threads,  # Number of threads to be used
	"{0}/concat_samples.bam".format(postprocesed_dir),
	' '.join(data),  # Input BAM files
	"2>>", os.path.join(temp, "samtools_merge_report.txt")])  # Output multiQC report
	subprocess.run(samtools_merge, shell=True)

	# print("Samtools index - Indexing the concatenated BAM file: in progress ..")
	samtools_index = " ".join([
	"samtools index",  # Indexing the concat_samples.bam file
	"-@", args.threads,  # Number of threads to be used
	"{0}/concat_samples.bam".format(postprocesed_dir),  # Input BAM file
	"2>>", os.path.join(temp, "samtools_index_report.txt")])  # Output samtools index report
	subprocess.run(samtools_index, shell=True)

	if not os.path.exists("{0}.dict".format(refGenome_GRCh38[:-3])):
		print("gatk - Building the reference dictionary for GATK: in progress ..")
		ref_dict = " ".join([
		"gatk CreateSequenceDictionary",  # Calling gatk CreateSequenceDictionary
		"--REFERENCE", refGenome_GRCh38,  # Input reference genome
		"--OUTPUT", "{0}.dict".format(refGenome_GRCh38[:-3]),  # Output file
		"2>>", os.path.join(temp, "ref_dict_report.txt")])  # Output report
		subprocess.run(ref_dict, shell=True)

	# Phred score recalibration
	known_sites = [f for f in glob.glob("{0}/*.vcf.gz".format(gatk_knownsites))]
	print("gatk BaseRecalibrator 1 - Base Quality Score Recalibration: in progress ..")
	preBaseRecalibrator = " ".join([
	"gatk BaseRecalibrator",  # Call gatk MarkDuplicates(v4.1.2.0)
	"--input", "{0}/concat_samples.bam".format(postprocesed_dir),
	"--reference", refGenome_GRCh38,
	"--output", "{0}/pre_recal_data.table".format(postprocesed_dir),
	"--known-sites", known_sites[0],  # databases of known polymorphic sites
	"--known-sites", known_sites[1],  # databases of known polymorphic sites
	"--known-sites", known_sites[2],  # databases of known polymorphic sites
	"2>>", os.path.join(postreports_dir, "preBaseRecalibrator_report.txt")])
	subprocess.run(preBaseRecalibrator, shell=True)

	# Generation of the recalibrated reads
	print("gatk applyBQSR - Generation of the recalibrated reads: in progress ..")
	applyBQSR = " ".join([
	"gatk ApplyBQSR",  # Call gatk ApplyBQSR (v4.1.2.0)
	"--input", "{0}/concat_samples.bam".format(postprocesed_dir),
	"--bqsr-recal-file", "{0}/pre_recal_data.table".format(postprocesed_dir),
	"--output", "{0}/recalibrated_concat_samples.bam".format(postprocesed_dir),	
	"2>>", os.path.join(postreports_dir, "applyBQSR_report.txt")])
	subprocess.run(applyBQSR, shell=True)

	# Phred score recalibration
	print("gatk BaseRecalibrator 2 - Base Quality Score Recalibration: in progress ..")
	afterBaseRecalibrator = " ".join([
	"gatk BaseRecalibrator",  # Call gatk MarkDuplicates(v4.1.2.0)
	"--input", "{0}/recalibrated_concat_samples.bam".format(postprocesed_dir),
	"--reference", refGenome_GRCh38,
	"--output", "{0}/after_recal_data.table".format(postprocesed_dir),
	"--known-sites", known_sites[0],  # databases of known polymorphic sites
	"--known-sites", known_sites[1],  # databases of known polymorphic sites
	"--known-sites", known_sites[2],  # databases of known polymorphic sites
	"2>>", os.path.join(postreports_dir, "afterBaseRecalibrator_report.txt")])
	subprocess.run(afterBaseRecalibrator, shell=True)

	# Evaluate and compare base quality score recalibration (BQSR) tables
	print("gatk analyzeCovariates - Evaluate and compare base quality score recalibration: in progress ..")
	analyzeCovariates = " ".join([
	"gatk AnalyzeCovariates",  # Call gatk ApplyBQSR (v4.1.2.0)
	"--before-report-file", "{0}/pre_recal_data.table".format(postprocesed_dir),
	"--after-report-file", "{0}/after_recal_data.table".format(postprocesed_dir),
	"--intermediate-csv-file", "{0}/concat_recal_samples_BQSR.csv".format(postreports_dir),	
	"--plots-report-file", "{0}/concat_recal_samples_BQSR.pdf".format(postreports_dir),	
	"2>>", os.path.join(postreports_dir, "analyzeCovariates_report.txt")])
	subprocess.run(analyzeCovariates, shell=True)

	os.system('rm {0}/concat_samples*'.format(postreports_dir))
	os.system('rm {0}/*data.table'.format(postreports_dir))
	return

def variant_calling():
	print("CALLING VARIANTS")
	# if not os.path.exists(varcall_dir): os.makedirs(varcall_dir)
	# aligned_data = glob.glob("{0}/*.qt.genome.bam".format(alignment_dir))
	
	
	# for i, file in enumerate(aligned_data):
	# 	file_name = os.path.basename(file).split(".")[0]
	# 	print("{0}/{1} | Variant Calling {2} using GATK4...".format(i, len(aligned_data), file_name))
		

	## 1. Generate OXOG metrics:
	# "gatk CollectSequencingArtifactMetrics"
	# -I Tumor_Sample_Alignment.bam \
	# -O <job_identifier> \
	# --FILE_EXTENSION .txt \
	# -R GRCh38.d1.vd1.fa  ## Only chr1-22 + XYM

	
	# ## 2. Generate pileup summaries on tumor sample:
	# "gatk GetPileupSummaries"
	# -I Tumor_Sample_Alignment.bam \
	# -O <job_identifier>.targeted_sequencing.table \
	# -V af-only-gnomad-common-biallelic.grch38.main.vcf.gz \ # Germline reference from gnomad
	# -L intervals.bed \ ## Only chr1-22 + XYM
	# -R GRCh38.d1.vd1.fa

	# ## 3. Calculate contamination on tumor sample
	# "gatk CalculateContamination" \
	# -I <job_identifier>.targeted_sequencing.table \ # From step 2
	# -O <job_identifier>.targeted_sequencing.contamination.table

	# ## 4. Find tumor sample name from BAM
	# "gatk GetSampleName" \
	# -I Tumor_Sample_Alignment.bam \
	# -O <job_identifier>.targeted_sequencing.sample_name

	# ## 5. Run MuTect2 using only tumor sample on chromosome level (25 commands with different intervals)
	# "gatk Mutect2" \
	# -R GRCh38.d1.vd1.fa \
	# -L chr4:1-190214555 \ # Specify chromosome
	# -I Tumor_Sample_Alignment.bam \
	# -O 3.mt2.vcf \
	# -tumor <tumor_sample_name> \ # From step 4
	# --af-of-alleles-not-in-resource 2.5e-06 \
	# --germline-resource af-only-gnomad.hg38.vcf.gz \ # Germline reference from gnomad
	# -pon gatk4_mutect2_4136_pon.vcf.gz # New panel of normal created by 4136 TCGA curated normal samples, using GATK4


	# ## After this step, all chromosome level VCFs are merged into one.

	# ## 6. Sort VCF with Picard
	# "gatk SortVcf" \
	# SEQUENCE_DICTIONARY=GRCh38.d1.vd1.dict \
	# OUTPUT=<job_identifier>.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz \
	# I=merged_multi_gatk4_mutect2_tumor_only_calling.vcf \ # From step 5
	# CREATE_INDEX=true

	# ## 7. Filter variant calls from MuTect
	# "gatk FilterMutectCalls" \
	# -O <job_identifier>.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \
	# -V <job_identifier>.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz \ # From step 6
	# --contamination-table <job_identifier>.targeted_sequencing.contamination.table \ # From step 3
	# -L intervals.bed

	# ## 8. Filter variants by orientation bias
	# "gatk FilterByOrientationBias" \
	# -O <job_identifier>.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.vcf.gz \ # final output
	# -P <job_identifier>.pre_adapter_detail_metrics.txt \ # From step 1
	# -V <job_identifier>.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \ # From step 7
	# -L intervals.bed \
	# -R GRCh38.d1.vd1.fa \
	# -AM G/T \
	# -AM C/T

	return

def main():
	

	# quality_control()
	
	# aligning()

	# post_processing()

	variant_calling()

if __name__ == "__main__": main()