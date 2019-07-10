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
whole_genome_intervals = "/home/stavros/playground/progs/reference_files/gatk_subfiles/intervals/wgs_calling_regions.hg38.interval_list"
gatk_subfiles = "/home/stavros/playground/progs/reference_files/gatk_subfiles"
gatk_functional_annot = "/home/stavros/playground/progs/reference_files/gatk_subfiles/functional_annot/funcotator_dataSources.v1.6.20190124s"
gatk3_8 ="gatk --java-options \"-Xmx8G\" /home/stavros/playground/progs/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar"


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
parser.add_argument('-th', '--threads', dest='threads', default=str(40), metavar='', 
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
prefilter_plots = os.path.join(varcall_dir, "prefilter_plots")
filtered_plots = os.path.join(varcall_dir, "filtered_plots")


def quality_control():
	if not os.path.exists(prereports_dir): os.makedirs(prereports_dir)
	if not os.path.exists(preprocesed_dir): os.makedirs(preprocesed_dir)
	if not os.path.exists(temp): os.makedirs(temp)
	print("PREPROCESSING & QUALITY CONTROL")
	for path, subdir, folder in os.walk(args.project_dir):
		for i, dirs in enumerate(subdir, 0):
			# if dirs == "Tanda_104":
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
	
	print("multiQC - Summing all QC reports: in progress ..")
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
	print("\nPOST-PROCESSING AND PREPARING DATA FOR VARIANT CALLING")
	data = glob.glob("{0}/*.rdqt.genome.bam".format(alignment_dir))

	if not os.path.exists("{0}.dict".format(refGenome_GRCh38[:-3])):
		print("gatk - Building the reference dictionary for GATK: in progress ..")
		ref_dict = " ".join([
		"gatk CreateSequenceDictionary",  # Calling gatk CreateSequenceDictionary
		"--REFERENCE", refGenome_GRCh38,  # Input reference genome
		"--OUTPUT", "{0}.dict".format(refGenome_GRCh38[:-3]),  # Output file
		"2>>", os.path.join(temp, "ref_dict_report.txt")])  # Output report
		subprocess.run(ref_dict, shell=True)

	print("1/4 | gatk BaseRecalibrator 1 - Base Quality Score Recalibration: in progress ..")
	print("1/4 | gatk BaseRecalibrator 1 - Base Quality Score Recalibration: in progress ..")
	print("2/4 | gatk applyBQSR - Generation of the recalibrated reads: in progress ..")
	print("3/4 | gatk BaseRecalibrator 2 - Base Quality Score Recalibration: in progress ..")
	print("4/4 | gatk analyzeCovariates - Evaluate and compare base quality score recalibration: in progress ..")
	known_sites = [f for f in glob.glob("{0}/known_sites/*.vcf.gz".format(gatk_subfiles))]
	for file in data:
		file_name = os.path.basename(file).split(".")[0]

		if not os.path.exists("{0}.bai".format(file)):
			samtools_index = " ".join([
			"samtools index",  # Indexing the concat_samples.bam file
			"-@", args.threads,  # Number of threads to be used
			file,  # Input BAM file
			"2>>", os.path.join(temp, "samtools_index_report.txt")])  # Output samtools index report
			subprocess.run(samtools_index, shell=True)
	
		# Phred score recalibration
		preBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk MarkDuplicates(v4.1.2.0)
		"--input", file,
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--reference", refGenome_GRCh38,
		"--intervals", whole_genome_intervals,
		"--output", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),
		"2>>", os.path.join(postreports_dir, "preBaseRecalibrator_report.txt")])
		subprocess.run(preBaseRecalibrator, shell=True)

		# Generation of the recalibrated reads
		applyBQSR = " ".join([
		"gatk ApplyBQSR",  # Call gatk ApplyBQSR (v4.1.2.0)
		"--input", file,
		"--intervals", whole_genome_intervals,
		"--output", "{0}/{1}.recalibrated.bam".format(postprocesed_dir, file_name),	
		"--bqsr-recal-file", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),
		"2>>", os.path.join(postreports_dir, "applyBQSR_report.txt")])
		subprocess.run(applyBQSR, shell=True)

		# Phred score recalibration
		afterBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk MarkDuplicates(v4.1.2.0)
		"--reference", refGenome_GRCh38,
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--intervals", whole_genome_intervals,
		"--input", "{0}/{1}.recalibrated.bam".format(postprocesed_dir, file_name),
		"--output", "{0}/{1}_post_recal_data.table".format(postprocesed_dir, file_name),
		"2>>", os.path.join(postreports_dir, "afterBaseRecalibrator_report.txt")])
		subprocess.run(afterBaseRecalibrator, shell=True)

		# Evaluate and compare base quality score recalibration (BQSR) tables
		analyzeCovariates = " ".join([
		"gatk AnalyzeCovariates",  # Call gatk ApplyBQSR (v4.1.2.0)
		"--before-report-file", "{0}/{1}_recal_data.table".format(postprocesed_dir, file_name),
		"--after-report-file", "{0}/{1}_post_recal_data.table".format(postprocesed_dir, file_name),
		"--intermediate-csv-file", "{0}/{1}.BQSR.csv".format(postreports_dir, file_name),	
		"--plots-report-file", "{0}/{1}BQSR.pdf".format(postreports_dir, file_name),	
		"2>>", os.path.join(postreports_dir, "analyzeCovariates_report.txt")])
		subprocess.run(analyzeCovariates, shell=True)

	os.system('rm {0}/concat_samples*'.format(postprocesed_dir))
	os.system('rm {0}/*data.table'.format(postprocesed_dir))
	return

def calling_variants():
	print("CALLING VARIANTS")
	if not os.path.exists(prefilter_plots): os.makedirs(prefilter_plots)
	if not os.path.exists(filtered_plots): os.makedirs(filtered_plots)
	calibrated_data = glob.glob("{0}/*.recalibrated.bam".format(postprocesed_dir))
	
	### CALL CANDIDATE VARIANTS
	## 1. Run freebayes to call variants
	print("1/9 | freebayes - Call somatic SNVs and indels: in progress ..")
	freebayes = " ".join([
	"freebayes",
	"--bam", ' '.join(calibrated_data),
	"--fasta-reference", refGenome_GRCh38,
	"|", "gzip >,", "{0}/variants.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "freebayes_report.txt")])
	subprocess.run(freebayes, shell=True)

	# 2. Indexing the vcf file
	subprocess.run("tabix -p vcf {0}/variants.vcf.gz".format(varcall_dir), shell=True) 

	### STATS AND FILTERING STEPS
	## 3. Exporting stats of the called variants
	print("2/9 | rtg stats - Exporting stats of the called variants: in progress ..")
	rtg_stats = " ".join([
	"rtg vcfstats",
	"{0}/variants.vcf.gz".format(varcall_dir),
	">", "{0}/rtg_stats.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "rtg_stats_report.txt")])
	subprocess.run(rtg_stats, shell=True)

	print("3/9 | Bcftools stats - Exporting stats of the called variants: in progress ..")
	bcftools_stats = " ".join([
	"bcftools stats",
	"--samples -",  # list of samples for sample stats (include all)
	"--fasta-ref", refGenome_GRCh38,  # reference to determine INDEL context
	"{0}/variants.vcf.gz".format(varcall_dir),  # input file
	">", "{0}/bcftools_stats.txt".format(varcall_dir),  # output file
	"2>>", os.path.join(postreports_dir, "bcftools_stats_report.txt")])
	subprocess.run(bcftools_stats, shell=True)

	print("4/9 | Plot vcfstats - Exporting plots based on the stats of the called variants: in progress ..")
	plot_vcfstats = " ".join([
	"plot-vcfstats",
	"--prefix", prefilter_plots,
	"{0}/bcftools_stats.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "plot_vcfstats_report.txt")])
	subprocess.run(plot_vcfstats, shell=True)
	
	## 4. Filtering variants
	print("5/9 | Bcftools filter - Applying several filters to the called variants: in progress ..")
	bcftools_filter = " ".join([
	"bcftools filter",
	"--output-type z",
	"--threads", args.threads,
	"-e", "\'%QUAL<=30 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2\'",  # %DP<=10 ||
	"{0}/variants.vcf.gz".format(varcall_dir),
	">", "{0}/filtered_variants.vcf.gz".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "bcftools_stats_report.txt")])
	subprocess.run(bcftools_filter, shell=True)

	# 5. Indexing the filtered vcf file
	subprocess.run("tabix -p vcf {0}/filtered_variants.vcf.gz".format(varcall_dir), shell=True) 

	## 6. Exporting stats of the filtered variants
	print("6/9 | rtg stats - Exporting stats of the filtered variants: in progress ..")
	filt_rtg_stats = " ".join([
	"rtg vcfstats",
	"{0}/filtered_variants.vcf.gz".format(varcall_dir),
	">", "{0}/rtg_stats.filtered.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "rtg_stats_filtered_report.txt")])
	subprocess.run(filt_rtg_stats, shell=True)

	print("7/9 | Bcftools stats - Exporting stats of the filtered variants: in progress ..")
	filt_bcftools_stats = " ".join([
	"bcftools stats",
	"--samples -",  # list of samples for sample stats (include all)
	"--fasta-ref", refGenome_GRCh38,  # reference to determine INDEL context
	"{0}/filtered_variants.vcf.gz".format(varcall_dir),  # input file
	">", "{0}/bcftools_stats.filtered.txt".format(varcall_dir),  # output file
	"2>>", os.path.join(postreports_dir, "bcftools_stats_filtered_report.txt")])
	subprocess.run(filt_bcftools_stats, shell=True)

	print("8/9 | Plot vcfstats - Exporting plots based on the stats of the filtered variants: in progress ..")
	filt_plot_vcfstats = " ".join([
	"plot-vcfstats",
	"--prefix", filtered_plots,
	"{0}/bcftools_stats.filtered.txt".format(varcall_dir),
	"2>>", os.path.join(postreports_dir, "plot_vcfstats_filtered_report.txt")])
	subprocess.run(filt_plot_vcfstats, shell=True)

	### ANNOTATE VARIANTS
	## 7. Adding information to the discovered variants
	## At this step we run tools to add information to the discovered variants in our dataset.
	## One of those tools, Funcotator, can be used to add gene-level information to each variant.
	print("9/9 | Funcotator - Adding information to the discovered variants in our dataset: in progress ..")
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
	

	quality_control()
	
	aligning()

	post_processing()

	calling_variants()

if __name__ == "__main__": main()