import glob, os, sys, json, pysam

WORKING_DIR = os.path.dirname(workflow.snakefile)
if config=={}:
	print("Default config file loaded, from " + WORKING_DIR + "/config.json")
	configfile: WORKING_DIR+"/config.json"

## creation of the logs subdirectory
if not os.path.exists(WORKING_DIR+"/log"):
	os.mkdir(WORKING_DIR+"/log")

#is the workflow running in paired mode ?
PAIRED_MODE =  config["IS_PAIRED_END"]
PRELOAD_GENOME = config["PRELOAD_GENOME"]
DEDUP_UMI = config["DEDUP_UMI"]

## test of the path provided in the config.json file
if not os.path.exists( config["FASTQ_PATH"] ):
	print("The directory "+config["FASTQ_PATH"]+" doesn't exist. Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)
else:
	## If the path ends by /, the / is suppressed
	if ( config["FASTQ_PATH"][-1:] == "/" ):
		config["FASTQ_PATH"] = config["FASTQ_PATH"][:-1]

INPUT_FASTQS = glob.glob(config["FASTQ_PATH"]+'/*.fastq.gz')

SAMPLESsplitted = [os.path.basename(f).split(".") for f in INPUT_FASTQS]
SAMPLES=[] 

#remove .fastq.gz to get sample names
for s in SAMPLESsplitted:
	SAMPLES.append(".".join(s[0:-2]))

OUTDIR = config["OUTPUT_PATH"]
if(OUTDIR[-1] == "/") : OUTDIR = OUTDIR[:-1]

## suppress the .R1. and .R2. elements for paired-end fastq files for the alignement processus in SAMPLES
if PAIRED_MODE:
	SAMPLES = {itemR2 for itemR2 in SAMPLES if (config["PAIR_END_FILE_PATTERN"]+"2") not in itemR2}
	SAMPLES = {itemR1.replace((config["PAIR_END_FILE_PATTERN"]+"1"),'') for itemR1 in SAMPLES}

#OVERHANG calculation
if config["READ_LENGTH"]<1:
	with os.popen("gunzip -c " + INPUT_FASTQS[0] + " | sed -n '2p'","r") as inputRead:
		OVERHANG = len(inputRead.readline()) - 2
else:
	OVERHANG = int(config["READ_LENGTH"]) - 1

#forward/reverse adaptator accounting
#if no adptator removal ("FORWARD_ADAPTATOR" and "REVERSE_ADAPTATOR" == ""), input of STAR is directly "FASTQ_PATH"
CUTADAPT_PARAMETER="" # -a forward -A reverse for pairedEnd or just -a forward for singleEnd

if config["FORWARD_ADAPTATOR"] == "" and  config["REVERSE_ADAPTATOR"] =="":
	STAR_FASTQ_FOLDER = config["FASTQ_PATH"]
	USE_CUTADAPT = False
else:
	STAR_FASTQ_FOLDER = OUTDIR + "/FASTQ_CUTADAPT"
	USE_CUTADAPT = True

if PAIRED_MODE:
	PAIR_SUFFIX = [config["PAIR_END_FILE_PATTERN"]+"1",config["PAIR_END_FILE_PATTERN"]+"2"]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + config["FORWARD_ADAPTATOR"] + " -A " + config["REVERSE_ADAPTATOR"]
else:
	PAIR_SUFFIX = [""]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + config["FORWARD_ADAPTATOR"]

STAR_INDEX_DIR = os.path.dirname(config["FASTA_REFERENCE"])+"/STAR_INDEX"

if PAIRED_MODE:  print("Workflow set on paired end mode")
else : print("Workflow set on single end mode")
if USE_CUTADAPT: print("cutadapt will be used on fastq. To disable it, set FORWARD_ADAPTATOR and REVERSE_ADAPTATOR as an empty string in config file")

if os.path.isfile(OUTDIR+"/log/genomeIsLoaded"): os.remove(OUTDIR+"/log/genomeIsLoaded")



##############
rule all: 
	input: OUTDIR+"/results/multiqc_report.html"

rule STAR_INDEX:
	input:
		gtf = config["GTF_REFERENCE"],
		fasta = config["FASTA_REFERENCE"]
	output: directory(STAR_INDEX_DIR)
	params:
		cpu=config["STAR_GENOME_THREADS"],
		ram=int(config["GENOME_INDEX_RAM"]),
	log:
		out=OUTDIR+"/log/STAR_INDEX.out",
		err=OUTDIR+"/log/STAR_INDEX.err",
		starLog=OUTDIR+"/log/STAR_INDEX.Log.out"
	shell: """
	mkdir {output}
	STAR --runThreadN {params.cpu} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fasta} \\
	--outFileNamePrefix {OUTDIR}/log/STAR_INDEX. \\
	--sjdbGTFfile {input.gtf} --sjdbOverhang {OVERHANG}  --limitGenomeGenerateRAM {params.ram} 1> {log.out} 2> {log.err}
	"""

rule CUTADAPT:
	input: expand(config["FASTQ_PATH"]+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	output: expand(OUTDIR+"/FASTQ_CUTADAPT/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	params:
		minLen = config["CUTADAPT_MIN_READ_LEN"],
		cpu = config["THREAD_PER_SAMPLE"]
	log:
		err=OUTDIR+"/log/CUTADAPT_{sample}.err"
	shell: """
	cutadapt {CUTADAPT_PARAMETER} -j {params.cpu} --minimum-length {params.minLen}  {input} 2> {log.err} | gzip -c > {output}
	"""	

rule FASTQC:
	input: config["FASTQ_PATH"]+"/{sample}{pair}.fastq.gz"
	output: multiext(OUTDIR+"/fastQC/{sample}{pair}_fastqc",".zip",".html")
	params:
		outpath=OUTDIR+"/fastQC",
		cpu = 1
	log:
		out=OUTDIR+"/log/FASTQC_{sample}{pair}.out",
		err=OUTDIR+"/log/FASTQC_{sample}{pair}.err"
	shell: """
	fastqc -o {params.outpath} {input} 1> {log.out} 2> {log.err}
	"""

if PRELOAD_GENOME:
	STAR_LOAD_BEHAVIOR = "--genomeLoad LoadAndKeep"
	genomeIsLoadedFile = temp(OUTDIR+"/log/genomeIsLoaded")
	rule STAR_LOAD_GENOME:
		input: 
			starIndexDir=STAR_INDEX_DIR
		output:
			OUTDIR+"/log/genomeIsLoaded"
		log:
			out=OUTDIR+"/log/STAR_LOAD_GENOME.out",
			err=OUTDIR+"/log/STAR_LOAD_GENOME.err",
			starLog=OUTDIR+"/log/STAR_LOAD_GENOME.Log.out",
			starProgressLog=OUTDIR+"/log/STAR_LOAD_GENOME.progress.out"
		params: cpu=config["STAR_GENOME_THREADS"]
		shell: """
		STAR --genomeLoad Remove  --genomeDir {input.starIndexDir} --runThreadN {params.cpu} \\
		--outFileNamePrefix {OUTDIR}/log/STAR_LOAD_GENOME. 1> {log.out} 2> {log.err}
		
		STAR --genomeLoad LoadAndExit  --runThreadN {params.cpu}  \\
		--genomeDir {input.starIndexDir} \\
		--outFileNamePrefix {OUTDIR}/log/STAR_LOAD_GENOME.  1> {log.out} 2> {log.err}
		
		touch {output}
		rm {OUTDIR}/log/STAR_LOAD_GENOME.Aligned.out.sam
		"""
else:
	STAR_LOAD_BEHAVIOR=""
	genomeIsLoadedFile=[]

if config["STAR_ALIGN_MULTIPLE_FILE"]:
	
	rule BUILD_MANIFEST:
		input:
			fastq=expand(STAR_FASTQ_FOLDER+"/{sample}{pair}.fastq.gz",sample=SAMPLES,pair=PAIR_SUFFIX)
		output:
			OUTDIR+"/log/star_manifest.txt"
		run:
			with open(OUTDIR+"/log/star_manifest.txt","w") as writer:
				for sample in SAMPLES:
					base_name = STAR_FASTQ_FOLDER+"/"+sample
					if PAIRED_MODE:
						writer.write(base_name+PAIR[0]+".fastq.gz\t"+base_name+PAIR[1]+".fastq.gz\t"+sample+"\n")
					else:
						writer.write(base_name+".fastq.gz\t-\t"+sample+"\n")

	rule STAR_ALIGN:
		input: 
			fastq=expand(STAR_FASTQ_FOLDER+"/{sample}{pair}.fastq.gz",sample=SAMPLES,pair=PAIR_SUFFIX),
			starIndexDir=STAR_INDEX_DIR,
			manifest=OUTDIR+"/log/star_manifest.txt"
		output:
			temp(OUTDIR+"/log/STAR_ALIGN.Aligned.out.bam")
		log:
			out=OUTDIR+"/log/STAR_ALIGN.out",
			err=OUTDIR+"/log/STAR_ALIGN.err",
			starLog=OUTDIR+"/log/STAR_ALIGN.Log.out",
			starProgressLog=OUTDIR+"/log/STAR_ALIGN.progress.out"
		params: cpu = config["STAR_GENOME_THREADS"]
		shell: """
		STAR --runThreadN {params.cpu} --genomeDir {input.starIndexDir} --readFilesCommand zcat \\
		--readFilesManifest {input.manifest} --outSAMattributes NH HI AS nM RG \\
		--outFileNamePrefix {OUTDIR}/log/STAR_ALIGN. --outSAMtype BAM Unsorted \\
		{STAR_LOAD_BEHAVIOR} 1> {log.out} 2> {log.err}
		"""

	rule SPLIT_BAM:
		input: OUTDIR+"/log/STAR_ALIGN.Aligned.out.bam"
		output: expand(OUTDIR+"/BAM/{sample}.bam",sample=SAMPLES)
		params: cpu = 1
		run:
			reader = pysam.AlignmentFile(str(input), "rb")
			for outbam in output:
				writer = pysam.AlignmentFile(outbam, "wb", template=reader)
				writer.close()
			
			writer = pysam.AlignmentFile("/dev/null", "wb", template=reader)
			outFile=""
			for read in reader:
				rgtag = read.get_tag("RG:Z")
				if rgtag != outFile  :
					writer.close()
					outFile = rgtag
					writer = pysam.AlignmentFile(OUTDIR+"/BAM/"+outFile+".bam", "wb", template=reader)
				read.set_tag("RG", None)
				f=writer.write(read)

			writer.close()
			reader.close()


else:
	rule STAR_ALIGN:
		input: 
			fastq=expand(STAR_FASTQ_FOLDER+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX),
			starIndexDir=STAR_INDEX_DIR,
			genomeIsLoaded=expand("{file}", file=genomeIsLoadedFile) #optional input
		output:
			OUTDIR+"/log/STAR_ALIGN_{sample}.Aligned.out.bam"
		log:
			out=OUTDIR+"/log/STAR_ALIGN_{sample}.out",
			err=OUTDIR+"/log/STAR_ALIGN_{sample}.err",
			starLog=OUTDIR+"/log/STAR_ALIGN_{sample}.Log.out",
			starProgressLog=OUTDIR+"/log/STAR_ALIGN_{sample}.progress.out"
		params: cpu = config["THREAD_PER_SAMPLE"]
		shell: """
		
		STAR --runThreadN {params.cpu} --genomeDir {input.starIndexDir} --readFilesCommand zcat --readFilesIn {input.fastq} \\
		--outFileNamePrefix {OUTDIR}/log/STAR_ALIGN_{wildcards.sample}. --outSAMtype BAM Unsorted \\
		{STAR_LOAD_BEHAVIOR} 1> {log.out} 2> {log.err}

		"""
	
	rule MOVE_BAM:
		input: OUTDIR+"/log/STAR_ALIGN_{sample}.Aligned.out.bam"
		output: OUTDIR+"/BAM/{sample}.bam"
		params: cpu = 1
		shell: """
			mv {input} {output}
		"""

if PRELOAD_GENOME:
	rule STAR_UNLOAD_GENOME:
		input:
			bams = expand(OUTDIR+"/BAM/{sample}.bam", sample=SAMPLES),
			genomeIsLoaded = OUTDIR+"/log/genomeIsLoaded"
		log:
			out=OUTDIR+"/log/STAR_UNLOAD_GENOME.out",
			err=OUTDIR+"/log/STAR_UNLOAD_GENOME.err",
			starLog=OUTDIR+"/log/STAR_UNLOAD_GENOME.Log.out",
			starProgressLog=OUTDIR+"/log/STAR_UNLOAD_GENOME.progress.out"
		params: cpu=config["STAR_GENOME_THREADS"]
		shell: """
		STAR --genomeLoad Remove  --runThreadN {params.cpu}  --genomeDir {input.starIndexDir} \\
		--outFileNamePrefix {OUTDIR}/log/STAR_UNLOAD_GENOME.  1> {log.out} 2> {log.err}
		"""


if DEDUP_UMI:
	rule INDEX_BAM:
		input: OUTDIR+"/BAM/{sample}.bam"
		output: OUTDIR+"/BAM/{sample}.bam.bai"
		params: cpu = config["THREAD_PER_SAMPLE"]
		shell: """
		samtools index -@{params.cpu} {input}
		"""

	BAM_ALIGN_FOLDER = "DEDUP_BAM"
	rule DEDUP_UMI:
		input:
			bam=OUTDIR+"/BAM/{sample}.bam",
			bai=OUTDIR+"/BAM/{sample}.bam.bai"
		output: OUTDIR+"/DEDUP_BAM/{sample}.bam"
		log:
			out=OUTDIR+"/log/DEDUP_UMI_{sample}.out",
			err=OUTDIR+"/log/DEDUP_UMI_{sample}.err"
		params:
			paired="--paired" if PAIRED_MODE else ""
		shell: """
			umi_tools dedup --no-sort-output --stdin={input.bam} k:{params.paired} --log={log.out} > {output} 2> {log.err}
		"""
	
else: BAM_ALIGN_FOLDER = "BAM"
	
#rule SORT_BAM_READ_NAME:
#	input: OUTDIR+"/"+BAM_ALIGN_FOLDER+"/{sample}.bam"
#	output:  OUTDIR+"/SORTED_BAM/{sample}.bam"
#	params: cpu = config["THREAD_PER_SAMPLE"]
#	shell:"""
#	samtools sort -@{params.cpu} -n {input} -o {output}
#	"""

rule HTSEQ_COUNT:
	input: OUTDIR+"/"+BAM_ALIGN_FOLDER+"/{sample}.bam"
	output: OUTDIR+"/counts/{sample}.counts"
	params:
		gtf = config["GTF_REFERENCE"],
		featureID = config["FEATURE_ID"],
		featureType = config["FEATURE_TYPE"],
		cpu = 1
	log: err=OUTDIR+"/log/HTSEQ_COUNT_{sample}.err"
	shell: """
	htseq-count -f bam -t {params.featureType} -i {params.featureID} -s no {input} {params.gtf} 1> {output} 2> {log.err}
	"""
	
rule COUNTS_TABLE:
	input: expand(OUTDIR+"/counts/{sample}.counts",sample=SAMPLES)
	output:
		table = OUTDIR+"/results/rawCountsTable.tsv",
		stat = OUTDIR+"/results/alignStatTable.tsv"
	params: cpu = 1
	log:
		out=OUTDIR+"/log/COUNTS_TABLE.out",
		err=OUTDIR+"/log/COUNTS_TABLE.err"
	shell: "Rscript {WORKING_DIR}/countsTable.R {OUTDIR}  1> {log.out} 2> {log.err}"

rule MULTIQC:
	input:
		fastqc=expand(OUTDIR+"/fastQC/{sample}{pair}_fastqc{ext}", sample=SAMPLES,pair=PAIR_SUFFIX,ext=[".zip",".html"]),
		table = OUTDIR+"/results/rawCountsTable.tsv",
		stat = OUTDIR+"/results/alignStatTable.tsv",
	output: OUTDIR+"/results/multiqc_report.html"
	params:
		outpath = OUTDIR + "/results",
		cpu = 1
	shell: """
	multiqc -f -e general_stats -e tophat -e bowtie2 {OUTDIR} -o {params.outpath}
	"""
