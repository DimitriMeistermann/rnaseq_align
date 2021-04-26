import glob, os, sys, json, pysam, shutil


WORKING_DIR = os.path.dirname(workflow.snakefile)
if config=={}:
	print("Default config file loaded, from " + WORKING_DIR + "/config.json")
	configfile: WORKING_DIR+"/config.json"

## creation of the logs subdirectory
if not os.path.exists(WORKING_DIR+"/log"):
	os.mkdir(WORKING_DIR+"/log")

#put all config variable as variable in the snakefile
for configVar in config:
	if isinstance(config[configVar], str): exec(configVar+"= '"+config[configVar]+"'")
	else: exec(configVar+"="+str(config[configVar]))


## test of the path provided in the config.json file
if not os.path.exists(FASTQ_PATH):
	print("The directory " + FASTQ_PATH + " doesn't exist. Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)
else:
	## If the path ends by /, the / is suppressed
	if ( FASTQ_PATH[-1:] == "/" ):
		FASTQ_PATH = FASTQ_PATH[:-1]

INPUT_FASTQS = glob.glob(FASTQ_PATH+'/*.fastq.gz')

SAMPLESsplitted = [os.path.basename(f).split(".") for f in INPUT_FASTQS]
SAMPLES=[]

#remove .fastq.gz to get sample names
for s in SAMPLESsplitted:
	SAMPLES.append(".".join(s[0:-2]))

if(OUTPUT_PATH[-1] == "/") : OUTPUT_PATH = OUTPUT_PATH[:-1]

## suppress the .R1. and .R2. elements for paired-end fastq files for the alignement processus in SAMPLES
if IS_PAIRED_END:
	SAMPLES = {itemR2 for itemR2 in SAMPLES if (PAIR_END_FILE_PATTERN+"2") not in itemR2}
	SAMPLES = {itemR1.replace((PAIR_END_FILE_PATTERN+"1"),'') for itemR1 in SAMPLES}

#OVERHANG calculation
if READ_LENGTH<1:
	with os.popen("gunzip -c " + INPUT_FASTQS[0] + " | sed -n '2p'","r") as inputRead:
		OVERHANG = len(inputRead.readline()) - 2
else:
	OVERHANG = int(READ_LENGTH) - 1

#forward/reverse adaptator accounting
#if no adptator removal ("FORWARD_ADAPTATOR" and "REVERSE_ADAPTATOR" == ""), input of STAR is directly "FASTQ_PATH"
CUTADAPT_PARAMETER="" # -a forward -A reverse for pairedEnd or just -a forward for singleEnd

if FORWARD_ADAPTATOR == "" and  REVERSE_ADAPTATOR =="":
	STAR_FASTQ_FOLDER = FASTQ_PATH
	USE_CUTADAPT = False
else:
	STAR_FASTQ_FOLDER = OUTPUT_PATH + "/FASTQ_CUTADAPT"
	USE_CUTADAPT = True

if IS_PAIRED_END:
	PAIR_SUFFIX = [PAIR_END_FILE_PATTERN+"1",PAIR_END_FILE_PATTERN+"2"]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + FORWARD_ADAPTATOR + " -A " + REVERSE_ADAPTATOR
else:
	PAIR_SUFFIX = [""]
	if USE_CUTADAPT: CUTADAPT_PARAMETER = "-a " + FORWARD_ADAPTATOR

STAR_INDEX_DIR = os.path.dirname(FASTA_REFERENCE)+"/STAR_INDEX"

if IS_PAIRED_END:  print("Workflow set on paired end mode")
else : print("Workflow set on single end mode")
if USE_CUTADAPT: print("cutadapt will be used on fastq. To disable it, set FORWARD_ADAPTATOR and REVERSE_ADAPTATOR as an empty string in config file")

if os.path.isfile(OUTPUT_PATH+"/log/genomeIsLoaded"): os.remove(OUTPUT_PATH+"/log/genomeIsLoaded")

#for multiple files, create batch of files

if STAR_ALIGN_MULTIPLE_FILE:
	batchLim = MAX_SIZE_PER_MULTIPLE_ALIGN * 10**9
	batchNum=1
	currentSize=0
	currentBatch="batch"+str(batchNum)

	ALIGN_BATCHES={currentBatch:[]}

	for sample in SAMPLES:
		for pair in PAIR_SUFFIX:
			currentSize+=os.path.getsize(FASTQ_PATH+"/"+sample+pair+".fastq.gz")
		ALIGN_BATCHES[currentBatch].append(sample)
		#number of opened file in a process can not exess 1024, samtools split open all possible file in one process
		if (currentSize > batchLim or len(ALIGN_BATCHES[currentBatch]) > 1000)  and sample != SAMPLES[len(SAMPLES)-1]:
			batchNum+=1
			currentSize=0
			currentBatch="batch"+str(batchNum)
			ALIGN_BATCHES[currentBatch]=[]

##############
rule all: 
	input: OUTPUT_PATH+"/results/multiqc_report.html"

rule STAR_INDEX:
	input:
		gtf = GTF_REFERENCE,
		fasta = FASTA_REFERENCE
	output: directory(STAR_INDEX_DIR)
	params:
		cpu=STAR_GENOME_THREADS,
		ram=int(GENOME_INDEX_RAM),
	log:
		out=OUTPUT_PATH+"/log/STAR_INDEX.out",
		err=OUTPUT_PATH+"/log/STAR_INDEX.err",
		starLog=OUTPUT_PATH+"/log/STAR_INDEX.Log.out"
	shell: """
	mkdir {output}
	STAR --runThreadN {params.cpu} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fasta} \\
	--outFileNamePrefix {OUTPUT_PATH}/log/STAR_INDEX. \\
	--sjdbGTFfile {input.gtf} --sjdbOverhang {OVERHANG}  --limitGenomeGenerateRAM {params.ram} 1> {log.out} 2> {log.err}
	"""

rule CUTADAPT:
	input: expand(FASTQ_PATH+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	output: expand(OUTPUT_PATH+"/FASTQ_CUTADAPT/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
	params:
		minLen = CUTADAPT_MIN_READ_LEN,
		cpu = THREAD_PER_SAMPLE
	log:
		err=OUTPUT_PATH+"/log/CUTADAPT_{sample}.err"
	shell: """
	cutadapt {CUTADAPT_PARAMETER} -j {params.cpu} --minimum-length {params.minLen}  {input} 2> {log.err} | gzip -c > {output}
	"""	

rule FASTQC:
	input: FASTQ_PATH+"/{sample}{pair}.fastq.gz"
	output: multiext(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc",".zip",".html")
	params:
		outpath=OUTPUT_PATH+"/fastQC",
		cpu = 1
	log:
		out=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.out",
		err=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.err"
	shell: """
	fastqc -o {params.outpath} {input} 1> {log.out} 2> {log.err}
	"""

if PRELOAD_GENOME:
	STAR_LOAD_BEHAVIOR = "--genomeLoad LoadAndKeep"
	genomeIsLoadedFile = temp(OUTPUT_PATH+"/log/genomeIsLoaded")
	rule STAR_LOAD_GENOME:
		input: 
			starIndexDir=STAR_INDEX_DIR
		output:
			OUTPUT_PATH+"/log/genomeIsLoaded"
		log:
			out=OUTPUT_PATH+"/log/STAR_LOAD_GENOME.out",
			err=OUTPUT_PATH+"/log/STAR_LOAD_GENOME.err",
			starLog=OUTPUT_PATH+"/log/STAR_LOAD_GENOME.Log.out",
			starProgressLog=OUTPUT_PATH+"/log/STAR_LOAD_GENOME.progress.out"
		params: cpu=STAR_GENOME_THREADS
		shell: """
		STAR --genomeLoad Remove  --genomeDir {input.starIndexDir} --runThreadN {params.cpu} \\
		--outFileNamePrefix {OUTPUT_PATH}/log/STAR_LOAD_GENOME. 1> {log.out} 2> {log.err}
		
		STAR --genomeLoad LoadAndExit  --runThreadN {params.cpu}  \\
		--genomeDir {input.starIndexDir} \\
		--outFileNamePrefix {OUTPUT_PATH}/log/STAR_LOAD_GENOME.  1> {log.out} 2> {log.err}
		
		touch {output}
		rm {OUTPUT_PATH}/log/STAR_LOAD_GENOME.Aligned.out.sam
		"""
else:
	STAR_LOAD_BEHAVIOR=""
	genomeIsLoadedFile=[]

if STAR_ALIGN_MULTIPLE_FILE:

	def getFASTQSofAlignBatch(wildcards): 
		return expand(STAR_FASTQ_FOLDER+"/{sample}{pair}.fastq.gz",\
		sample=ALIGN_BATCHES[wildcards.alignBatch],pair=PAIR_SUFFIX)

	rule BUILD_MANIFEST:
		input:
			getFASTQSofAlignBatch
		output:
			OUTPUT_PATH+"/log/star_manifest_{alignBatch}.txt"
		run:
			with open(str(output),"w") as writer:
				for sample in ALIGN_BATCHES[wildcards.alignBatch]:
					base_name = STAR_FASTQ_FOLDER+"/"+sample
					if IS_PAIRED_END:
						writer.write(base_name+PAIR_SUFFIX[0]+".fastq.gz\t"+base_name+PAIR_SUFFIX[1]+".fastq.gz\t"+sample+"\n")
					else:
						writer.write(base_name+".fastq.gz\t-\t"+sample+"\n")

	rule STAR_ALIGN:
		input: 
			starIndexDir=STAR_INDEX_DIR,
			manifest=OUTPUT_PATH+"/log/star_manifest_{alignBatch}.txt"
		output:
			temp(OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.Aligned.out.bam")
		log:
			out=OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.out",
			err=OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.err",
			starLog=OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.Log.out",
			starProgressLog=OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.progress.out"
		params: cpu = STAR_GENOME_THREADS
		shell: """
		STAR --runThreadN {params.cpu} --genomeDir {input.starIndexDir} --readFilesCommand zcat \\
		--readFilesManifest {input.manifest} --outSAMattributes NH HI AS nM RG \\
		--outFileNamePrefix {OUTPUT_PATH}/log/STAR_ALIGN_{wildcards.alignBatch}. --outSAMtype BAM Unsorted \\
		{STAR_LOAD_BEHAVIOR} 1> {log.out} 2> {log.err}
		"""

	rule SPLIT_BAM:
		input: OUTPUT_PATH+"/log/STAR_ALIGN_{alignBatch}.Aligned.out.bam"
		output: touch(OUTPUT_PATH+"/log/{alignBatch}_splitIsDone"), 
		params: cpu = THREAD_PER_SAMPLE
		shell: """			
		samtools split -f '{OUTPUT_PATH}/log/%!.bam' -@{params.cpu} {input}
		"""

	rule FAKE_RULE_MULTIPLE_FILE:
		input: expand(OUTPUT_PATH+"/log/{alignBatch}_splitIsDone",alignBatch=ALIGN_BATCHES)
		output: expand(OUTPUT_PATH+"/BAM/{sample}.bam", sample=SAMPLES)
		run:
			missingFile=[]
			for alignBatch in ALIGN_BATCHES:
				missingFileInBAM = False

				for sample in ALIGN_BATCHES[alignBatch]:
					bamFile = OUTPUT_PATH+"/log/"+sample+".bam"
					if not os.path.exists(bamFile):
						missingFileInBAM=True
						missingFile.append(bamFile)
				
				if missingFileInBAM: os.remove(OUTPUT_PATH+"/log/"+alignBatch+"_splitIsDone")
			
			if len(missingFile)>0:
				exit("BAM splitting step is a failure, the following BAM are missing: "+ " ".join(missingFile))
			
			#no error
			for sample in SAMPLES:
				os.rename(OUTPUT_PATH+"/log/"+sample+".bam", OUTPUT_PATH+"/BAM/"+sample+".bam")

else:
	rule STAR_ALIGN:
		input: 
			fastq=expand(STAR_FASTQ_FOLDER+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX),
			starIndexDir=STAR_INDEX_DIR,
			genomeIsLoaded=expand("{file}", file=genomeIsLoadedFile) #optional input
		output:
			OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.Aligned.out.bam"
		log:
			out=OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.out",
			err=OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.err",
			starLog=OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.Log.out",
			starProgressLog=OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.progress.out"
		params: cpu = THREAD_PER_SAMPLE
		shell: """
		
		STAR --runThreadN {params.cpu} --genomeDir {input.starIndexDir} --readFilesCommand zcat --readFilesIn {input.fastq} \\
		--outFileNamePrefix {OUTPUT_PATH}/log/STAR_ALIGN_{wildcards.sample}. --outSAMtype BAM Unsorted \\
		{STAR_LOAD_BEHAVIOR} 1> {log.out} 2> {log.err}

		"""
	
	rule MOVE_BAM:
		input: OUTPUT_PATH+"/log/STAR_ALIGN_{sample}.Aligned.out.bam"
		output: OUTPUT_PATH+"/BAM/{sample}.bam"
		params: cpu = 1
		shell: """
			mv {input} {output}
		"""

if PRELOAD_GENOME:
	rule STAR_UNLOAD_GENOME:
		input:
			bams = expand(OUTPUT_PATH+"/BAM/{sample}.bam", sample=SAMPLES),
			genomeIsLoaded = OUTPUT_PATH+"/log/genomeIsLoaded"
		log:
			out=OUTPUT_PATH+"/log/STAR_UNLOAD_GENOME.out",
			err=OUTPUT_PATH+"/log/STAR_UNLOAD_GENOME.err",
			starLog=OUTPUT_PATH+"/log/STAR_UNLOAD_GENOME.Log.out",
			starProgressLog=OUTPUT_PATH+"/log/STAR_UNLOAD_GENOME.progress.out"
		params: cpu=STAR_GENOME_THREADS
		shell: """
		STAR --genomeLoad Remove  --runThreadN {params.cpu}  --genomeDir {input.starIndexDir} \\
		--outFileNamePrefix {OUTPUT_PATH}/log/STAR_UNLOAD_GENOME.  1> {log.out} 2> {log.err}
		"""

if DEDUP_UMI:
	rule SORT_BAM:
		input: OUTPUT_PATH+"/BAM/{sample}.bam"
		output: temp(OUTPUT_PATH+"/SORTED_BAM/{sample}.bam")
		params: cpu = THREAD_PER_SAMPLE
		shell: """
		samtools sort -@{params.cpu} -o {output} {input} 
		"""

	rule INDEX_BAM:
		input: OUTPUT_PATH+"/SORTED_BAM/{sample}.bam"
		output: temp(OUTPUT_PATH+"/SORTED_BAM/{sample}.bam.bai")
		params: cpu = THREAD_PER_SAMPLE
		shell: """
		samtools index -@{params.cpu} {input}
		"""

	BAM_ALIGN_FOLDER = "DEDUP_BAM"
	rule DEDUP_UMI:
		input:
			bam=OUTPUT_PATH+"/SORTED_BAM/{sample}.bam",
			bai=OUTPUT_PATH+"/SORTED_BAM/{sample}.bam.bai"
		output: OUTPUT_PATH+"/DEDUP_BAM/{sample}.bam"
		log:
			out=OUTPUT_PATH+"/log/DEDUP_UMI_{sample}.out",
			err=OUTPUT_PATH+"/log/DEDUP_UMI_{sample}.err"
		params:
			paired="--paired" if IS_PAIRED_END else "",
			cpu = THREAD_PER_SAMPLE
		shell: """
			umi_tools dedup --no-sort-output --stdin={input.bam} k:{params.paired} --log={log.out}  --output-stats={OUTPUT_PATH}/log/DEDUP_UMI_{wildcards.sample} | samtools sort -n -@{params.cpu} -o {output} 2> {log.err}
		"""
	
else: BAM_ALIGN_FOLDER = "BAM"
	
rule HTSEQ_COUNT:
	input: OUTPUT_PATH+"/"+BAM_ALIGN_FOLDER+"/{sample}.bam"
	output: OUTPUT_PATH+"/counts/{sample}.counts"
	params:
		gtf = GTF_REFERENCE,
		featureID = FEATURE_ID,
		featureType = FEATURE_TYPE,
		cpu = 1
	log: err=OUTPUT_PATH+"/log/HTSEQ_COUNT_{sample}.err"
	shell: """
	htseq-count -f bam -t {params.featureType} -i {params.featureID} -s no {input} {params.gtf} 1> {output} 2> {log.err}
	"""
	
rule COUNTS_TABLE:
	input: expand(OUTPUT_PATH+"/counts/{sample}.counts",sample=SAMPLES)
	output:
		table = OUTPUT_PATH+"/results/rawCountsTable.tsv",
		stat = OUTPUT_PATH+"/results/alignStatTable.tsv"
	params: cpu = 1
	log:
		out=OUTPUT_PATH+"/log/COUNTS_TABLE.out",
		err=OUTPUT_PATH+"/log/COUNTS_TABLE.err"
	shell: "Rscript {WORKING_DIR}/countsTable.R {OUTPUT_PATH}  1> {log.out} 2> {log.err}"

rule MULTIQC:
	input:
		fastqc=expand(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc{ext}", sample=SAMPLES,pair=PAIR_SUFFIX,ext=[".zip",".html"]),
		table = OUTPUT_PATH+"/results/rawCountsTable.tsv",
		stat = OUTPUT_PATH+"/results/alignStatTable.tsv",
	output: OUTPUT_PATH+"/results/multiqc_report.html"
	params:
		outpath = OUTPUT_PATH + "/results",
		cpu = 1
	shell: """
	multiqc -f -e general_stats -e tophat -e bowtie2 {OUTPUT_PATH} -o {params.outpath}
	"""
