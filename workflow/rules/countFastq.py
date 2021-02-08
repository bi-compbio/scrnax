# Generate count matrix based on new GTF
# Input: Fastq files
# Output: count matrix in H5 format

rule mapFastq:
  input:
    R1 = config["fastqDir"] + "{fastq}_R1.fastq.gz",
    R2 = config["fastqDir"] + "{fastq}_R2.fastq.gz"
  output:
    countMtx = config["outputDir"] + "starSolo/{fastq}/{fastq}.Solo.out/Gene/raw/matrix.mtx",
    sortBam = config["outputDir"] + "starSolo/{fastq}/{fastq}.Aligned.sortedByCoord.out.bam"
  params:
    outputPrefix = lambda wildcards: config["outputDir"] + "starSolo/" + wildcards.fastq + "/" + wildcards.fastq + ".",
    whitelist = config["whitelist"],soloType = config["STARArgs"].get("soloType","Droplet"),
    barcodeLen = 0, CBLen = config["CBLen"], UMIStart = config["UMIStart"], UMILen = config["UMILen"],
    genome = config["StarsoloGenome"],
    readFilesCommand = config["STARArgs"].get("readFilesCommand","zcat"),
    outSAMtype = config["STARArgs"].get("outSAMtype","BAM SortedByCoordinate"),
    soloStrand = config["STARArgs"].get("soloStrand","Forward"),
    winAnchorMultimapNmax = config["STARArgs"].get("winAnchorMultimapNmax",2000),
    outFilterMultimapNmax = config["STARArgs"].get("outFilterMultimapNmax",2000),
    outSAMprimaryFlag = config["STARArgs"].get("outSAMprimaryFlag","AllBestScore"),
    outSAMmultNmax = config["STARArgs"].get("outSAMmultNmax",1),
    limitBAMsortRAM = config["STARArgs"].get("limitBAMsortRAM",40000000000),
    limitOutSJoneRead = config["STARArgs"].get("limitOutSJoneRead",10000),
    limitOutSJcollapsed = config["STARArgs"].get("limitOutSJcollapsed",3000000),
    limitIObufferSize = config["STARArgs"].get("limitIObufferSize",300000000),
    outSAMattributes = config["STARArgs"].get("outSAMattributes"),
    additionalArguments = config["STARArgs"].get("additionalArguments", ""),
    starOptions = starsoloParameters
  threads: math.ceil(config["STARArgs"].get("threads",16) * scaleDownThreads)
  conda:
    "../envs/generateH5.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * config["STARArgs"].get("memory",40000)
  shell:
    """
      STAR \
        --genomeDir {params.genome} \
        --runThreadN {threads} \
        --readFilesIn {input.R2} {input.R1}  \
        --readFilesCommand {params.readFilesCommand} \
        --outSAMtype {params.outSAMtype} \
	      --winAnchorMultimapNmax {params.winAnchorMultimapNmax} \
	      --outFilterMultimapNmax {params.outFilterMultimapNmax} \
	      --outSAMprimaryFlag {params.outSAMprimaryFlag} \
	      --outSAMmultNmax {params.outSAMmultNmax} \
	      --limitBAMsortRAM {params.limitBAMsortRAM} \
	      --limitOutSJoneRead {params.limitOutSJoneRead} \
        --limitOutSJcollapsed {params.limitOutSJcollapsed} \
        --limitIObufferSize {params.limitIObufferSize} \
        --outFileNamePrefix {params.outputPrefix} \
        --outSAMattributes {params.outSAMattributes} \
	      {starsoloParameters} \
        {params.additionalArguments}
    """
    
rule featureCount:
  input:
    bamfile = rules.mapFastq.output.sortBam,
    gtf = config["featureCountGTF"]
  output:
    temp(config["outputDir"] + "featureCount/" + config["countBy"] + "/{fastq}.Aligned.sortedByCoord.out.bam.featureCounts.bam")
  params:
    outDir = config["outputDir"] + "featureCount/" + config["countBy"] + "/", strand = config["featureCountArgs"].get("strand",1)
  threads: math.ceil(config["featureCountArgs"].get("threads",16) * scaleDownThreads)
  conda:
    "../envs/generateH5.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * config["featureCountArgs"].get("memory",10000)
  shell:
    """
      featureCounts --fraction -M -s {params.strand} -T {threads} -t exon -g gene_name -a {input.gtf} -o {params.outDir}{wildcards.fastq}.counts.txt {input.bamfile} -R BAM --Rpath {params.outDir}
    """


rule featureCountToDT:
  input:
    rules.featureCount.output
  output:
    config["outputDir"] + "featureCount/" + config["countBy"] + "/" + "{fastq}.tsv.gz"
  params:
    countScript = srcdir("../scripts/countBamFromSTARSolo.sh"),
    outDir = config["outputDir"] + "featureCount/",
    numHit = config['numHit']
  threads: math.ceil(config["featureCountArgs"].get("threads",8) * scaleDownThreads)
  conda:
    "../envs/generateH5.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * config["featureCountArgs"].get("memory",4000)
  shell:
    """
      samtools view {input} | NUMHIT={params.numHit} bash {params.countScript} | pigz -p {threads} -c > {output}
    """

rule DTToH5:
  input:
    rules.featureCountToDT.output
  output:
    countMtx = config["outputDir"] + "featureCount/" +  config["countBy"] + "/{fastq}.h5",
    molInfo = config["outputDir"] + "featureCount/" +  config["countBy"] + "/{fastq}.molInfo.h5"
  params:
    writeToh5 = srcdir("../scripts/writeHDF5.R"), whitelist = config["whitelist"]
  threads: math.ceil(config["DTToH5Args"].get("threads",16) * scaleDownThreads)
  conda:
    "../envs/generateH5.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * config["DTToH5Args"].get("memory",4000)
  shell:
    """
      Rscript {params.writeToh5} {input} {params.whitelist} {output.countMtx} {threads} {wildcards.fastq}
    """
