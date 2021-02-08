from snakemake.utils import R
from itertools import repeat
from os.path import join
import os
import shutil
import pandas
import jsonschema
import math

# report: "report/workflow.rst"

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://512303535801.dkr.ecr.eu-west-1.amazonaws.com/scrnax_conda_env:v1.1.1"

container: "docker://continuumio/miniconda3:4.7.12"

##### load config and sample sheets #####


# validate(config, schema="../schemas/config.schema.yaml")

# samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
# samples.index.names = ["sample_id"]
# validate(samples, schema="../schemas/samples.schema.yaml")

# Globals ---------------------------------------------------------------------

scaleDownThreads = config.get("scaleDownThreads", 1)

# Selftest -------------------------------------------------------------------
if config.get('example', "") == "true":
  configfile: srcdir("../config/conf.selftest.json")
else:
  configfile: srcdir("../config/conf.default.json")
  
# Initial optional parameters -------------------------------------------------
# optionalKeys = {
#   "numHit": 10000,
#   "featureCountArgs": {},
#   "featureCountToDTArgs": {},
#   "STARArgs": {},
#   "DTToH5Args": {}
# }
# for tmpKey, tmpValue in optionalKeys.items():
#   if tmpKey not in config:
#     config[tmpKey] = tmpValue

# Prepare outputs -------------------------------------------------------------
resultFiles = [config["outputDir"] + "ortholog/fixed.gtf"]

if "bamDir" in config:
  if config['bamDir'][len(config['bamDir'])-1] != "/":
    config['bamDir'] = config['bamDir'] + "/"
  IDS, = glob_wildcards(config['bamDir'] + "{sample}.bam")
  resultFiles = [
    resultFiles,
    config["outputDir"] + "denovo/fixed.gtf",
    config["outputDir"] + "combined.fixed.gtf"
  ]
else:
  IDS = ""
  config["bamDir"] = ""

# How the fastqs shall be counted, i.e., by original GTF, or from denovo assembled GTF, or from ortholog GTF, or combined GTF.

if "countBy" in config:
  if config["countBy"] == "combined":
    config["featureCountGTF"] = config["outputDir"] + "combined.fixed.gtf"
  else:
    config["featureCountGTF"] = config["outputDir"] + config["countBy"] + "/fixed.gtf"
elif config["bamDir"] != "":
    config["countBy"] = "combined"
    config["featureCountGTF"] = config["outputDir"] + "combined.fixed.gtf"
else:
  config["countBy"] = "ortholog"
  config["featureCountGTF"] = config["outputDir"] + "ortholog/fixed.gtf"

# define the parameters for STARSolo
# Generate count matrix if fastq files are provided.

if "fastqDir" in config:
  if config['fastqDir'][len(config['fastqDir'])-1] != "/":
    config['fastqDir'] = config['fastqDir'] + "/"
  FASTQS, = glob_wildcards(config['fastqDir'] + "{fastq}_R1.fastq.gz")
  resultFiles = resultFiles + [config["outputDir"] + "featureCount/" + config["countBy"] + "/" + str(i) + ".h5" for i in FASTQS]
  jsS = open(srcdir("schemas/jsonSchemaStarSolo.json")).read()
  # And validate the json config file
else:
  jsS = open(srcdir("schemas/jsonSchemaWithoutStarSolo.json")).read()
  config["fastqDir"] = ""
# validate json config file
jsSchema = json.loads(jsS)
try:
  jsonschema.validate(config, jsSchema)
except jsonschema.exceptions.ValidationError as e:
  print("well-formed but invalid JSON:", e)
  sys.exit()
except json.decoder.JSONDecodeError as e:
  print("poorly-formed text, not JSON:", e)
  sys.exit()

# set default parameters for STARSolo
starsoloParameters = '''  --soloType {solotype} \
        --soloCBwhitelist {whitelist} \
        --soloCBlen {CBLen} \
        --soloUMIstart {UMIStart} \
        --soloUMIlen {UMILen} \
        --soloBarcodeReadLength 0 \
        --soloStrand {soloStrand} \
	      --soloFeatures {soloFeatures} \
    '''.format(solotype = config["STARArgs"].get("soloType","CB_UMI_Simple"), whitelist = config['whitelist'], CBLen = config['CBLen'], UMIStart = config["UMIStart"], UMILen = config["UMILen"], soloStrand = config["STARArgs"].get("soloStrand","Forward"), soloFeatures = config["STARArgs"].get("soloFeatures","Gene GeneFull SJ Velocyto"))

config["STARArgs"]["outSAMattributes"] = "NH HI nM AS CR UR CB UB GX GN sS sQ sM"

# print(resultFiles)
rule all:
	input: 
	  resultFiles
	  
include: "rules/prepareOriginalGTF.py"
include: "rules/assembleGTF.py"
include: "rules/orthologGTF.py"
include: "rules/countFastq.py"


rule fixOriginGtf:
  input:
    config['refGTF']
  output:
    fixedGtf = config["outputDir"] + "origin/" + "fixed.gtf",
  params:
    fixGTFRscript = srcdir("scripts/fixGffcompareGTF.R")
  conda: 
    "envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      Rscript {params.fixGTFRscript} {input} {output.fixedGtf}
    """

rule fixOrthologGtf:
  input:
    annoGtf = rules.annoOrthologGtf.output.annoGtf,
    mtGtf = config["outputDir"] + "mt.{strand}.gtf"
  output:
    fixedGtf = config["outputDir"] + "ortholog/" + "fixed.{strand}.gtf",
    tmpGtf = temp(config["outputDir"] + "ortholog/anno.{strand}.combined.fixed.tmp.gtf")
  params:
    fixGTFRscript = srcdir("scripts/fixGffcompareGTF.R")
  conda: 
    "envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      cat {input.annoGtf} > {output.tmpGtf}
      cat {input.mtGtf} >> {output.tmpGtf}
      Rscript {params.fixGTFRscript} {output.tmpGtf} {output.fixedGtf}
    """

rule mergeStrandedOrthologGtf:
  input:
    config["outputDir"] + "ortholog/" + "fixed.pos.gtf",
    config["outputDir"] + "ortholog/" + "fixed.neg.gtf"
  output:
    config["outputDir"] + "ortholog/" + "fixed.gtf"
  conda: 
    "envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      cat {input} > {output}
    """
    
rule fixDenovoGtf:
  input:
    annoGtf = rules.annoDevoGtf.output,
    mtGtf = config["outputDir"] + "mt.{strand}.gtf"
  output:
    fixedGtf = config["outputDir"] + "denovo/fixed.{strand}.gtf",
    tmpGtf = temp(config["outputDir"] + "denovo/anno.{strand}.combined.fixed.tmp.gtf")
  params:
    fixGTFRscript = srcdir("scripts/fixGffcompareGTF.R")
  conda: 
    "envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      cat {input.annoGtf} > {output.tmpGtf}
      cat {input.mtGtf} >> {output.tmpGtf}
      Rscript {params.fixGTFRscript} {output.tmpGtf} {output.fixedGtf}
    """

rule mergeStrandedDenovoGtf:
  input:
    config["outputDir"] + "denovo/" + "fixed.pos.gtf",
    config["outputDir"] + "denovo/" + "fixed.neg.gtf"
  output:
    config["outputDir"] + "denovo/" + "fixed.gtf"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      cat {input} > {output}
    """

rule combineAndFixGtf:
  input:
    #rules.fixDenovoGtf.output.fixedGtf,
    #rules.fixOrthologGtf.output.fixedGtf
    denoAnnoGtf = rules.annoDevoGtf.output,
    orthAnnoGtf = rules.annoOrthologGtf.output,
    refGtf = config["outputDir"] + "reference.{strand}.gtf",
    mtGtf = config["outputDir"] + "mt.{strand}.gtf"
  output:
    combinedMergedGtf = config["outputDir"] + "merged.{strand}.gtf",
    combinedMergedGtfNR = config["outputDir"] + "merged.{strand}.nr.gtf",
    tmpGtf = config["outputDir"] + "tmp.{strand}.annotated.gtf",
    fixedGtf = config["outputDir"] + "combined.fixed.{strand}.gtf"
    # gtflist = config["outputDir"] + "gtflist.txt"
  params:
    fixGTFRscript = srcdir("scripts/fixGffcompareGTF.R"),
    # annoOrthPath = rules.annoOrthologGtf.output,
    # annoDevoPath = rules.annoDevoGtf.output,
    outprefix = lambda wildcards: config["outputDir"] + "tmp." + wildcards.strand,
    mergeGtfDir = config["outputDir"]
  conda: 
    "envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads:
    8
  shell:
    """
      cat {input.denoAnnoGtf} {input.orthAnnoGtf} > {output.combinedMergedGtf}
      stringtie --merge -l COMB -g -1000000000 -i -F 0 -T 0 -f 0 {output.combinedMergedGtf} > {output.combinedMergedGtfNR}
      gffcompare -o {params.outprefix} -r {input.refGtf} {output.combinedMergedGtfNR} -K
      cat {input.mtGtf} >> {output.tmpGtf}
      Rscript {params.fixGTFRscript} {output.tmpGtf} {output.fixedGtf}
   """
    
rule mergeStrandedCombinedGtf:
  input:
    config["outputDir"] + "combined.fixed.pos.gtf",
    config["outputDir"] + "combined.fixed.neg.gtf"
  output:
    config["outputDir"] + "combined.fixed.gtf"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 1000
  threads: 1
  shell:
    """
      cat {input} > {output}
    """
