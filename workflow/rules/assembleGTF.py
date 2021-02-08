# Denovo assembling and merge assembled GTF with original GTF
# Input: Bam files, fixed original GTF
# Output: merged GTF


rule assemblGtf:
  input:
    config['bamDir'] + "{sample}.bam"
  output:
    fullgtf = config["outputDir"] + "denovo/{sample}.scallop.gtf",
    posGtf = config["outputDir"] + "denovo/{sample}.scallop.pos.gtf",
    negGtf = config["outputDir"] + "denovo/{sample}.scallop.neg.gtf"
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      scallop --min_transcript_coverage 10 --min_num_hits_in_bundle 100 --min_splice_bundary_hits 10 \
        --library_type first  -i {input}  -o {output.fullgtf} --verbose 0  > /dev/null
      mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="-") print; }}' {output.fullgtf} > {output.negGtf}
      mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="+") print; }}' {output.fullgtf} > {output.posGtf}
    """

rule mergeDevoGtf:
  input:
    inGtfs = expand(config["outputDir"] + "denovo/{sample}.scallop.{{strand}}.gtf", sample = IDS),
    refGtf = config["outputDir"] + "reference.{strand}.gtf"
  output:
    mergedDenovoGtf = config["outputDir"] + "denovo/merged.{strand}.gtf"
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      stringtie --merge -g 250 {input.refGtf} {input.inGtfs} -l DENO | \
        grep -v '^MT' > {output}
    """

rule annoDevoGtf:
  input:
    inGtf = rules.mergeDevoGtf.output.mergedDenovoGtf,
    refGtf = config["outputDir"] + "reference.{strand}.gtf"
  output:
    config["outputDir"] + "denovo/anno.{strand}.annotated.gtf"
  params:
    outprefix = lambda wildcards: config["outputDir"] + "denovo/anno." + wildcards.strand
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      gffcompare -o {params.outprefix} -r {input.refGtf} {input.inGtf}
    """
