# Prepare original GTF annotation
# Input: original GTF
# Output: GTFs splitted based on MT genes and strandness

rule extractUTR:
  input: 
    config["goodGTF"]
  output: 
    temp(config["outputDir"] + "ortholog/" + "threeUTR.gtf")
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($3=="three_prime_utr") {{$3 = "exon"; print "chr"$0;}} }}' {input} > {output}
    """
        
rule mapUTR:
  input:
    utrgtf = rules.extractUTR.output,
    liftChain = config["liftChain"]
  output:
    tmpMapGtf = temp(config["outputDir"] + "ortholog/" + "mapUTR.tmp.gtf"),
    mapGtfPos = config["outputDir"] + "ortholog/" + "mapUTR.pos.gtf",
    mapGtfNeg = config["outputDir"] + "ortholog/" + "mapUTR.neg.gtf"
  params:
    minMatch = config['liftMinMatch'],
    rtracklayer = srcdir("../scripts/rtracklayer.R")
  threads: 16
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 30000
  shell:
    """
      # liftOver -gff -minMatch={params.minMatch} {input.utrgtf} {input.liftChain} {output.tmpMapGtf} unmap.bed
      Rscript {params.rtracklayer} {input.utrgtf} {input.liftChain} {output.tmpMapGtf}
      # fix short exons
      mawk 'BEGIN{{FS="\\t";OFS="\\t"}}{{if($5-$4 > 50 && $7 == "+") {{if($1 ~ /^chr/) $1 = substr($1, 4, length($1)); print $0; $3 = "transcript"; print $0;}}}}' {output.tmpMapGtf} >  {output.mapGtfPos}
      mawk 'BEGIN{{FS="\\t";OFS="\\t"}}{{if($5-$4 > 50 && $7 == "-") {{if($1 ~ /^chr/) $1 = substr($1, 4, length($1)); print $0; $3 = "transcript"; print $0;}}}}' {output.tmpMapGtf} >  {output.mapGtfNeg}
    """

rule mergeOrthologGtf:
  input:
    config["outputDir"] + "ortholog/" + "mapUTR.{strand}.gtf",
    config["outputDir"] + "reference.{strand}.gtf"
  output:
    mergedGtf = config["outputDir"] + "ortholog/" + "merged.{strand}.gtf"
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      stringtie -f 0 -T 0 -F 0 -l ORTH --merge -g 250 -o {output.mergedGtf} {input}
    """
    
rule annoOrthologGtf:
  input:
    mergedGtf = rules.mergeOrthologGtf.output.mergedGtf,
    refGtf = config["outputDir"] + "reference.{strand}.gtf"
  output:
    annoGtf = config["outputDir"] + "ortholog/" + "anno.{strand}.annotated.gtf"
  params:
    outprefix = lambda wildcards: config["outputDir"] + "ortholog/" + "anno." + wildcards.strand
  conda: 
    "../envs/refineUTR.yaml"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      gffcompare -o {params.outprefix} -r {input.refGtf} {input.mergedGtf} -K
    """
