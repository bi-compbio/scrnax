# Prepare original GTF annotation
# Input: original GTF
# Output: GTFs splitted based on MT genes and strandness

rule prepareNoneMTRef:
  input:
    config["refGTF"]
  output:
    posGtf = config["outputDir"] + "reference.pos.gtf",
    negGtf = config["outputDir"] + "reference.neg.gtf"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      grep -v '^MT' {input} | mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="+") print; }}' > {output.posGtf}
      grep -v '^MT' {input} | mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="-") print; }}' > {output.negGtf}
    """

rule prepareMTRef:
  input:
    config["refGTF"]
  output:
    posGtf = config["outputDir"] + "mt.pos.gtf",
    negGtf = config["outputDir"] + "mt.neg.gtf"
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 4000
  threads: 1
  shell:
    """
      grep '^MT' {input} | mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="+") print; }}' > {output.posGtf} || true
      grep '^MT' {input} | mawk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if($7=="-") print; }}' > {output.negGtf} || true
    """
