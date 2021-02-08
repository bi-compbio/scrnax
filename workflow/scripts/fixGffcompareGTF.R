library(data.table)

args = commandArgs(T)
ingtf = args[1]
outgtf = args[2]

gtf = fread(ingtf, header = F, sep = "\t", quote = '"')


# filter out transcripts that 1) opposite strand or 2) not annotated --------

#gtf = gtf[!grepl('class_code "s"', V9) & !grepl('class_code "x"', V9) & !grepl('class_code "u"', V9) & !grepl('class_code "k"', V9), ]
gtf = gtf[!grepl('class_code "s"', V9) & !grepl('class_code "x"', V9) & !grepl('class_code "u"', V9), ]
gtf[, geneID := gsub(".*gene_id (.*?);.*", "\\1", V9, perl = T)]
idToName = unique(gtf[V3 == "transcript", .(geneName = ifelse( grepl("gene_name", V9), gsub(".*gene_name (.*?);.*", "\\1", V9, perl = T), "NULL"), geneID) ])
idToName[, geneName := ifelse(geneName == "NULL", geneID, geneName)]
idToName = unique(idToName)
gtf = merge(gtf, idToName, by = "geneID", allow.cartesian=TRUE)

gtf[, V9 := ifelse(grepl("gene_id", V9) & !grepl("gene_name", V9), paste0(V9, "gene_name ", geneName), V9)]

fwrite(gtf[, !c("geneName", "geneID"), with = F], file = outgtf, col.names = F, row.names = F, sep = "\t", quote = F)

