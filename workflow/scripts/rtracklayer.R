library(rtracklayer)
library(data.table)
library(magrittr)


# Inputs ------------------------------------------------------------------

input = commandArgs(T)
inputUTR = input[1]
chainfile = input[2]
outUTR = input[3]

# Prepare original 3UTR GTF -----------------------------------------------

utr = fread(inputUTR, header = F)
utr[, locstr := paste0(V1, ":", V4, "-", V5, ":", V7)]
utr = utr[V5 - V4 >= 200, ]


# Convert into GRange object ----------------------------------------------

grs = utr %>%
  .[, .("Chrom" = V1, "start" = V4, "end" = V5, "strand" = V7, "name" = locstr)] %>%
  unique %>%
  as(., "GRanges")
seqlevelsStyle(grs) = "UCSC"
names(grs) = mcols(grs)$name

# Map GRange --------------------------------------------------------------

ch = import.chain(chainfile)

newgrs = liftOver(grs, ch)
result = newgrs %>% # list of GRange objects
  sort %>% # sort elements in each GRange object
  reduce(., min.gapwidth = 500L, ignore.strand = F) %>% # collapse neibouring objects < 500bp
  as.data.frame() %>% # convert it to data.frame
  data.table() %>% # convert it to data.table
  # .[, start := start - 1] %>% 
  .[, numLoci := .N, by = group_name] %>% 
  .[numLoci == 1, ]

utrInfo = utr %>% 
  .[, oriWidth := V5 - V4 + 1] %>% 
  .[, c("locstr", "V3", "V9", "oriWidth"), with = F]

options(scipen = 999)
refinedUTR = merge(result, utr, by.x = "group_name", by.y = "locstr") %>% 
  .[abs(log2(oriWidth/width)) < log2(1.2), ] %>% 
  .[, c("seqnames", "V2", "V3", "start", "end", "V6", "strand", "V8", "V9"), with = F] %>% # replace the original coordinates with mapped ones
  fwrite(file = outUTR, col.names = F, sep = "\t", quote = F)
