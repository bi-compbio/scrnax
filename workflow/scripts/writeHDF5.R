suppressPackageStartupMessages(library("rprojroot"))
suppressPackageStartupMessages(library("data.table"))
scriptPath = dirname(thisfile());
suppressPackageStartupMessages(source(file.path(scriptPath, "rlib.R")))
input = commandArgs(T)
infile = input[1]
whitelist = input[2]
outfile = input[3]
nthread = as.numeric(input[4])
sampleID = input[5]
setDTthreads(threads = nthread)
result = readStarSoloCount(infile, sampleID = sampleID, nthread = nthread, whitelist = whitelist)
write10xh5(path = outfile, x = result$countMtx, version = "3")
molInfoFile = ifelse(grepl("h5$", outfile), gsub("h5$", "molInfo.h5", outfile), paste0(outfile, ".molInfo.h5"))
writeMolInfoIn10xh5(path = molInfoFile, UMIInfo = result$UMIInfo)