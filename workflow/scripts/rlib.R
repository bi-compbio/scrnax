DTToSparseMatrix = function(x){
  x[, c("feature", "cell") := .(factor(feature), factor(cell))]
  tmpM = Matrix::sparseMatrix(i = as.numeric(x$feature), j = as.numeric(x$cell), x = x$value)
  rownames(tmpM) = levels(x$feature)
  colnames(tmpM) = levels(x$cell)
  return(tmpM)
}

toSparseMatrix = function(x, nColSubset = 1000){
  # tmp = list();
  result = NULL;
  i = 1;
  k = 1;
  while(i <= dim(x)[2]){
    startCol = i;
    endCol = ifelse(i + nColSubset > dim(x)[2], dim(x)[2], nColSubset + i);
    print(paste0("Processing column ", startCol, " to ", endCol));
    # tmp[[k]] = Matrix(as.matrix(x[, startCol:endCol]), sparse = T);
    tmp = Matrix(as.matrix(x[, startCol:endCol]), sparse = T);
    k = k + 1;
    i = endCol + 1;
    if(is.null(result)){
      result = tmp
    } else {
      result = cbind(result, tmp);
    }
  }
  # result = do.call("cBind", tmp);
  if(!is.null(rownames(x))){
    rownames(result) = rownames(x);  
  }
  if(!is.null(colnames(x))){
    colnames(result) = colnames(x);  
  }
  return(result);
}

readStarSoloCount = function(x, sampleID = "starsolo", names.delim = "__", whitelist = NULL, nthread = 1, umidist = 1, cbdist = 1, topCBs = 5000){
  if(!require("Rncc",character.only = TRUE)){
    devtools::install_github(repo = c('bowhan/Rncc'), upgrade = F, quiet = T, force = F)
  } 
  cat(paste0("Reading ", x, "\n"))
  countDT = fread(x, sep = "\t")
  # if the CB or UMI has been corrected by STARsolo
  cbdist = ifelse(any(grepl("CB", names(countDT))), 0, cbdist)
  
  umidist = ifelse(any(grepl("UB", names(countDT))), 0, umidist)
  
  names(countDT) = c("cell", "UMI", "feature", "value")
  if(nchar(countDT[, cell][1])<4){
    # for smartseq2 data, the dataID would be assigned as cell barcode
    countDT[, c("cell", "UMI") := .(NULL, NULL)]
    countDT[, cell := sampleID]
    countDT[, UMI := "AAAAAAAAAA"]
    cbdist = 0
    umidist = 0
  }
  cat(paste0("Correcting cell barcodes ...\n"))
  if(!is.null(whitelist) && file.exists(whitelist) && cbdist != 0){
    wl = fread(whitelist, header = F)
    
    tmp = data.table("cell" = unique(countDT[, cell]))
    cbindex = stringdist::amatch(tmp[, cell], wl[, V1], method = "hamming", maxDist = cbdist, nthread = nthread)
    tmp[, corCell := ifelse(is.na(cbindex), "None", wl[, V1][cbindex])]
    tmp = tmp[corCell != "None", ]
    countDT = merge(tmp, countDT, by = "cell")
    countDT[, cell := corCell]
    countDT[, corCell := NULL]
    # countDT[, cell:= ifelse(is.na(cbindex), "None", wl[, V1][cbindex])]
    # countDT = countDT[cell != "None", ]
    if(dim(countDT)[1]==0){
      stop("No correct cell barcode has been identified!")
    }
  }
  cat(paste0("Correcting UMI ...\n"))
  if(umidist == 0){
    countDT[, UMIGroup := 0:(.N-1), by = c("cell", "feature")]
    # countDTUMI = countDT[, .("value" = .N), by = c("cell", "feature")]
  } else {
    countDT[, UMIGroup := Rncc::RcppConnectedComponents(UMI, umidist), by = c("cell", "feature")]
    # countDTUMI = countDT[, .("value" = Rncc::RcppNumberConnectedComponents(UMI, umidist)), by = c("cell", "feature")]  
  }
  countDTUMI = countDT[, .("value" = max(UMIGroup) + 1), by = c("cell", "feature")]
  countMtx = DTToSparseMatrix(countDTUMI)
  
  CBRank = countDTUMI[, .("totUMI" = sum(value)), by = c("cell")]
  CBRank[, UMIRank := frank(totUMI * -1, ties.method = "random")]
  countDTTopCB = countDT[cell %in% CBRank[UMIRank < topCBs, cell], ]
  # countDTTopCB[, UMI := factor(UMI)]
  # countDTTopCB[, c("cell", "feature", "UMI") := lapply(.SD, as.factor), .SDcol = c("cell", "feature", "UMI")]
  
  # assign the representive UMI to each UMI group
  countDTTopCB = countDTTopCB[, .("UMI" = UMI[1], "value" = sum(value)), by = c("cell", "UMIGroup", "feature")]
  
  # if(!is.null(sampleID)){
  #   colnames(countMtx) = paste0(sampleID, names.delim, colnames(countMtx))
  # }
  cat(paste0("Done.\n"))
  result = list(
    "countMtx" = countMtx,
    "UMIInfo" = countDTTopCB
  )
  return(result)
}

#' @export
write10xh5 <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, gene.type="Gene Expression", overwrite=FALSE, genome="unknown", version=c("2", "3"))
  # Writes a count matrix to the specified path in the 10x style.   
  # This allows us to easily create things for testing read10xCounts. 
  # 
  # written by Aaron Lun
  # created 6 January 2018
{
  # Doing all the work on a temporary location next to 'path', as we have permissions there.
  # This avoids problems with 'path' already existing.
  temp.path <- tempfile(tmpdir=dirname(path)) 
  on.exit({ 
    if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
  })
  
  # Checking the values.
  if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
    stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
  }
  if (ncol(x)!=length(barcodes)) { 
    stop("'barcodes' must of of the same length as 'ncol(x)'")
  }
  
  # Determining what format to save in.
  version <- match.arg(version)
  .write_hdf5(temp.path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
  
  # We don't put this at the top as the write functions might fail; 
  # in which case, we would have deleted the existing 'path' for nothing.
  if (overwrite) {
    unlink(path, recursive=TRUE)
  } else if (file.exists(path)) { 
    stop("specified 'path' already exists")
  }
  file.rename(temp.path, path)
  return(invisible(TRUE))
}

#' @importFrom rhdf5 h5createFile h5createGroup h5write
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.write_hdf5 <- function(path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version="3") {
  rhdf5::h5createFile(path)
  
  if (version=="3") {
    group <- "matrix"
  } else {
    group <- genome
  }
  rhdf5::h5createGroup(path, group)
  
  rhdf5::h5write(barcodes, file=path, name=paste0(group, "/barcodes"))
  
  # Saving feature information.
  if (version=="3") {
    rhdf5::h5createGroup(path, file.path(group, "features"))
    rhdf5::h5write(gene.id, file=path, name=paste0(group, "/features/id"))
    rhdf5::h5write(gene.symbol, file=path, name=paste0(group, "/features/name"))
    
    rhdf5::h5write(rep(gene.type, length.out=length(gene.id)),
                   file=path, name=paste0(group, "/features/feature_type"))
    
    rhdf5::h5write(rep(genome, length.out=length(gene.id)),
                   file=path, name=paste0(group, "/features/genome"))
    
  } else {
    rhdf5::h5write(gene.id, file=path, name=paste0(group, "/genes"))
    rhdf5::h5write(gene.symbol, file=path, name=paste0(group, "/gene_names"))
  }
  
  # Saving matrix information.
  x <- as(x, "dgCMatrix")
  rhdf5::h5write(x@x, file=path, name=paste0(group, "/data"))
  rhdf5::h5write(dim(x), file=path, name=paste0(group, "/shape"))
  rhdf5::h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
  rhdf5::h5write(x@p, file=path, name=paste0(group, "/indptr"))
  
  return(NULL)
}

# Writing them to 10X-like molecular info HDF5 files.
writeMolInfoIn10xh5 = function(path, UMIInfo = NULL){
  # Writes a UMI count DT to the specified path in the 10x molecule info style.   
  # written by Aaron Lun
  # created 6 January 2018
  outFile = path
  h5 <- rhdf5::h5createFile(outFile)
  actual.barcodes <- factor(UMIInfo$cell)
  actual.feature = factor(UMIInfo$feature)
  if(class(UMIInfo$UMI) == "factor"){
    actual.umi = as.character(UMIInfo$UMI)
  } else {
    actual.umi = UMIInfo$UMI
  }
  rhdf5::h5write(as.integer(actual.barcodes) - 1L, outFile, "barcode_idx") # technically should be saved as 64-bit, but not possible here.
  rhdf5::h5write(levels(actual.barcodes), outFile, "barcodes")
  
  rhdf5::h5write(as.integer(actual.feature), outFile, "feature_idx")
  rhdf5::h5write(UMIInfo$value, outFile, "count")
  rhdf5::h5createGroup(outFile, "features")
  rhdf5::h5write(levels(actual.feature), outFile, "features/id")
  rhdf5::h5write(actual.umi, outFile, "umi")
  rhdf5::h5write(rep(1L, nrow(UMIInfo)), outFile, "gem_group")
  return(outFile)
}

geneLen = function(gtfDT, geneID = "gene_name"){
  gtfDT = gtfDT[V3 == "exon", ]
  gtfDT[, geneID := gsub(paste0(".*", geneID, " (.*?);.*"),"\\1",V9,perl=T)]
  gtfDT[, geneID := gsub('"', "", geneID, perl = T)]
  gr=GenomicRanges::GRanges(
    seqnames=paste(gtfDT[, V1], gtfDT[, geneID],sep="___"), 
    ranges = IRanges::IRanges(gtfDT[, V4], width = gtfDT[,V5] - gtfDT[, V4] + 1), strand = gtfDT[,V7]
  )
  grSlim=GenomicRanges::reduce(gr);
  grSlimDT=data.table(as.data.frame(grSlim))
  # sum up exons per locus
  result = grSlimDT[, sum(end - start + 1), by = seqnames]
  result[, c("chr", "gene") := tstrsplit(seqnames, split = "___")]
  
  # average length per gene
  result = result[, .("geneLen" = mean(V1), "numLoci" = .N), by = gene]
  
  return(result)
}