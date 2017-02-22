writeHapsToVcf <- function(matHapMat, patHapMat, map, idNames, fileName) {

if (is.null( dim(matHapMat))) {
  noSnps <- length(matHapMat)
} else {
  noSnps <- dim(matHapMat)[2]
}
  
# recode alleles to vcf
recodeV <- function(v) {
    vRecoded <- vector()
    iMaj <- (v==3); vRecoded[iMaj] <- "0"
    iMin <- (v==1); vRecoded[iMin] <- "1"
    vRecoded
  }

# form one row of the vcf file
parseLocus <- function(vH1, vH2) {

  vR <- vector();
  noInds <- length(vH1);

  phasedLoci <- (vH1>0)
  # assign arbitrary phase to no-phased heterozygous loci
  vH1[!phasedLoci] <- 3
  vH2[!phasedLoci] <- 1
  

  h1Ind <- seq(from=1, by=4, length.out=noInds)
  h2Ind <- seq(from=3, by=4, length.out=noInds)
  symbolInd <- seq(from=2, by=4, length.out=noInds)
  spaceInd <- seq(from=4, by=4, length.out=noInds)

  vR[h1Ind] <- recodeV(vH1)
  vR[h2Ind] <- recodeV(vH2)
  vR[spaceInd] <- "\t"

  # indicate phase
  phaseSymbolVector<-vector()
  phaseSymbolVector[phasedLoci] <- "|"
  phaseSymbolVector[!phasedLoci] <- "/"
  vR[symbolInd] <- phaseSymbolVector
  
  #
  vR
}

# initialize the vcf file
cat("##fileformat=VCFv4.0\n##source=pseq\n##FILTER=<ID=PASS,Description=\"Passed variant FILTERs\">\n", file=fileName, append=FALSE)

headerString1 <- paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', sep='\t')
headerString2 <- paste(idNames,sep='\t', collapse="\t")
headerString <- paste(headerString1, headerString2 ,sep='\t')
cat(headerString, file=fileName, fill=FALSE, append=TRUE)
cat('\n', file=fileName, fill=FALSE, append=TRUE)

  
 # loop over SNPs
for (s in 1:noSnps) {
#for (s in 1:100) {
  print(paste("Outputting snp", s, "out of", noSnps))

  if (!is.null( dim(matHapMat))) {
  h1Vec <- matHapMat[,s]; h2Vec <- patHapMat[,s];
  } else {
  h1Vec <- matHapMat[s]; h2Vec <- patHapMat[s];
  }
  
  
  chromosomeString <- paste('chr', allData$map[s,'chromosome'],sep="")
  descriptionString <- paste(chromosomeString, allData$map[s,'position'],  allData$map[s,'snp.names'], 'A', 'G', '.', '.', '.', 'GT', sep='\t')
  
  snpString <- parseLocus(h1Vec, h2Vec)

  cat(paste(descriptionString,'\t',sep=""), file=fileName, append=TRUE, fill=FALSE)
  cat(snpString, file=fileName, append=TRUE, fill=FALSE, sep="")
    cat('\n', file=fileName, append=TRUE, fill=FALSE)
  # work out the row header

}


}
