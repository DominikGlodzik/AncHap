getGtHaplotypes <- function(allData, pedigree, markerRegion) {

# get the ids of genotyped individuals
sequedIds <- allData$idNames

indsSequed <- (as.character(pedigree$id) %in% sequedIds) 
indsWithSequedFathers <- (as.character(pedigree$fa) %in% sequedIds) 
indsWithSequedMothers <- (as.character(pedigree$mo) %in% sequedIds)


######
# check the pedigree
# loop over lines of pedigree with sequenced parents
for (i in which(indsSequed & indsWithSequedFathers)) {
  
  oi <- which(as.character(pedigree[i,'id'])==sequedIds)  
  # initially the two haplotypes are the same
  probandsGenotype <- allData$iGene[[oi]][markerRegion]
  if (indsWithSequedFathers[i]) {
    fatherId <- which(sequedIds == pedigree$fa[i])
    fathersGenotype <- allData$iGene[[fatherId]][markerRegion]
    # count opposing alleles between proband and parents
    oppProbFa <- sum(((probandsGenotype==1) & (fathersGenotype==3))|((probandsGenotype==3) & (fathersGenotype==1)))
    #if (oppProbFa>100) {
    #  kimmoWithSeqedFathers[i] <- FALSE
    #}
    print(paste('oppProbFa',oppProbFa))
  }
}
######


p <- 1
matHapList <- list()
patHapList <- list()
idArray <- vector()
phasingSuccess <- vector()
datasetId <- vector()
globalInd <- vector()

                                        # for individuals with both parents genotyped
# (by pedigree row)
for (i in which(indsSequed & indsWithSequedFathers & indsWithSequedMothers)) {

  print(as.character(pedigree[i,'id']))
  print(paste(i,"th proband"))

                                        # oi^{th} index among the 749 individuals
  oi <- which(as.character(pedigree[i,'id'])==sequedIds)
  globalInd[i] <- oi
                                        # initially the two haplotypes are the same
  probandsGenotype <- allData$iGene[[oi]][markerRegion]
  probandHomozygous <- (probandsGenotype==1) | (probandsGenotype==3)
  probandHeterozygous <- (probandsGenotype==2)

  # prepare the two haplotypes
  patHap <- integer(length(probandsGenotype))
  matHap <- integer(length(probandsGenotype))
    
  patHap[probandHomozygous] <- probandsGenotype[probandHomozygous]
  matHap[probandHomozygous] <- probandsGenotype[probandHomozygous]
   
    # recovering the parents genotypes
  fatherId <- which(sequedIds == pedigree$fa[i])
  fathersGenotype <- allData$iGene[[fatherId]][markerRegion]
  motherId <- which(sequedIds == pedigree$mo[i])
  mothersGenotype <- allData$iGene[[motherId]][markerRegion]

    # count opposing alleles between proband and parents
  oppProbFa <- sum(((probandsGenotype==1) & (fathersGenotype==3))|((probandsGenotype==3) & (fathersGenotype==1)))
  oppProbMo <- sum(((probandsGenotype==1) & (mothersGenotype==3))|((probandsGenotype==3) & (mothersGenotype==1)))
  print(paste('oppProbFa',oppProbFa))
  print(paste('oppProbMo',oppProbMo))
  
  fatherHomozygous <- (fathersGenotype==1)|(fathersGenotype==3)
  patHap[fatherHomozygous] <- fathersGenotype[fatherHomozygous]
  matHap[probandHeterozygous & fatherHomozygous] <- oppAllele(fathersGenotype[probandHeterozygous & fatherHomozygous])
    
  motherHomozygous <- (mothersGenotype==1)|(mothersGenotype==3)
  matHap[motherHomozygous] <- mothersGenotype[motherHomozygous]
  patHap[probandHeterozygous & motherHomozygous] <- oppAllele(mothersGenotype[probandHeterozygous & motherHomozygous])
    
  # evaluate the phasing wrt to the individual's genotype
  homoLoci <- ((probandsGenotype==1)|(probandsGenotype==3))
    
                                        # check the haplotypes are homozygous where they need to be
  haploHomo <- (matHap[homoLoci]==patHap[homoLoci])
                                        # check that at homozygous loci, the haplotypes are consistent with the genotype
  consistentHaplo <- (probandsGenotype[homoLoci]==patHap[homoLoci])
  consistentHaploMat <- (probandsGenotype[homoLoci]==matHap[homoLoci])
  
    #print(sum(!haploHomo)/length(homoLoci))
    #print(sum(!consistentHaplo))
    #print(sum(!consistentHaploMat))
    

  matHapList[[p]] <- matHap
  patHapList[[p]] <- patHap
  idArray[p] <- as.character(pedigree[i,'id'])
  p <- p + 1 
                                          # discrepancies of haplotypes from original genotypes   
  phasingSuccess[i] <- (1 - sum(matHap==0)/length(markerRegion))


}

result <- list()
result$matHapList <- matHapList;
result$patHapList <- patHapList;
result$idArray <- idArray;
# in 1:749
result$globalInd <- globalInd;
result
} # end the function

# aux function
oppAllele <- function (v) {
  r = v;
  r[v==1]<-3
  r[v==3]<-1
  r
}
