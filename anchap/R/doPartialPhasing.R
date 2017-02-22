doPartialPhasing <- function (probandId, firstSAnalysed, lastSAnalysed, selectedSequencesH1, selectedSequencesH2, chromEnds, map, iGene, chromList) {


  
  getChromStart <- function(chromosomeNo) { if (chromosomeNo==1)
    {value<-as.integer(1)} else {value <- as.integer(chromEnds[chromosomeNo-1]+1);}
  value }

  getChromEnd <- function(chromosomeNo)
  {as.integer(chromEnds[chromosomeNo]);}

  getChromNo <- function(firstS, lastS) {

  chromMap <- map$chromosome[firstS:lastS]; chromNo <- max(chromMap)
  if (chromNo!=min(chromMap))
  print("Surrogate seq on different chromosomes!"); chromNo }

  
    chromNo <- getChromNo(firstSAnalysed, lastSAnalysed);

    
    firstSChrom <- as.integer(getChromStart(chromNo))
    lastSChrom <- as.integer(getChromEnd(chromNo));
    
    noSnpsChrom <- as.integer(lastSChrom - firstSChrom +1);
    firstSAnalysed <- as.integer(firstSAnalysed);
  
    probandId <- as.integer(probandId);
    noInds <- as.integer(length(iGene));
    
    genesVec <-chromList[[chromNo]];

    noH1Seqs <<- dim(selectedSequencesH1)[1]; h1Starts <<- as.numeric(selectedSequencesH1$start); h1Ends <<-  as.numeric(selectedSequencesH1$end); h1Ids <<-  (selectedSequencesH1$id);
    noH2Seqs <<- dim(selectedSequencesH2)[1]; h2Starts <<- as.numeric(selectedSequencesH2$start); h2Ends <<-  as.numeric(selectedSequencesH2$end); h2Ids <<-  (selectedSequencesH2$id);
    
    noH1atLocus <- integer(lastSAnalysed-firstSAnalysed)
    noH2atLocus <- integer(lastSAnalysed-firstSAnalysed)
    haplotype1 <- integer(lastSAnalysed-firstSAnalysed)
    haplotype2 <- integer(lastSAnalysed-firstSAnalysed)
    majority1 <- integer(lastSAnalysed-firstSAnalysed)
    majority2 <- integer(lastSAnalysed-firstSAnalysed)
    genotype <- integer(lastSAnalysed-firstSAnalysed)
    incons1 <- integer(lastSAnalysed-firstSAnalysed)
    incons2  <- integer(lastSAnalysed-firstSAnalysed)   

    lastSAnalysed <- as.integer(lastSAnalysed-1);

    
    #dyn.load("fastComputation/doFastPhasing.so")
    r <- .C("doFastPhasing", noSnpsChrom, firstSChrom, lastSChrom, firstSAnalysed, lastSAnalysed, probandId, noInds, # 1-7
            genesVec, # 8
            noH1Seqs, h1Starts, h1Ends, h1Ids, # 9 - 12
            noH2Seqs, h2Starts, h2Ends, h2Ids, # 13-16
            noH1atLocus, noH2atLocus, # 17-18
            haplotype1, haplotype2, majority1, majority2, genotype, incons1, incons2) # 19-25


    result <- list();    
    result[[1]]<-r[[19]]
    result[[2]]<-r[[20]]

    rm(r)


    result

}
