doStitching2 <- function(individualPhased, surrMat, iGene, chromList, surrogateMargin, excludedParents, overlapT, matchT, map, chromEnds) {
cat("Fast stitching method")
  
draftH1 <- integer(length(iGene[[1]]))
draftH2 <- integer(length(iGene[[1]]))
  
getSeqsWithin <- function(seq, selectedSequences) {
  sStart <- seq$start
  sEnd <- seq$end
                                        # start before end
                                        # end after start
  chosenSeqs <- selectedSequences$start<=sEnd & selectedSequences$end>=sStart
  chosenSeqs <- selectedSequences[chosenSeqs, ]
}

countOppHomoR <- function(s1,s2) {
  if (length(s1)!=length(s2)) {
    stop();
  }
  zeros1 <- (s1==0)
  zeros2 <- (s2==0)
  oppHomos <- (abs(s1-s2)==2) & (!zeros1) & (!zeros2)
  sum(oppHomos)
}

getOverlapLength <- function (individualPhased, firstH, secondH) {
  overlapStart <- max(firstH$start, secondH$start)
  overlapEnd <- min(firstH$end, secondH$end)-1

  noHeteroInOverlap <- sum(iGene[[individualPhased]][overlapStart:overlapEnd]==2)  
  noHeteroInOverlap
}

assignToHaplotypes <- function(allSurrogate) {

  
  countEmpty <- 0; countMatchBoth <- 0; countMatchH1 <- 0; countMatchH2 <- 0; countMatchNone <- 0; countNoOverlap <- 0;
  hapMatching <- data.frame()

  allSurrCount <- dim(allSurrogates)[1]

  ## LOOP OVER THE SURROGATES, starting with the longest one
  for (s in 1:allSurrCount) {
    if ((s %% 50)==0) {print(paste("surrogate sequence ", s, " of ",allSurrCount)) }

    # the sequence currently considered 
    cSurr <- allSurrogates[s,]

    sH1 <- getSeqsWithin(cSurr, selectedSequencesH1);
    sH2 <- getSeqsWithin(cSurr, selectedSequencesH2);  

                                        # phase the haplotypes with the sequences classified up till now
    haplotypes <- doPartialPhasing(individualPhased,cSurr$start, cSurr$end, sH1, sH2,  chromEnds, map, iGene, chromList)  
    genotype <- iGene[[individualPhased]][cSurr$start : (cSurr$end-1)]
    # key loci - heterozygous loci that had been phased
    keyLoci <- sum((genotype==2) & (haplotypes[[1]]!=0))
    
    if (keyLoci == 0) {
      selectedSequencesH1<<-rbind(selectedSequencesH1, cSurr)
      countEmpty <-  countEmpty + 1;
    }
    else if (keyLoci>overlapT){

        # try matching
        # matching haplotype 1
        # match the haplotype phased with seqs so far against the genotype of the potential surrogate
        oppHomo1 <- countOppHomoR(iGene[[cSurr$id]][cSurr$start: (cSurr$end-1)], haplotypes[[1]] )
        matchH1 <- ((oppHomo1/keyLoci)<matchT)
        # matching haplotype 2
        # match the haplotype phased with seqs so far against the genotype of the potential surrogate
        oppHomo2 <- countOppHomoR(iGene[[cSurr$id]][cSurr$start: (cSurr$end-1)], haplotypes[[2]])
        matchH2 <- ((oppHomo2/keyLoci)<matchT)

                  # how the new sequence matches the haplotypes
        aMatch <- data.frame(matchH1Rate=(oppHomo1/keyLoci), matchH2Rate=(oppHomo2/keyLoci), matchH1=matchH1, matchH2=matchH2, keyLoci = keyLoci, start=cSurr$start, end=cSurr$end )
        hapMatching <- rbind(hapMatching, aMatch)

          
          if (matchH1 & matchH2) {
            # both matching: error
            countMatchBoth <- countMatchBoth +1
            matchingBothSequences <<- rbind(matchingBothSequences, cSurr)
          }
          else {
            if (matchH1) {
              # matching H1
              # sequence accepted
              selectedSequencesH1<<-rbind(selectedSequencesH1, cSurr)
              countMatchH1 <- countMatchH1 + 1;
            }
            else if (matchH2 & (keyLoci>overlapT)) {
              # matching H2
              # sequence accepted
              selectedSequencesH2<<-rbind(selectedSequencesH2, cSurr)
              countMatchH2 <- countMatchH2 + 1;
            }
            else {
              # not matching any
              rejectedSequences <<- rbind(rejectedSequences, cSurr)
              countMatchNone <- countMatchNone + 1;

              cSeqH1 <- extractIBDnoChrom(iGene[[cSurr$id]][cSurr$start: (cSurr$end-1)], haplotypes[[1]], 100)
              if (dim(cSeqH1)[1]>0) {
                cSeqH1$start <- cSeqH1$start + cSurr$start -1
                cSeqH1$end <- cSeqH1$end + cSurr$start -1
                cSeqH1 <- cbind(cSeqH1, id=cSurr$id)
                partialSequencesH1 <<- rbind(partialSequencesH1, cSeqH1)
              }

              cSeqH2 <- extractIBDnoChrom(iGene[[cSurr$id]][cSurr$start: (cSurr$end-1)], haplotypes[[2]], 100)
              if (dim(cSeqH2)[1]>0) {
                cSeqH2$start <- cSeqH2$start + cSurr$start -1
                cSeqH2$end <- cSeqH2$end + cSurr$start -1
                cSeqH2 <- cbind(cSeqH2, id=cSurr$id)
                partialSequencesH2 <<- rbind(partialSequencesH2, cSeqH2)
              }
            }
          }

        

    }
    else {
        # not enough overlap
        # gray area
        # no matching
        countNoOverlap <- countNoOverlap + 1;
        noOverlapSequences <<- rbind(noOverlapSequences,cSurr)
    }

  }

  result<-list();
  result$matchingSummary <- data.frame ( allSurrCount =  allSurrCount,  countEmpty = (countEmpty/ allSurrCount), countMatchBoth = (countMatchBoth/ allSurrCount), countMatchH1 = (countMatchH1/ allSurrCount), countMatchH2 =  (countMatchH2/ allSurrCount), countMatchNone = (countMatchNone/ allSurrCount),  countNoOverlap = (countNoOverlap/ allSurrCount))
  result$sequenceMatching <- hapMatching

  if ((dim(hapMatching)[1])>0) {
    rejected <- (!hapMatching$matchH1) & (!hapMatching$matchH2)
    result$rejectedSeqs <- hapMatching[rejected,]
  }


  result
}


  
      selectedSequencesH1 <<- data.frame()
      selectedSequencesH2 <<- data.frame()
      rejectedSequences <<- data.frame();
      matchingBothSequences <<- data.frame();
      noOverlapSequences <<- data.frame();

      partialSequencesH1 <<- data.frame();
      partialSequencesH2 <<- data.frame();

      source("utils/getAllSurrogatesOf.R")

      allSurrogates <<- getAllSurrogatesOf(surrMat, individualPhased, surrogateMargin)   

      notExcludedParents <- !(allSurrogates$id %in% excludedParents);
      allSurrogates <<- allSurrogates[notExcludedParents,]

 

      sortedSurrogates <- data.frame()
      if ((dim(allSurrogates)[1])>0) {
        sortedSurrogates <- allSurrogates[order(-allSurrogates$end+allSurrogates$start) , ]
      }


    gc()

      result <- list();
      result[['summary']] <- assignToHaplotypes(sortedSurrogates)         
      result[['selectedSequencesH1']] <- selectedSequencesH1
      result[['selectedSequencesH2']] <- selectedSequencesH2


      result


    }

