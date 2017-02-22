doStitching3 <- function(individualPhased, surrMat, iGene, chromList, excludedParents, overlapT, matchT, map, chromEnds) {

# finds sequences with a chosen sequence
# Inputs
# seq - the reference sequence
# selectedSequences - available sequences
# Outputs
# sequences within  seq
getSeqsWithin <- function(seq, selectedSequences) {
  sStart <- seq$start
  sEnd <- seq$end
                                        # start before end
                                        # end after start
  chosenSeqs <- selectedSequences$start<=sEnd & selectedSequences$end>=sStart
  chosenSeqs <- selectedSequences[chosenSeqs, ]
}

# counts opposing homozygotes between sequences of same length
# Inputs
# s1, s2 - two sequences to compare
# Outputs
# count of of opposing homozygotes between the sequences
countOppHomoR <- function(s1,s2) {
  if (length(s1)!=length(s2)) {
    stop();
  }
  zeros1 <- (s1==0)
  zeros2 <- (s2==0)
  oppHomos <- (abs(s1-s2)==2) & (!zeros1) & (!zeros2)
  sum(oppHomos)
}

# finds the overlap between two surrogate sequencies,
# - count of heterozygous loci of the surrogate parents
# Inputs
# individualPhased - id of the individual being phased
# firstH - first surrogate sequence
# secondH - second surrogate sequence
getOverlapLength <- function (individualPhased, firstH, secondH) {
  overlapStart <- max(firstH$start, secondH$start)
  overlapEnd <- min(firstH$end, secondH$end)-1

  noHeteroInOverlap <- sum(iGene[[individualPhased]][overlapStart:overlapEnd]==2)  
  noHeteroInOverlap
}



assignToHaplotypes <- function(allSurrogate) {

  # alles from the surrogates, voted on haplotype 1
  h1Voting <<- matrix(0,length(iGene[[1]]), 2)

  # add the seq to the hap 1 seqs
  # update the voting based on the new haplotype
  #  haploID = which haplotype the seq had been added to 
  addSeqVoting <- function(seq, haploID) {

    genotype <- iGene[[seq$id]][cSurr$start : (cSurr$end-1)]
    homoLoci1 <- which(genotype==1);
    homoLoci3 <- which(genotype==3);

    if (haploID==1) {
      h1Voting[seq$start-1+homoLoci1, 1] <<- h1Voting[seq$start-1+homoLoci1, 1] + 1;
      h1Voting[seq$start-1+homoLoci3, 2] <<- h1Voting[seq$start-1+homoLoci3, 2] + 1;
    }
    else if (haploID==2){
      h1Voting[seq$start-1+homoLoci1, 2] <<- h1Voting[seq$start-1+homoLoci1, 2] + 1;
      h1Voting[seq$start-1+homoLoci3, 1] <<- h1Voting[seq$start-1+homoLoci3, 1] + 1;
    }

  }


  # converts voting to haplotypes
  # seq = start and endpoint for the haplotype
  # genotype = genotype of the proband
  # allele voting table for the whole genome
  getHaploFromVoting <- function(seq, genotype) {

    haplotypes <- list();
    voteOfInterest <- h1Voting[cSurr$start : (cSurr$end-1), ]
    presentVotes <- (apply(voteOfInterest, 1,sum)>0);
    
    winningVotes <- apply(voteOfInterest, 1,which.max);

    winningVotesH1 <- winningVotes;
    winningVotesH1[winningVotesH1==2] <- 3

    winningVotesH2 <- winningVotes;
    winningVotesH2[winningVotesH2==1] <- 3
    winningVotesH2[winningVotesH2==2] <- 1

    haplotype <- integer(length(genotype))
    probandsHomoLoci <- which((genotype==1)|(genotype==3));
 

    haplotypes[[1]] <- haplotype;
    haplotypes[[2]] <- haplotype;
    
    haplotypes[[1]][presentVotes] <- winningVotesH1[presentVotes];
    haplotypes[[2]][presentVotes] <- winningVotesH2[presentVotes];

    haplotypes[[1]][probandsHomoLoci] <- genotype[probandsHomoLoci]
    haplotypes[[2]][probandsHomoLoci] <- genotype[probandsHomoLoci]
    
    haplotypes
  }


  

  
  countEmpty <- 0; countMatchBoth <- 0; countMatchH1 <- 0; countMatchH2 <- 0; countMatchNone <- 0; countNoOverlap <- 0;
  hapMatching <- data.frame()

  allSurrCount <- dim(allSurrogates)[1]


  #######
  ## LOOP OVER THE SURROGATES, starting with the longest one
  if (allSurrCount==0) {surrIds <- c()} else {surrIds <- 1:allSurrCount}
  for (s in surrIds ) {

    # the sequence currently considered 
    cSurr <- allSurrogates[s,]

    # get other sequences within this one
    #sH1 <- getSeqsWithin(cSurr, selectedSequencesH1);
    #sH2 <- getSeqsWithin(cSurr, selectedSequencesH2);  

       
                                        # phase the haplotypes with the sequences classified up till now
    genotype <- iGene[[individualPhased]][cSurr$start : (cSurr$end-1)]
    haplotypes <- getHaploFromVoting(cSurr, genotype)
    #haplotypes <- doPartialPhasing(individualPhased,cSurr$start, cSurr$end, sH1, sH2,  chromEnds, map, iGene, chromList)  
   
    # key loci - heterozygous loci that had been phased
    keyLoci <- sum((genotype==2) & (haplotypes[[1]]!=0))

    if (keyLoci == 0) {
      selectedSequencesH1<<-rbind(selectedSequencesH1, cSurr)
      addSeqVoting(cSurr, 1)
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
              addSeqVoting(cSurr, 1)
              countMatchH1 <- countMatchH1 + 1;
            }
            else if (matchH2) {
              # matching H2
              # sequence accepted
              selectedSequencesH2<<-rbind(selectedSequencesH2, cSurr)
              addSeqVoting(cSurr, 2)
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

# above are all the functions, which use the data structures defined below
# here the body of the function starts
  
      selectedSequencesH1 <<- data.frame()
      selectedSequencesH2 <<- data.frame()
      rejectedSequences <<- data.frame();
      matchingBothSequences <<- data.frame();
      noOverlapSequences <<- data.frame();

      partialSequencesH1 <<- data.frame();
      partialSequencesH2 <<- data.frame();

      # the borders had been trimmed earlier, no need to do this again 
      surrogateMargin <- 0
      allSurrogates <<- getAllSurrogatesOf(surrMat, individualPhased, surrogateMargin)   

      # at this stage, according to the settings, the children/parents are excluded
      notExcludedParents <- !(allSurrogates$id %in% excludedParents);
      allSurrogates <<- allSurrogates[notExcludedParents,]

      # sort the sequence by their length
      sortedSurrogates <- data.frame()
      if ((dim(allSurrogates)[1])>0) {
        sortedSurrogates <- allSurrogates[order(-allSurrogates$end+allSurrogates$start) , ]
      }

      result <- list();
      result[['summary']] <- assignToHaplotypes(sortedSurrogates)         
      result[['selectedSequencesH1']] <- selectedSequencesH1
      result[['selectedSequencesH2']] <- selectedSequencesH2


      result


    }

