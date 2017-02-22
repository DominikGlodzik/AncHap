extractIBDc <- function(id1,id2, g2ID, lThresh, chromEnds,geneticMap=NULL, gDistT=NULL) {
 
# Extracts IBD sequences betweem id1 and id2, which are longer than lThresh.     
#
# Args:
#   id1 : Genotype of the first individual.
#   id2 : Genotype of the second individual.
#   lThresh : The threshold above which an non-opposing homozygote is considered IBD.
#   chromEnds : For each chromosome, indices of last SNPs on every chromosome.
  
# Returns:
#   A data frame with sound IBD sequences between id1 and id2.

  ###### prepare the input for the function call

  # the genotypes of two individuals
  iid1<- as.integer(id1)
  iid2<- as.integer(id2)
  # genotype length
  n <- length(id1)

  # the maximum number of sequences you could have
  if (!is.null(lThresh)) {
  resultC = as.integer(320000/lThresh);
  } else {
    resultC <- 1000
  }
  
  # prepare matrices to store the found sequences
  starts <- integer(resultC);
  ends <- integer(resultC);
  # prepare the count of sequences
  noSeqs <- 0;
  resList <- data.frame();
  ###### end prepare the input for the function call

  ###### call the C function

                                        # if the genetic distance is NOT given
  if (is.null(geneticMap)) {
    
    r <- .C("ibdThreshMiss",n, iid1, iid2, as.integer(lThresh), as.integer(noSeqs), starts, ends, as.integer(chromEnds))

###### format the results of the call
    noSeqs <- r[[5]];
    starts <- r[[6]]
    ends <- r[[7]]
    
  } else {
    # the genetic map has been given

    r <- .C("ibdThreshMiss",n, iid1, iid2, as.integer(lThresh+1), as.integer(noSeqs), starts, ends, as.integer(chromEnds))

###### format the results of the call
    noSeqs <- r[[5]];
    starts <- r[[6]];
    ends <- r[[7]];

    # filter only the sequences longer than a given genetic length threshold
    validSeq <- (geneticMap[ends[1:noSeqs]] - geneticMap[starts[1:noSeqs]])>gDistT;

    starts <- starts[validSeq]
    ends <- ends[validSeq]
    noSeqs <- sum(validSeq)


    
    if (noSeqs>0) {
      invalidSeq <- ((geneticMap[ends[1:noSeqs]] - geneticMap[starts[1:noSeqs]])<gDistT);
      if ((length(invalidSeq)>0)&&(sum(invalidSeq)>0)) {browser()}
    }

    

     
  }
  
  
  # convert the vectors into  data frame
  

  
  if (noSeqs>0) {
    resList <- data.frame(start=starts[1:noSeqs],
                          end = (ends[1:noSeqs]),
                          length=(ends[1:noSeqs]-starts[1:noSeqs]),
                          id = g2ID
                          )
  }
  ###### end format the results of the call

  iid1 <- NULL
  iid2 <- NULL
  starts <- NULL
  ends <- NULL
  #gc()
  

 # return 
 resList
} 
