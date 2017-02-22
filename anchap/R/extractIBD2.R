extractIBD2 <- function(id1,id2, g2ID, lThresh, chromEnds) {
 
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
  resultC = as.integer(320000/lThresh);
  # prepare matrices to store the found sequences
  starts <- integer(resultC);
  ends <- integer(resultC);
  # prepare the count of sequences
  noSeqs <- 0;
  resList <- data.frame();
  ###### end prepare the input for the function call

  ###### call the C function
  r <- .C("ibd2",n, iid1, iid2, as.integer(lThresh), as.integer(noSeqs), starts, ends, as.integer(chromEnds))

  
  ###### format the results of the call
  noSeqs <- r[[5]];
  starts <- r[[6]]
  ends <- r[[7]]

  if (noSeqs>0) {
    for (i in 1:noSeqs) {
      aResult <- data.frame(start=starts[[i]],end=ends[[i]]-1, length= ends[[i]]-starts[[i]], id= g2ID);
      resList <- rbind(resList, aResult);
    }
  }
  ###### end format the results of the call

 # return 
 resList
} 
