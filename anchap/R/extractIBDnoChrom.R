extractIBDnoChrom <- function(id1,id2, lThresh) {
  # extract IBD sequences betweem id1 and id2, which are longer than lThresh 
  
 
  iid1<- as.integer(id1)
  iid2<- as.integer(id2)

  n <- length(id1)
  
  resultC = as.integer(length(id1)/lThresh);
  
  starts <- integer(resultC);
  ends <- integer(resultC);
  noSeqs <- 0;
  resList <- data.frame();
  
  r <- .C("ibdFastExtract",n, iid1, iid2, as.integer(lThresh), as.integer(noSeqs), starts, ends)

  
  noSeqs <- r[[5]];
  starts <- r[[6]]
  ends <- r[[7]]

  if (noSeqs>0) {
  for (i in 1:noSeqs) {
    aResult <- data.frame(start=starts[[i]],end=ends[[i]]);
    resList <- rbind(resList, aResult);
  }
}
 resList
  
} 
