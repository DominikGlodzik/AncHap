evalIBDvsReference <- function(testSharing, testNames, refSharing, idArray, noSnps){
  # testSharing is a list of data frames
  # refSharing is a data frame

  dfToVector <- function(maxM, df) {
    v <- integer(maxM)
    if ((dim(df)[1])>0) {
      for (r in 1:(dim(df)[1])) { v[df$start[r]:df$end[r]] <- 1; }
    }
    v
  }

  TP <- integer(length(testSharing));
  TN <- integer(length(testSharing));
  FP <- integer(length(testSharing));
  FN <- integer(length(testSharing));
  
  
  # for all individuals with parents seq'd
  for (i in 1:length(idArray)) {
    print(paste("individual",i))
    label1 <- idArray[i]
                                        #for all individuals with parents seq'd
    for (j in 1:length(idArray)) {
      if (i<j) {
        label2 <-idArray[j]

        # extract reference sharing
        
        refRows <- (refSharing$label1==label1 & refSharing$label2==label2) | (refSharing$label1==label2 & refSharing$label2==label1)
        refPairSharing <- refSharing[refRows,]
        refSharingVector <- dfToVector(noSnps, refPairSharing)
        
        # extract test sharing
        # testSharing is a list
        for (tc in 1:length(testSharing)) {
          testPairRows <- (testSharing[[tc]]$label1==label1 & testSharing[[tc]]$label2==label2) | (testSharing[[tc]]$label1==label2 & testSharing[[tc]]$label2==label1)
          testPairSharing <- testSharing[[tc]][testPairRows,]
          testSharingVector <- dfToVector(noSnps, testPairSharing)
          
          TP[tc] <- TP[tc] + sum(refSharingVector & testSharingVector);
          TN[tc] <- TN[tc] + sum(!refSharingVector & !testSharingVector);
          FP[tc] <- FP[tc] + sum(!refSharingVector & testSharingVector);
          FN[tc] <- FN[tc] + sum(refSharingVector & !testSharingVector);  
        }       
        
      } # end if
    } # end inner for
  } # end outer for
  
  result <- data.frame(testNames=testNames, sensitivity=TP/(TP+FN), falsDiscRate=FP/(FP+TP), TP=TP, TN=TN, FP=FP, FN=FN)
  result
  
} # end function

getPairHaploSharing <- function(id1, id2, maxSnp) {
  result <- integer(maxSnp)
  pairRows <- (((haploSharing$id1==id1) & (haploSharing$id2==id2))|((haploSharing$id1==id2) & (haploSharing$id2==id1)))
  haploSharingPair <-  haploSharing[pairRows,]
  if (sum(pairRows)>0) {
    for (i in 1:sum(pairRows)) {
    region <- haploSharingPair[i,'start']:(haploSharingPair[i,'end']-1)
    result[region] <- result[region] + 1;
  }
  }
  result
}

getPairAnchap2Sharing <- function(id1, id2, maxSnp) {
  result <- integer(maxSnp)
  pairRows <- (((surrMatHapDfCut$id1==id1) & (surrMatHapDfCut$id2==id2))|((surrMatHapDfCut$id1==id2) & (surrMatHapDfCut$id2==id1)))
  haploSharingPair <-  surrMatHapDfCut[pairRows,]
  if (sum(pairRows)>0) {
    for (i in 1:sum(pairRows)) {
    region <- haploSharingPair[i,'start']:(haploSharingPair[i,'end']-1)
    result[region] <- result[region] + 1;
  }
  }
  result
}

computeSurrDensity2IBD<- function(surrMatDf, noSnps=NA) {
  
  if (is.na(noSnps)) {noSnps<-max(surrMatDf$end)}
  idMap <- integer(597)
  t <- as.integer(names(table(c(surrMatDf$id1,surrMatDf$id2))))
  # mapping between the identifiers and array indices
  idMap[t] <- 1:length(t)
  # 3d matrix: no fully phased individuals, no fully phased individuals x no Snps
  surrCount <- array(0, dim=c(length(t),length(t),noSnps))
  for (s in 1:(dim(surrMatDf)[1])) {
    surrCount[idMap[surrMatDf$id1[s]], idMap[surrMatDf$id2[s]], surrMatDf$start[s]:(surrMatDf$end[s]-1)] <-   surrCount[idMap[surrMatDf$id1[s]], idMap[surrMatDf$id2[s]], surrMatDf$start[s]:(surrMatDf$end[s]-1)] + 1
  }
  # make independent from IBD2
  surrCount <- (surrCount>0)
  # sum over all pairs of indinviduals
  # colSums removes the first dimension
  # so overall removes the first and the second dimension
  result <- colSums(colSums(surrCount))
  result
}

# blue - always less than the other one
computeSurrDensity2IBDLoop<- function(surrMatDf, noSnps=NA) {
  if (is.na(noSnps)) {noSnps<-max(surrMatDf$end)}
  
  # 1d matrix: no fully phased individuals, no fully phased individuals x no Snps
  surrCount <- integer(noSnps)
  indices <- 1:max(c(surrMatDf$id1,surrMatDf$id2))


  t <- as.integer(names(table(c(surrMatDf$id1,surrMatDf$id2))))

  jointCount <- integer(noSnps)
  for (i in indices) {
    cat(paste(i," "))
    for (j in indices) {
        if ((i<j)) {
        #if ((i<j) && (i %in% t) && (j %in% t)) {
          surrCount <- integer(noSnps)
          pairSurrMatDf <- surrMatDf[(surrMatDf$id1==i) & ((surrMatDf$id2==j)), ]
          if ((dim(pairSurrMatDf)[1])>0) {
            
            for (s in 1:(dim(pairSurrMatDf)[1])) {
              surrCount[pairSurrMatDf$start[s]:(pairSurrMatDf$end[s]-1)] <-  surrCount[pairSurrMatDf$start[s]:(pairSurrMatDf$end[s]-1)] + 1
            } # end loop over pair sequences
            # discard IBD2: only binary IBD status
            surrCount <- (surrCount>0)
            jointCount <- jointCount + surrCount
          } # end if
        } # end if
      } # end inner for
  } # end outer for
  cat("\n")

  result <- jointCount
  result
}
