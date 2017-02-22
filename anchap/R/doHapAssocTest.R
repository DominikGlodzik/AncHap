doHapAssocTest <- function(
                           clusterAssignmentH1,
                           clusterAssignmentH2,
                           clusterSizes,
                           phenotypes,
                           MIN.CLUSTER.SIZE) {

  
  clusterSizes10 <- clusterSizes[clusterSizes>=MIN.CLUSTER.SIZE]
  
  if (length(clusterSizes10)>0) {
    
    topClusters <- as.integer(unlist(dimnames(clusterSizes10)));
    
    haplotypeId1ofN <- matrix(0,length(clusterAssignmentH1), length(topClusters));
    colnames(haplotypeId1ofN ) <-  c(topClusters) 
    for (indiv in 1:(length(clusterAssignmentH1))) {
      
      hapId1 <- clusterAssignmentH1[indiv];    
      if (hapId1 %in% topClusters) {    haplotypeId1ofN[indiv,as.character(hapId1)] <- haplotypeId1ofN[indiv,as.character(hapId1)] + 1; }
      hapId2 <- clusterAssignmentH2[indiv];    
      if (hapId2 %in% topClusters) {    haplotypeId1ofN[indiv,as.character(hapId2)] <- haplotypeId1ofN[indiv,as.character(hapId2)] + 1; }
      
  }
    colnames(haplotypeId1ofN ) <- paste("h_",c(topClusters),sep="")
    
    dataHapPh <- data.frame(haplotypeId1ofN, phenotypes)
    
    xnam <- paste(colnames(haplotypeId1ofN ), sep="")
    fmla <- as.formula(paste("phenotypes ~ ", paste(xnam, collapse= "+")))
    
                                        # fit the linear model
                                        # use the anova table: sum of squares and degrees of freedom
    anovaTest <- anova(lm(fmla,dataHapPh))
    
    ## test the hypothesis that there is no
                                        # f test
    
    noTopHaps <- length(topClusters)

  
    x <- (sum(anovaTest$"Sum Sq"[1:noTopHaps])/sum(anovaTest$"Df"[1:noTopHaps]))/anovaTest$"Mean Sq"[noTopHaps+1]

  
  # get a p-Value
    pValue <- 1 - pf(x, noTopHaps, anovaTest$"Df"[1+noTopHaps])
    
    logP <- -log(pValue)
    
    logP
}
  else {
    print("The association test failed: clusters are too small.")
  }
  
}
