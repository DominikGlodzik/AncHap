analyseTrios <- function(
                         iGene,
                         chromEnds,
                         pedigree,
                         orcaIds,
                         resultPath,
                         phasingsFolder,
                         plotFolder
                         ) {

 genomeLength <- length(iGene[[1]])

# check how many opposing homozygotes there are between a child and its parents
# there shouldn't be any
checkTrio <- function(iProb, iFa, iMo) {

    countOppHomo <- function(id1,id2) {
                                        # extract IBD sequences betweem id1 and id2, which are longer than lThresh 
                                        # use dyn.load("ibdFastExtract.so") to load the objec
      iid1<- as.integer(id1)
      iid2<- as.integer(id2)
      n <- length(id1)
      r <- .C("oppHCounter",n, iid1, iid2, as.integer(0))
      oppHCount <- r[[4]]
      oppHCount
    } 

    
  oppHomo <- countOppHomo(iGene[[iProb]][1:chromEnds[22]],iGene[[iFa]][1:chromEnds[22]])
  inconsistencies[length(inconsistencies)+1] <<- oppHomo

  oppProbFa <-  oppHomo;
  
  oppHomo <- countOppHomo(iGene[[iProb]][1:chromEnds[22]],iGene[[iMo]][1:chromEnds[22]])
  inconsistencies[length(inconsistencies)+1] <<- oppHomo


  oppProbMo <- oppHomo;

  data.frame(oppProbFa=oppProbFa, oppProbMo=oppProbMo)
}
#

# make the fields of pedigree accesible

inconsistencies <- vector()

# stores the ids (indices of myMatrix) oof trios
trios<-data.frame();

 
print("Looping over the pedigree")
for (p in 1:length(pedigree$id)) {


probId<-pedigree$id[p]
probMo<-pedigree$mo[p]
probFa<-pedigree$fa[p]

                                        # a trio detected
  if ((probId %in% orcaIds)&&(probMo %in%  orcaIds)&&(probFa %in%  orcaIds)) {

    aTrio <- data.frame(prob=probId, fa=probFa, mo=probMo)
   
    iProb <- which(orcaIds==probId)
  
    iFa <- which(orcaIds==probFa)

    iMo <- which(orcaIds==probMo)
    print(paste("trio", iProb, iFa, iMo))

    oppHomosParents <- checkTrio(iProb, iFa,  iMo )

    r <-  analyseTrioHaplotyped(iProb,  iMo, iFa, phasingsFolder, plotFolder, genomeLength, chromEnds);

 
    aResult<-cbind(data.frame(probId=probId, faId=probFa, moId = probMo), r)
    print(aResult)
    trios<-rbind(trios, aResult)

  }

}
print("-- estimated error rate")
 estErrorRate<- sum(trios$errorCount)/sum(trios$errorProneLoci)
print(estErrorRate)
 
# save the results to a file
resultPath<-paste(plotFolder, 'trioResults.Rdata')


 save(trios, file=resultPath)
trios
}
