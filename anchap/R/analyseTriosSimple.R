# gives the error rate 
analyseTriosSimple <- function(
                         allData,
                         H
                         ) {

# length of the region where haplotypes are compared
blockLength <- 20;
noSnps <- dim(H)[2]
noBlocks <- floor(noSnps/blockLength)

pedigree <- allData$pedigree;
orcaIds <- allData$orcaIds[1:floor(dim(H)[1]/2)];

# stores the ids (indices of myMatrix) oof trios
trios<-data.frame();

trioId <- 1; pMoDisc <- vector(); pFaDisc <- vector();
print("Looping over the pedigree")
for (p in 1:length(pedigree$id)) {

probId <- pedigree$id[p]
probMo <- pedigree$mo[p]
probFa <- pedigree$fa[p]

                                        # a trio detected
  if ((probId %in% orcaIds)&&(probMo %in%  orcaIds)&&(probFa %in%  orcaIds)) {

    aTrio <- data.frame(prob=probId, fa=probFa, mo=probMo)
   
    iProb <- which(orcaIds==probId) 
    iFa <- which(orcaIds==probFa)
    iMo <- which(orcaIds==probMo)

    probandHaps <- list(); fatherHaps <- list(); motherHaps <- list();

    # prepare matrices
    # blocks x block position

    probandHaps[[1]] <- t(matrix(H[(iProb-1)*2+1,1:(noBlocks*blockLength)], ncol=noBlocks))
    probandHaps[[2]] <- t(matrix(H[iProb*2,1:(noBlocks*blockLength)], ncol=noBlocks))

    motherHaps[[1]] <- t(matrix(H[(iMo-1)*2+1,1:(noBlocks*blockLength)], ncol=noBlocks))
    motherHaps[[2]] <- t(matrix(H[iMo*2,1:(noBlocks*blockLength)], ncol=noBlocks))

    fatherHaps[[1]] <- t(matrix(H[(iFa-1)*2+1,1:(noBlocks*blockLength)], ncol=noBlocks))
    fatherHaps[[2]] <- t(matrix(H[iFa*2,1:(noBlocks*blockLength)], ncol=noBlocks))
    
    # compare proband, mother
    rMo <- comparePairHaps(probandHaps, motherHaps)
    # compare proband, father
    rFa <- comparePairHaps(probandHaps, fatherHaps)

    aTrioResult<- data.frame(moInconsCount=rMo[[1]], moSignifLociCount=rMo[[2]],faInconsCount=rFa[[1]],  faSignifLociCount=rFa[[2]])
    trios<-rbind(trios, aTrioResult)

    print(trioId)
    trioId <- trioId + 1;
 
  }
}

trioErrorRate <- (sum(trios$moInconsCount) + sum(trios$faInconsCount))/sum(sum(trios$moSignifLociCount) + sum(trios$faSignifLociCount))

result <- list();
result[['trios']] <- trios
result[['trioErrorRate']] <- trioErrorRate
result



}
