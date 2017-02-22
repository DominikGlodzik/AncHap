\name{analyseTrios}
\alias{analyseTrios}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find trios in the cohort and compare their haplotypes.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analyseTrios(iGene, chromEnds, pedigree, orcaIds, resultPath, phasingsFolder, plotFolder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iGene}{
%%     ~~Describe \code{iGene} here~~
}
  \item{chromEnds}{
%%     ~~Describe \code{chromEnds} here~~
}
  \item{pedigree}{
%%     ~~Describe \code{pedigree} here~~
}
  \item{orcaIds}{
%%     ~~Describe \code{orcaIds} here~~
}
  \item{resultPath}{
%%     ~~Describe \code{resultPath} here~~
}
  \item{phasingsFolder}{
%%     ~~Describe \code{phasingsFolder} here~~
}
  \item{plotFolder}{
%%     ~~Describe \code{plotFolder} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(
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

  
  oppHomo <- countOppHomo(iGene[[iProb]][1:max(chromEnds)],iGene[[iFa]][1:max(chromEnds)])
  inconsistencies[length(inconsistencies)+1] <<- oppHomo
  #shared <- getSurrMatEntry(iProb,iFa)
  #sharedLength <- sum(shared$end - shared$start)/chromEnds[22]

  oppProbFa <-  oppHomo;
  #print(paste("Proband - father",  oppHomo, "IBD sharing", sharedLength ))
  
  oppHomo <- countOppHomo(iGene[[iProb]][1:max(chromEnds)],iGene[[iMo]][1:max(chromEnds)])
  inconsistencies[length(inconsistencies)+1] <<- oppHomo
  #shared <- getSurrMatEntry(iProb,iMo)
  #sharedLength <- sum(shared$end - shared$start)/chromEnds[22]

  oppProbMo <- oppHomo;

  data.frame(oppProbFa=oppProbFa, oppProbMo=oppProbMo)
  }





# make the fields of pedigree accesible

inconsistencies <- vector()

# stores the ids (indices of myMatrix) oof trios
trios<-data.frame();

print("Looping over the pedigree")
for (p in 1:length(id)) {


probId<-pedigree$id[p]
probMo<-pedigree$mo[p]
probFa<-pedigree$fa[p]

  if ((probId \%in\% orcaIds)&&(probMo \%in\%  orcaIds)&&(probFa \%in\%  orcaIds)) {

    iProb <- which(orcaIds==probId)
  
    iFa <- which(orcaIds==probFa)

    iMo <- which(orcaIds==probMo)
    print(paste("trio", iProb, iFa, iMo))
    

    oppHomosParents <- checkTrio(iProb, iFa,  iMo )
    hapComp <- analyseTrioHaplotyped(iProb,  iMo, iFa, phasingsFolder, plotFolder, genomeLength, chromEnds)
    aTrio <- data.frame(prob=probId, fa=probFa, mo=probMo)

    aResult<-cbind(data.frame(probId=probId, faId=probFa, moId = probMo), oppHomosParents, hapComp)
    trios<-rbind(trios, aResult)
  }

  }
print(paste("-- estimated error rate", mean(inconsistencies)/ genomeLength))

# save the results to a file
resultPath<-paste(plotFolder,'trioResults.Rdata')


save(trios,file=resultPath)
trios
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
