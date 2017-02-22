\name{analyseTrioHaplotyped}
\alias{analyseTrioHaplotyped}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare haplotypes of a parent-offspring trio.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analyseTrioHaplotyped(iProb, iMo, iFa, phasingsFolder, plotFolder, genomeLength, chromEnds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iProb}{
%%     ~~Describe \code{iProb} here~~
}
  \item{iMo}{
%%     ~~Describe \code{iMo} here~~
}
  \item{iFa}{
%%     ~~Describe \code{iFa} here~~
}
  \item{phasingsFolder}{
%%     ~~Describe \code{phasingsFolder} here~~
}
  \item{plotFolder}{
%%     ~~Describe \code{plotFolder} here~~
}
  \item{genomeLength}{
%%     ~~Describe \code{genomeLength} here~~
}
  \item{chromEnds}{
%%     ~~Describe \code{chromEnds} here~~
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
function (iProb, iMo, iFa, phasingsFolder, plotFolder, genomeLength, chromEnds) {

  compareHaplos <- function (id1, id2, iParent, height) {
    
    plotDisagree <- function (pnts, height, c) {
      ys <- integer(length(pnts)) + height
      points(pnts,  ys, lwd=0.05, pch="x", cex=0.2, col=c) 
    }
    
    plotSurrs <- function (surrs, height, c) {
      ys <- c(height, height)
      
      if ((dim(surrs)[1])>0) {
        for (s in 1:(dim(surrs)[1])) {
          xs <- c(surrs[s,'start'], surrs[s,'end'])
          lines(xs,ys, col=c)
          iTaken[surrs[s,'start']:surrs[s,'end']] <<- 1
        }
      }
    
    }


    iTaken <- integer(genomeLength)

    
    SURrTHRESHOLD=100;
    
    surrs11<- extractIBDc(id1[[1]], id2[[1]], iParent, SURrTHRESHOLD, chromEnds)
    h1h1 <- sum(surrs11$end - surrs11$start)/genomeLength
    
    plotSurrs(surrs11, height, 'red')
    
    surrs12<- extractIBDc(id1[[1]], id2[[2]], iParent, SURrTHRESHOLD, chromEnds)
    h1h2 <- sum(surrs12$end - surrs12$start)/genomeLength
    plotSurrs(surrs12, height+1, 'red')
    
    surrs21<- extractIBDc(id1[[2]], id2[[1]], iParent, SURrTHRESHOLD, chromEnds)
    h2h1 <- sum(surrs21$end - surrs21$start)/genomeLength
    plotSurrs(surrs21, height+2, 'red')
  
    surrs22<- extractIBDc(id1[[2]], id2[[2]], iParent, SURrTHRESHOLD, chromEnds)
    h2h2<-sum(surrs22$end - surrs22$start)/genomeLength
    plotSurrs(surrs22, height+3, 'red')
    
    allSurrs <- rbind(surrs11, surrs12, surrs21, surrs22)
    plotSurrs(allSurrs, height+4, 'blue')

    totalHapSharing <- h1h1+h1h2 + h2h1+ h2h2;
  }

  
  samples <- sample(1:genomeLength, 20000,replace=FALSE)
  fileName <- paste(plotFolder, "trio","_", iProb ,'.pdf', sep="")
  pdf(fileName)
  
  plot(c(1,max(chromEnds)), c(0 ,11),'n', ylab='', xlab='', bty='o',  yaxt="n")
  
  prob <- list()
  filename <- paste(phasingsFolder, "results",iProb, ".Rdata", sep="")
  load(filename)
  prob[[1]] <- h1
  prob[[2]] <- h2
  
  probPhased <-  1 -(sum(h1==as.raw(0)))/genomeLength
  

  mo <- list()
  filename <- paste(phasingsFolder, "results",iMo, ".Rdata", sep="");
  load(filename)
  mo[[1]] <- h1
  mo[[2]] <- h2
  
  moPhased <-  1 - (sum(h1==as.raw(0)))/genomeLength
  moSharing <- compareHaplos(prob,mo, iMo, 1)
  
  fa <- list()
  filename <- paste(phasingsFolder, "results",iFa, ".Rdata", sep="")
  load(filename)
  fa[[1]] <- h1
  fa[[2]] <- h2
  faPhased <-  1 - ( sum(h1==as.raw(0)))/genomeLength
  faSharing <- compareHaplos(prob,fa, iFa, 7)
  
  dev.off()
  print(paste("Generated plot ", fileName))
  data.frame( probPhased= probPhased, moPhased=moPhased, faPhased=faPhased, moSharing=moSharing, faSharing=faSharing)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
