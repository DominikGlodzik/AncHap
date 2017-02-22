\name{getAllSurrogatesOf}
\alias{getAllSurrogatesOf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ggggggggggg
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractIBDc(id1, id2, g2ID, lThresh, chromEnds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id1}{
%%     ~~Describe \code{id1} here~~
}
  \item{id2}{
%%     ~~Describe \code{id2} here~~
}
  \item{g2ID}{
%%     ~~Describe \code{g2ID} here~~
}
  \item{lThresh}{
%%     ~~Describe \code{lThresh} here~~
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

# finds all ibd sequences for an individual, using the surrogacy matrix
getAllSurrogatesOf <- function (surrMat, proband, margin) {

# recovers ibs sequences between two individuals, from the surrogacy matrix
getIbs <- function(i1, i2) {if (i1<i2) {ibsSeqs <- surrMat[[i1,i2]]}  else if (i1>i2) {ibsSeqs <- surrMat[[i2, i1]]}; ibsSeqs}
  
  # requires surrMat
  # requires
  # requires function filterSurrogates
  noInds <- dim(surrMat)[1]
  allSurrogates <- data.frame();
  noInds <- dim(surrMat)[1]
  for (i in 1:noInds) {
   if(i!=proband) { 
    surrogates<-getIbs(i,proband)

    if (!is.null(surrogates)) {
    
      noSurrogates <- dim(surrogates)[1]
      if (noSurrogates>0) {
        batch <- data.frame(id=i, start=surrogates$start, end=surrogates$end)
                                        #print(batch)
        allSurrogates <- rbind(allSurrogates,batch)
      }
    }

    
   }   
  }

  if ((dim(allSurrogates)[1])>0) {
  allSurrogates<-allSurrogates[allSurrogates$end>0,]
  is.vector(allSurrogates$start)
  allSurrogates.sorted <- allSurrogates[order(allSurrogates$start) , ]
  allSurrogates <-allSurrogates.sorted

  allSurrogates$start <- allSurrogates$start+margin
  allSurrogates$end <- allSurrogates$end-margin
}
  allSurrogates
}





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
