\name{extractIBD2}
\alias{extractIBD2}
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
function (id1, id2, g2ID, lThresh, chromEnds) 
{
    iid1 <- as.integer(id1)
    iid2 <- as.integer(id2)
    n <- length(id1)
    resultC = as.integer(320000/lThresh)
    starts <- integer(resultC)
    ends <- integer(resultC)
    noSeqs <- 0
    resList <- data.frame()
    r <- .C("ibdThresh", n, iid1, iid2, as.integer(lThresh), 
        as.integer(noSeqs), starts, ends, as.integer(chromEnds))
    noSeqs <- r[[5]]
    starts <- r[[6]]
    ends <- r[[7]]
    if (noSeqs > 0) {
        for (i in 1:noSeqs) {
            aResult <- data.frame(start = starts[[i]], end = ends[[i]] - 
                1, length = ends[[i]] - starts[[i]], id = g2ID)
            resList <- rbind(resList, aResult)
        }
    }
    resList
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
