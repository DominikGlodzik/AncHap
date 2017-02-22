\name{readGenMap}
\alias{readGenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
aaaa
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcCoverage(sequences, allSnps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
%%     ~~Describe \code{sequences} here~~
}
  \item{allSnps}{
%%     ~~Describe \code{allSnps} here~~
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
function (sequences, allSnps) 
{
    result <- list()
    noSels <- dim(sequences)[1]
    if (noSels > 0) {
        sequences <- sequences[order(sequences$start), ]
        stitches <<- data.frame()
        stitchStart <- sequences[1, "start"]
        stitchEnd <- sequences[1, "end"]
        i <- 1
        while (i <= noSels) {
            currStart <- sequences[i, "start"]
            currEnd <- sequences[i, "end"]
            if ((stitchEnd - 1) < currStart) {
                stitch <- data.frame(start = stitchStart, end = stitchEnd)
                stitches <- rbind(stitches, stitch)
                stitchStart <- currStart
                stitchEnd <- currEnd
            }
            else {
                stitchEnd <- max(stitchEnd, currEnd)
            }
            i <- i + 1
        }
        coverage <- -1
        stitch <- data.frame(start = stitchStart, end = stitchEnd)
        stitches <<- rbind(stitches, stitch)
        coverage <- round(sum(stitches$end - stitches$start)/allSnps, 
            digits = 5)
        data.frame(coverage = coverage, noStitches = dim(stitches)[1])
        result[["coverage"]] <- coverage
        result[["noStitches"]] <- dim(stitches)[1]
        result[["stitches"]] <- stitches
    }
    else {
        result[["coverage"]] <- 0
        result[["noStitches"]] <- 0
        result[["stitches"]] <- 0
    }
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
