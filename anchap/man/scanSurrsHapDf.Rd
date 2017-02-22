\name{scanSurrsHapDf}
\alias{scanSurrsHapDf}
\title{
Scanning for surrogate parents, from partially phased haplotypes. Stage III of ANCHAP. 
}
\description{
Scans partially phased haplotypes for IBD segments.
}
\usage{
scanSurrsHapDf(result, lengthThresh, chromEnds, noHaps=2*length(result$phasingResults), dataFraction=-1, geneticMap=NULL, gDistT=NULL)
}

\arguments{

  \item{result}{A result of the phasing function. A list.}
  \item{lengthThresh}{IBD threshold, expressed as a count of consecutive markers. Default 200 SNPs.}
  \item{chromEnds}{Identidiers of last SNP on each chromosome.}
  \item{noHaps}{Number of haplotypes to consider in the analysis.}
  \item{dataFraction}{Count of SNPs to analyse. -1 to analyse all of them.    }
  \item{geneticMap=NULL}{An array with each entry corresponding to a SNP. The value gives a genetic position of the SNP according to a genetic map.}
  \item{gDistT=NULL}{IBD threshold, expressed in cM. Default 2 cM.}
}

\value{
  \item{surrMatHapDf}{
	A surrogacy matrix expressed as a data frame. Has fields: start, end, id1, id2.
    }
}
\references{
}
\author{
    Dominik Glodzik
}
\examples{
}
\keyword{ haplotype }
\keyword{ surrogate }
\keyword{ parents }