\name{pickReseqIndsInteger}
\alias{pickReseqIndsInteger}
\title{
Optimisation of re-sequencing studies.
Represents the coverage matrix as a matrix of integers.
}
\description{Chooses individuals for re-sequencing according to the surrogacy matrix.}
\usage{
function (surrMatDf, noPick=noInds, plotPath=NULL, data, noInds=NULL, appFactor=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surrMatDf}{A data frame with regions of sharing betweeen indviduals.}
  \item{noPick}{ Number of individuals for re-sequencing to pick.}
  \item{plotPath}{Path to a plot file (.pdf).}
  \item{data}{ Cohort data in the standard anchap format.}
  \item{noInds}{ Number of individuals to analyse.}
  \item{apFactor}{Approximation factor. Number of markers in a marker block.}
}

\value{
  \item{pickedInds}{A data frame with a row for every picked individual. Has fields: id, cohortUsefulCoverage, cohortCoverage.}
}
\author{
    Dominik Glodzik
}
\examples{
}
\keyword{ resequencing}
\keyword{optimisation}