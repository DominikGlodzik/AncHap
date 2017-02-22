\name{surrMatToDf}
\alias{surrMatToDf}
\title{
Convertion between formats of surrogacy matrices
}
\description{

Converts the surrogacy matrix from a 2D list to a data frame, with fields: start, end, id1, id2.

}
\usage{
surrMatToDf <- function(surrMat)
}
\arguments{
  \item{surrMat}{
	A 2D list, with i,j th entry correpsonding to sharing bewtween between individuals i,j. Each entry is a data frame.
    }
}
\details{
}
\value{
A data frame, with fields: start, end, id1, id2.
}
\references{
}
\author{
Dominik Glodzik
}
\examples{
}

\keyword{ surrogacy}
\keyword{matrix}
\keyword{conversion}