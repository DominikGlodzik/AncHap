\name{pickReseqIndsFast}
\alias{pickReseqIndsFast}
\title{
Optimisation of re-sequencing studies.


(Represents the coverage matrix as binary variables.)

}
\description{
The uncovered IBD sharing within a cohort can be used for efficient selection of individuals to resequence, with a view to using them as a reference for imputation. Selection of individuals for resequencing is based on maximizing representation of haplotypes and minimizing multiple resequencing of the same haplotypes. As the first individual for resequencing we choose the one whose haplotypes have the most copies in the rest of the cohort. After excluding regions that have been covered by sharing with individuals already chosen, we repeat the procedure of selecting the individual with most copies until a target level of coverage has been achieved.
}

\usage{
pickReseqIndsFast(surrMatDf, noPick, plotPath=NULL, data, appFactor=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surrMatDf}{A data frame with regions of sharing betweeen indviduals.}
  \item{noPick}{Number of individuals for re-sequencing to pick. Cannot be more than there are distinct individuals in the surrogacy matrix. Defaults to all individuals whose SNP genotypes are provided.}
  \item{plotPath}{Path to a plot file (.pdf). NULL (default) if plotting not required.}
  \item{data}{Cohort data in the standard anchap format.}
  \item{appFactor}{Approximation factor. Number of markers in a marker block. Defaults to 10. No approximation if set to 1.}
}

\value{
  \item{pickedInds}{A data frame with a row for every picked individual.
 Has fields: id, cohortUsefulCoverage, cohortCoverage.
Id indicates the preference for resequencing. cohortUsefulCoverage tells what percentage of all genotypes is covered by haplotypes of corresponding individuals, given that all individual of higher priority are sequenced (incldues the genotype of corresponding individual). cohortCoverage for a given individual indicates what percentage of all genotypes in the cohort is spanned by haplotypes shared with the individual (does NOT incldue the genotype of corresponding individual).
}
}
\author{
    Dominik Glodzik
}
\examples{
}
\keyword{ resequencing}
\keyword{optimisation}