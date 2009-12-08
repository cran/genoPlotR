\name{range.dna_seg}
\alias{range.dna_seg}
\alias{range.comparison}
\alias{range.annotation}
\title{
  Range calculation
}
\description{
  Calculate the range of dna_seg and comparisons.
}
\usage{
\method{range}{dna_seg}(x, ...)
\method{range}{comparison}(x, overall=TRUE, ...)
\method{range}{annotation}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Object to calculate the range from.
  }
  \item{overall}{
    Logical, \code{TRUE} by default. Should the range be calculated over
    the whole object? If \code{FALSE}, a range is calculated on each
    side of the comparison.
  }
  \item{\dots}{
    Unused.
  }
}
\details{
  Calculate the overall range of a \code{dna_seg}, \code{comparison} or
  an \code{annotation} object. 
}
\value{
  A numeric of length 2. For \code{comparison}, if \code{overall} is
  \code{FALSE}, a data frame with two rows and two columns, \code{xlim1}
  and \code{xlim2}.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{dna_seg}}, \code{\link{comparison}}, \code{\link{trim}}
  for further examples.
}
\examples{
## Load data
data(three_genes)

## On dna_seg
dna_segs[[1]]
range(dna_segs[[1]])

## On comparison
comparisons[[2]]
range(comparisons[[2]])
range(comparisons[[2]], overall=FALSE)
}
\keyword{ data }
