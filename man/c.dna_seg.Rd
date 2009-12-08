\name{c.dna_seg}
\alias{c.dna_seg}
\title{
  Concatenate dna_seg objects
}
\description{
  Concatenate dna_seg objects.
}
\usage{
\method{c}{dna_seg}(...)
}
\arguments{
  \item{\dots}{
    \code{dna_seg}s to be concatenated.
  }
}
\value{
  A \code{dna_seg} object
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{dna_seg}}
}
\examples{
## load data
data(three_genes)

dna_segs[1:2]
c(dna_segs[[1]], dna_segs[[2]])

}
\keyword{ data }
