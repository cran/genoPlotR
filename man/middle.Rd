\name{middle}
\alias{middle}
\title{
  Middles of a dna_seg
}
\description{
  Returns a vector containing the middle of the genes of a
  dna_seg. Useful to prepare annotations, for example.
}
\usage{
middle(dna_seg)
}
\arguments{
  \item{dna_seg}{
    A \code{dna_seg} object.
  }
}
\value{
  A numeric vector.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{annotation}}, \code{\link{dna_seg}}
}
\examples{

## Load data
data(barto)

## Get middles of the first dna_seg
mid <- middle(barto$dna_segs[[1]])

}
\keyword{ data }
