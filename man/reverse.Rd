\name{reverse}
\alias{reverse}
\alias{reverse.default}
\alias{reverse.dna_seg}
\alias{reverse.comparison}
\title{
  Reverse objects
}
\description{
  Reverse objects, mainly \code{dna_seg} and \code{comparison}
}
\usage{
reverse(x, ...)
\method{reverse}{default}(x, ...)
\method{reverse}{dna_seg}(x, ...)
\method{reverse}{comparison}(x, side = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The object to reverse.
  }
  \item{\dots}{
    Unused.
  }
  \item{side}{
    In the case of comparisons, the side of the comparison that should
    be reversed. If \code{side=1}, the first side will be reversed. If
    \code{side=2}, the second side will be reversed. If \code{side<1},
    no side is reversed. If \code{side>2}, both sides are reversed.
  }
}
\value{
  The same object as input.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{dna_seg}}, \code{\link{comparison}}
}
\examples{
## load data
data(three_genes)

## on dna_seg
dna_segs[[1]]
reverse(dna_segs[[1]])
## on comparison
reverse(comparisons[[2]], side=1)
reverse(comparisons[[2]], side=3)

## With mauve backbone
data(mauve_bbone)
## Plot
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons)

## Reverse B_bacilliformis, and the corresponding comparison (first "side")
bbone$dna_segs[[1]] <- reverse(bbone$dna_segs[[1]])
bbone$comparisons[[1]] <- reverse(bbone$comparisons[[1]], 1)
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons)

}
\keyword{ data }
