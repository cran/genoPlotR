\name{gene_types}
\alias{gene_types}
\title{
  Gene types
}
\description{
  Returns a vector containing the available gene types.
}
\usage{
gene_types(auto = TRUE)
}
\arguments{
  \item{auto}{
    Logical. Should type "auto" be added?
  }
}
\value{
  A character vector.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{plot_gene_map}}
}
\examples{
gene_types()

## Load data
data(barto)
n <- length(gene_types(auto=FALSE))

## Get a small subset from the barto dataset
dna_seg <- barto$dna_segs[[3]][1:n,]
plot_gene_map(list(dna_seg))

# Change gene_types and plot again
dna_seg$gene_type <- gene_types(auto=FALSE)
dna_seg$col <- rainbow(n)
dna_seg_r <- dna_seg
dna_seg_r$strand <- -dna_seg$strand
plot_gene_map(list(dna_seg, dna_seg_r))

}
\keyword{ data }
