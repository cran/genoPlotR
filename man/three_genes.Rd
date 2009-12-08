\name{three_genes}
\docType{data}
\alias{three_genes}
\title{Three genes data set}
\description{
  A set of three made-up genes, compared in three chromosomes.
}
\usage{
data(three_genes)
}
\format{
  Two dataframes, representing the three genes in three DNA
  segments:
  \itemize{
    \item{\code{dna_segs}}{which is a list of three \code{dna_seg}
      objects, containing each three rows (or genes).}
    \item{\code{comparisons}}{which is a list of two \code{comparison}
      objects.}
  }
}
\examples{
data(three_genes)
plot_gene_map(dna_segs, comparisons)
}
\keyword{datasets}
