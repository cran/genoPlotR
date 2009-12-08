\name{mauve_bbone}
\docType{data}
\alias{mauve_bbone}
\title{Mauve backbone of 4 Bartonella genomes}
\description{
  The result of a multiple genome alignment with Mauve. 
}
\usage{
data(mauve_bbone)
}
\format{
  \code{bbone}, a list of two dataframes, representing the regions which
  are conserved in at least two genomes:
  \itemize{
    \item{\code{dna_segs}}{which is a list of 4 \code{dna_seg}
      objects, containing the mauve blocks for each genome.}
    \item{\code{comparisons}}{which is a list of 3 \code{comparison}
      objects.}
  }
  A bash script to obtain the same file as in the data is
  available in the \code{extdata} folder of the package. Find its
  location by running 
  \code{system.file('extdata/mauve.sh', package = 'genoPlotR')}.
    
  The resulting backone file can then be read with
  \code{\link{read_mauve_backbone}}.
}
\references{
  Mauve: http://asap.ahabs.wisc.edu/mauve/
}

\examples{
data(mauve_bbone)
plot_gene_map(bbone$dna_segs, bbone$comparisons)
}
\keyword{datasets}
