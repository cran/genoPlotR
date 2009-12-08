\name{barto}
\docType{data}
\alias{barto}
\title{Comparison of 4 Bartonella genomes}
\description{
  Comparison of 4 Bartonella genomes by BLAST. 
}
\usage{
data(barto)
}
\format{
  \code{barto}, a list of three dataframes, representing the four genomes
  and their pairwise comparisons:
  \itemize{
    \item{\code{dna_segs}}{which is a list of 4 \code{dna_seg}
      objects, containing all the protein genes for each genome. Obtained
      by reading ptt files downloaded from NCBI with
      \code{read_dna_seg_from_ptt}}.
    \item{\code{comparisons}}{which is a list of 3 \code{comparison}
      objects, obtained by doing genome-to-genome (fasta files) BLASTS,
      and then reading the resulting tab files with
      \code{read_comparison_from_blast}.}
    \item{\code{rnt_segs}}{which is a list of 4 \code{dna_seg} objects,
      containing all the RNA genes of the four genomes. Obtained by
      reading rnt files downloaded from NCBI with
      \code{read_dna_seg_from_ptt}}.
  }
  A bash script to obtain the same file as in the dataset is
  available in the \code{extdata} folder of the package. Find its
  location by running 
  \code{system.file('extdata/barto.sh', package = 'genoPlotR')}.
}
\references{
  BLAST: http://www.ncbi.nlm.nih.gov/blast/
}
\examples{
data(barto)
plot_gene_map(barto$rnt_segs, barto$comparisons, gene_type="blocks")
}
\keyword{datasets}
