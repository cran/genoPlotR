\name{chrY_subseg}
\docType{data}
\alias{chrY_subseg}
\title{Comparisons of subsegments of the Y chromosome in human and chimp}
\description{
  A subsegment of the Y chromosome in Homo sapiens and Pan troglodytes,
  to illustrate support for exons and introns.
}
\usage{
data(chrY_subseg)
}
\format{
  A list of two data frames, representing the Y segment in the two
  species, and containing:
  \itemize{
    \item{\code{dna_segs}}{which is a list of two \code{dna_seg}
      objects, containing each three rows (or genes).}
    \item{\code{comparison}}{which is a list of one \code{comparison}
      objects.}
  }
}
\details{
Header for the Homo sapiens genbank file:
LOCUS       NC_000023             220001 bp    DNA     linear   CON 10-JUN-2009
DEFINITION  Homo sapiens chromosome X, GRCh37 primary reference assembly.
ACCESSION   NC_000023 REGION: 2600000..2820000 GPC_000000047

Header for the Pan troglodytes file:
LOCUS       NC_006491             220001 bp    DNA     linear   CON 18-SEP-2006
DEFINITION  Pan troglodytes chromosome X, reference assembly (based on
            Pan_troglodytes-2.1).
ACCESSION   NC_006491 REGION: 2620000..2840000

}
\examples{
data(chrY_subseg)
plot_gene_map(chrY_subseg$dna_segs, chrY_subseg$comparison, dna_seg_scale=TRUE,
              scale=FALSE)
}
\keyword{datasets}
