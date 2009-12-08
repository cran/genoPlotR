\name{genoPlotR-package}
\alias{genoPlotR-package}
\alias{genoPlotR}
\docType{package}
\title{
  genoPlotR - a R framework to produce publication-grade maps of genes
  and genomes.
}
\description{
  A R framework to plot comparison of gene stretches or genomes, a la
  ACT (Artemis Comparison Tool), but with production-grade graphics, and
  a static interface. Reads directly from tabular files or from wide-spread
  biological formats such as BLAST and PTT (NCBI).
}
\details{
  \tabular{ll}{
    Package: \tab genoPlotR\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2009-12-08\cr
    License: \tab GPL (>=2)\cr
    LazyLoad: \tab yes\cr
  }
  The only plotting function is \code{\link{plot_gene_map}}, which
  produces \code{link[grid]{grid}} graphics. Data is composed mainly of
  DNA segments (\code{\link{dna_seg}}) objects, which represent
  collections of genes or segments of genomes, and of
  \code{\link{comparison}} objects, which are the pairwise comparisons
  between the \code{dna_seg}s. Data can be read from files (see
  \code{\link{read_functions}}) or from R objects like
  \code{data.frame}s or \code{list}s, with \code{\link{dna_seg}} and
  \code{\link{comparison}} conversion functions. 
}
\author{
  Lionel Guy
  
  Maintainer: Lionel Guy \email{lionel.guy@ebc.uu.se}
}
\keyword{ aplot }
\seealso{
  \code{\link{plot_gene_map}} for plotting. \code{\link{dna_seg}} and
  \code{\link{comparison}} for the base objects and conversion functions.
  \code{\link{read_dna_seg_from_tab}}, \code{\link{read_dna_seg_from_ptt}}, 
  \code{\link{read_comparison_from_tab}} and
  \code{\link{read_comparison_from_blast}} to read from files.
}
\examples{
## simple example
## dna segments
## data.frame with several genes
names1 <- c("feat1", "feat2", "feat3")
starts1 <- c(2, 1000, 1050)
ends1 <- c(600, 800, 1345)
strands1 <- c("-", -1, 1)
cols1 <- c("blue", "grey", "red")
df1 <- data.frame(name=names1, start=starts1, end=ends1,
                  strand=strands1, col=cols1)
dna_seg1 <- dna_seg(df1)
is.dna_seg(dna_seg1)

## with only one gene, or two, and merging
gene2a <- dna_seg(list(name="feat1", start=50, end=900, strand="-", col="blue"))
genes2b <- dna_seg(data.frame(name=c("feat2", "feat3"), start=c(800, 1200),
                              end=c(1100, 1322), strand=c("+", 1),
                              col=c("grey", "red")))
dna_seg2 <- c.dna_seg(gene2a, genes2b)
is.dna_seg(dna_seg2)

## reading from file
dna_seg3_file <- system.file('extdata/dna_seg3.tab', package = 'genoPlotR')
dna_seg3 <- read_dna_seg_from_tab(dna_seg3_file)
is.dna_seg(dna_seg3)

## comparison
## from a data.frame
comparison1 <- as.comparison(data.frame(start1=starts1, end1=ends1,
                                        start2=dna_seg2$start,
                                        end2=dna_seg2$end))
is.comparison(comparison1)

## from a file
comparison2_file <- system.file('extdata/comparison2.tab',
                                package = 'genoPlotR')
comparison2 <- read_comparison_from_tab(comparison2_file,
                                        color_scheme="red_blue")
is.comparison(comparison1)

## plot
plot_gene_map(dna_segs=list(dna_seg1, dna_seg2, dna_seg3),
              comparisons=list(comparison1, comparison2)) 
}
