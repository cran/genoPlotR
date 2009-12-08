\name{auto_annotate}
\alias{auto_annotate}
\title{
  Auto-annotate dna_segs
}
\description{
  Annotate dna_segs in a smart way. This is especially designed for
  \code{dna_seg}s read from genbank or embl files, but can be extended for
  other uses. In short, it produces annotations from \code{dna_seg}s, grouping
  the tags for operons (atpA, atpB, atC) into one tag (atpA-C), and
  similarly for numbered genes (bep1-9).
}
\usage{
auto_annotate(dna_seg, locus_tag_pattern=NULL, names=dna_seg$gene,
              keep_genes_only=TRUE, ...)
}
\arguments{
  \item{dna_seg}{
    A \code{dna_seg} object.
  }
  \item{locus_tag_pattern}{
    \code{NULL} by default. A character giving a pattern, that is used
    to simplify names. Specially useful to transform long locus tags
    into numbers (e.g. Eco003456 becomes 3456).
  }
  \item{names}{
    A character vector with as many elements as there are rows in the
    \code{dna_seg}. By default, the \code{gene} column of the
    \code{dna_seg} is taken. Gives the names to be summarized.
  }
  \item{keep_genes_only}{
    A logical, \code{TRUE} by default. If set, the genes that have a
    name that is "-" or empty are not annotated.
  }
  \item{\dots}{
    Further arguments to be passed to \code{annotation} function, like
    \code{rot} or \code{color}.
  }
}
\value{
  An \code{annotation} object.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{annotation}}, \code{\link{dna_seg}}.
}
\examples{
## Prepare dna_seg
names <- paste("Eco", sprintf("\%04d", 1:20), sep="")
gene <- c("-", "atpC", "atpB", "atpA", "atp2", 
          "-", "-", "cda1", "cda2", "cda3",
          "vcx23", "vcx22", "vcx21", "cde20",
          "-", "gfrU", "gfrT", "gfrY", "gfrX", "gfrW")
ds <- dna_seg(data.frame(name=names, start=(1:20)*3, end=(1:20)*3+2,
                         strand=rep(1, 20), gene=gene,
                         stringsAsFactors=FALSE))
## Original annotation
annot1 <- annotation(x1=middle(ds), text=ds$gene, rot=30)
## auto_annotate with various options
annot2 <- auto_annotate(ds)
annot3 <- auto_annotate(ds, keep_genes_only=FALSE, rot=45)
annot4 <- auto_annotate(ds, keep_genes_only=FALSE,
                        locus_tag_pattern="Eco", col="red")
## Plot
plot_gene_map(list(ds, ds, ds, ds),
              annotations=list(annot1, annot2, annot3, annot4))

}
\keyword{ data }
