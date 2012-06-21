\name{annotation}
\alias{annotation}
\alias{as.annotation}
\alias{is.annotation}
\title{
  Annotation class and class functions
}
\description{
  An annotation describes a DNA segment. It has labels attached to
  positions. Each label can be attached to a single position or to a
  range. 
}
\usage{
annotation(x1, x2 = NA, text, rot = 0, col = "black")
as.annotation(df, x2 = NA, rot = 0, col = "black")
is.annotation(annotation)
}
\arguments{
  \item{x1}{
    Numeric. A vector giving the first or only position of the
    label. Mandatory. 
  }
  \item{x2}{
    Numeric. A vector of the same length as \code{x1}. If a row (or the
    whole column is  \code{NA}, then the annotation(s) will be attached
    to \code{x0}. Else, the annotation will be attached to the range
    between both positions. \code{NA} by default.
  }
  \item{text}{
    Character of the same length as \code{x0}. Gives the text of the
    labels. Mandatory.
  }
  \item{rot}{
    Numeric of the same length as \code{x0}. Gives the rotation, in
    degrees, of the labels. 0 by default.
  }
  \item{col}{
    Vector of the same length as \code{x0}. The color of the
    labels. \code{black} by default.
  }
  \item{df}{
    A data frame to convert to an annotation object. Should have at
    least columns \code{x1} and \code{text}.
  }
  \item{annotation}{
    An object to test.
  }
}
\details{
  An \code{annotation} object is a data frame with columns \code{x0},
  \code{x1}, \code{text}, \code{col} and \code{rot}. They give,
  respectively, the first (or only) position, eventually the second
  position, the text, the color and the rotation of the annotation. When
  plotted with \code{plot_gene_map}, it will add an annotation row on
  top of the first \code{dna_seg}. Labels for which only one position is
  given will be centered on that position. Labels for which two
  positions are given are linked by an horizontal square bracket and the
  label is plotted in the middle of the positions.
}
\value{
  \code{annotation} and  \code{as.annotation} return an annotation object.
  \code{is.annotation} returns a logical.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{plot_gene_map}}, \code{\link{middle}}.
}
\examples{
## loading data
data(three_genes)

## Calculating middle positions
mid_pos <- middle(dna_segs[[1]])

# Create first annotation
annot1 <- annotation(x1=mid_pos, text=dna_segs[[1]]$name)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, annotations=annot1)

## Exploring options
annot2 <- annotation(x1=c(mid_pos[1], dna_segs[[1]]$end[2]),
                     x2=c(NA, dna_segs[[1]]$end[3]),
                     text=c(dna_segs[[1]]$name[1], "region1"),
                     rot=c(30, 0), col=c("grey", "black"))
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot2, annotation_height=1.3)

## Annotations on all the segments
annots <- lapply(dna_segs, function(x){
  mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=30)
})
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annots, annotation_height=1.8, annotation_cex=1)

##
## Using a bigger dataset from a 4-genome comparison
##
data(barto)
## Adding a tree
tree <- newick2phylog("(BB:2.5,(BG:1.8,(BH:1,BQ:0.8):1.9):3);")
## Showing several subsegments
xlims2 <- list(c(1445000, 1415000, 1380000, 1412000),
               c(  10000,   45000,   50000,   83000, 90000, 120000),
               c(  15000,   36000,   90000,  120000, 74000,  98000),
               c(   5000,    82000))
## Adding annotations for all genomes, allow segments to be placed out
## of the longest segment
annots <- lapply(barto$dna_segs, function(x){
  mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=30)
  # removing gene names starting with "B" and keeping 1 in 4
  idx <- grep("^[^B]", annot$text, perl=TRUE)
  annot[idx[idx \%\% 4 == 0],] 
})
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree,
              annotations=annots,
              xlims=xlims2,
              limit_to_longest_dna_seg=FALSE,
              dna_seg_scale=TRUE)


}
\keyword{ data }
