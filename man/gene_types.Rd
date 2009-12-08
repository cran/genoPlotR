\name{gene_types}
\alias{gene_types}
\title{
  Gene types
}
\description{
  Returns a vector containing the available gene types. In addition to
  these gene types, the user can provide graphical functions that return
  a list or a single \code{grob} object. 
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
\details{
  \code{dna_seg}s may contain one character column
  \code{gene_type}. Elements in this column should either be one of the
  predefined gene types, or refer to a graphical function that has
  exactly the same name and that returns a \code{grob} or a \code{gList}
  object.

  A gene object (i.e. a single row of a \code{dna_seg}) is passed to the
  graphical function, as well as the contents of the \dots. The start
  and line width of an element can thus be accessed via
  \code{gene$start} and \code{gene$lwd}. Extra columns that would be
  added in the \code{dna_seg} can be used similarly. Extra arguments can
  also be globally passed via \dots when calling \code{plot_gene_map}.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{plot_gene_map}}, \code{\link{dna_seg}} 
}
\examples{
## To view pre-coded gene types:
gene_types()

## Load data
data(barto)
n <- length(gene_types(auto=FALSE))

## Get a small subset from the barto dataset
dna_seg <- barto$dna_segs[[3]][1:n,]
plot_gene_map(list(dna_seg))

## Change gene_types and plot again
dna_seg$gene_type <- gene_types(auto=FALSE)
dna_seg$col <- rainbow(n)
dna_seg_r <- dna_seg
dna_seg_r$strand <- -dna_seg$strand

## Add an annotation
annot <- annotation(middle(dna_seg), text=dna_seg$gene_type, rot=45,
                    col=dna_seg$col)

## Plot
plot_gene_map(list(dna_seg, dna_seg_r), annotations=list(annot, annot),
              annotation_height=5, dna_seg_line=grey(0.7))


## Using home-made graphical functions
## Data
data(three_genes)

## Functions returning grobs.
## Creates a triangle
triangleGrob <- function(gene, ...) {
  x <- c(gene$start, (gene$start+gene$end)/2, gene$end)
  y1 <- 0.5 + 0.4*gene$strand
  y <- c(y1, 0.5, y1)
  polygonGrob(x, y, gp=gpar(fill=gene$col, col=gene$col, lty=gene$lty,
                      lwd=gene$lwd), default.units="native")
}

## Draws a star. Note that the limits of the dna_seg region are
## voluntarily not respected
starGrob <- function(gene, ...){
  ## Coordinates for the star
  x <- sin(((0:5)/2.5)*pi)*(gene$end-gene$start)/2 + (gene$end+gene$start)/2
  y <- cos(((0:5)/2.5)*pi)*gene$strand*2 + 0.5
  idx <- c(1, 3, 5, 2, 4, 1)
  ## Attribute line_col only if present in the gene
  line_col <- if (!is.null(gene$line_col)) gene$line_col else gene$col
  ## Having a conditional transparency, depending on a length cut-off
  ## passed via dots
  length_cutoff <- list(...)$length_cutoff
  if (!is.null(length_cutoff)){
    alpha <- if ((gene$end-gene$start) < length_cutoff)  0.3 else  0.8
  } else alpha <- 1
  
  ## Grobs
  g <- polygonGrob(x[idx], y[idx], gp=gpar(fill=gene$col, col=line_col,
                                     lty=gene$lty, lwd=gene$lwd, alpha=alpha),
                   default.units="native")
  t <- textGrob(label="***", x=(gene$end+gene$start)/2, y=0.5,
                default.units="native")
  gList(g, t)
}

## Replacing the standard types
dna_segs[[1]]$gene_type <- "triangleGrob"
dna_segs[[2]]$gene_type <- "starGrob"
## Adding more variables
dna_segs[[2]]$line_col <- c("black", grey(0.3), "blue")
## Mix of several types on the same line
dna_segs[[3]]$gene_type <- c("starGrob", "triangleGrob", "arrows")

## Plot
plot_gene_map(dna_segs, comparisons, length_cutoff=600)


}
\keyword{ data }
