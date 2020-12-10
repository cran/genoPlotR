\name{plot_gene_map}
\alias{plot_gene_map}
\title{
  Plot gene and genome maps
}
\description{
  This plotting function represents linearly DNA segments and their
  comparisons. It will plot one line per DNA segment, eventually
  separated by the comparisons. In addition, a tree can be plotted on
  the left of the plot, and annotations on the top row. Since this is a
  grid plot, it can be placed into other graphics, or modified subsequently.
}
\usage{
plot_gene_map(dna_segs,
              comparisons = NULL,
              tree = NULL,
              tree_width = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL, 
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = !scale,
              n_scale_ticks=7,
              scale_cex=0.6,
              global_color_scheme = c("auto", "auto", "blue_red", 0.5),
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0,
              ...)
}
\arguments{
  \item{dna_segs}{
    A list of \code{dna_seg} objects. Mandatory.
  }
  \item{comparisons}{
    A list of \code{comparison} objects. Optional. 
  }
  \item{tree}{
    A tree, under the form of a \code{\link[ade4]{phylog}} object. If
    specified, takes place at the left of the tags. See details below
    for more information.
  }
  \item{tree_width}{
    Numeric. The width of the tree area in the plot, in inches. By
    default, takes 20 percent of the total plot.
  }
  \item{tree_branch_labels_cex}{
    Numeric or \code{NULL} (default). If the tree provided contains
    node annotations, they will be displayed with this \code{cex}. If
    equal to 0, node annotations are not displayed.
  }
  \item{tree_scale}{
    Logical. Plot a scale for the tree? Default is FALSE.
  }
  \item{legend}{
    Yet unimplemented.
  }
  \item{annotations}{
    An \code{annotation} object or a list of annotation objects. See
    details. Optional.
  }
  \item{annotation_height}{
    Numeric. The height, in lines, of the annotation line. One by
    default, if \code{annotation is defined}.
  }
  \item{annotation_cex}{
    Numeric. The \code{cex} (i.e. the character expansion) of the
    annotation line. 
  }
  \item{seg_plots}{
    A list of \code{seg_plot} objects of the length as \code{dna_segs},
    a single \code{seg_plot} or \code{NULL} (default). To draw plots
    associated to a \code{dna_seg}. See \code{\link{seg_plot}} for more
    information and some examples.
  }
  \item{seg_plot_height}{
    The height of the \code{seg_plot} regions. By default, equals to 3
    (lines, see next argument).
  }
  \item{seg_plot_height_unit}{
    The unit of the height of the \code{seg_plot} regions. Should be a
    valid \code{unit}. See the \code{grid} documentation for more
    information. If equals to "null", then the height will be calculated
    as a proportion of the comparison region (i.e. 0.5 means the
    \code{seg_plot} region will be half the size of the \code{comparison}. 
  }
  \item{seg_plot_yaxis}{
    Can be \code{NULL}, \code{FALSE} or a numercic. In the first two
    cases, no y-axis is drawn for the \code{seg_plot}s. If numeric, a
    axis is drawn with approximately that number of ticks. 
  }
  \item{seg_plot_yaxis_cex}{
    The character expansion of the \code{seg_plot_yaxis} labels. Equals
    to \code{scale_cex} by default.
  }
  \item{xlims}{
    A list with as many elements as there are \code{dna_seg}s, or
    \code{NULL}. If \code{NULL}, the whole segment will be
    represented. If a list, each element of the list is a numeric
    vector, representing pairs of  left and right limits for each
    subsegment. See details. 
  }
  \item{offsets}{
    A list or a vector with as many elements as there are
    \code{dna_seg}s, or \code{NULL}. If is a numeric vector, gives the
    offset of the first subsegment. If is a list, each element should
    have the same length as there are subsegments (see
    \code{xlims}). Gives then the length of each gap. If \code{NULL},
    the size of the gaps is optimized to minimize the lengths of the
    comparisons. See details.
  }
  \item{minimum_gap_size}{
    A numeric. How much of the plotting region should a gap be, at
    least. Default is 0.05 (20\% of the plotting region).
  }
  \item{fixed_gap_length}{
    Should the gaps have a fixed length? Otherwise, the gap length will
    be optimized to minimize the size of comparisons. \code{FALSE} by
    default. 
  }
  \item{limit_to_longest_dna_seg}{
    A logical. Should the plot be restricted to the longest
    \code{dna_seg}? If no, the other segments can be extended to better
    fit comparisons.
  }
  \item{main}{
    A character. Main title of the plot.
  }
  \item{main_pos}{
    Position of the main title. One of \code{centre}, \code{left} or
    \code{right}.
  }
  \item{dna_seg_labels}{
    A character, same length as \code{dna_segs}. The names of the
    segments. If \code{NULL}, the names of \code{dna_segs} will be
    used, if available. Else, no name are plotted. If a \code{tree} is
    given, names must exist either in \code{dna_seg_labels} or in the
    names of \code{dna_segs}.
  }
  \item{dna_seg_label_cex}{
    A numeric. The character size for the DNA segments labels, or tree
    labels. Default is 1.
  }
  \item{dna_seg_label_col}{
    A color, of length 1 or of the same length as \code{dna_segs}. Gives
    the color of the labels. Default is black.
  }
  \item{gene_type}{
    A character. Describes the type of representation of genes or
    \code{dna_seg} elements. See details.
  }
  \item{arrow_head_len}{
    A numeric. Gives the length of arrow heads for gene type
    "arrows". The arrow head extends at maximum at half of the gene. Set
    to \code{Inf} to have all arrow heads covering the half of the
    gene. 200 by default.
  }
  \item{dna_seg_line}{
    A vector, either logical or giving colors, of length 1 or of same
    length as \code{dna_segs}. Should the line in the middle of the
    segments be drawn, and if yes, in what color. \code{TRUE} by
    default, which gives black lines. \code{FALSE} (logical, or as a
    string) results in no plotting.
  }
  \item{scale}{
    A logical. Should the scale be displayed on the plot. \code{TRUE} by
    default. 
  }
  \item{dna_seg_scale}{
    A logical, of length one or of the same length as
    \code{dna_segs}. Should a scale be displayed below each or all dna
    segments, respectively. \code{!scale} by default. 
  }
  \item{n_scale_ticks}{
    A integer. The (approximate) number of ticks on the longest
    segment. Default: 7. 
  }
  \item{scale_cex}{
    A numeric. The character size for the scale labels. Default is 1.
  }
  \item{global_color_scheme}{
    A character of length 4. If no \code{col} column is present on any
    comparison or is \code{override_color_schemes} is set, apply a
    global color scheme over all comparions. See below for more details.
    \code{c("auto", "auto", "blue_red")} by default. 
  }
  \item{override_color_schemes}{
    A logical. If \code{TRUE}, apply a global color scheme even if there
    are comparisons that have \code{col} columns. \code{FALSE} by
    default. 
  }
  \item{plot_new}{
    Logical. Produce a new plot? If \code{TRUE}, uses
    \code{grid.newpage} before plotting.
  }
  \item{debug}{
    A numeric. If > 0, only that number of element will be plotted for
    each \code{dna_seg} and \code{comparison}.
  }
  \item{\dots}{
    Further arguments to be passed to user-defined graphical functions.
  }
}
\details{
  One line is plotted per \code{\link{dna_seg}}. Eventually, the space 
  between the lines will be filled with the
  \code{\link{comparison}}s. \code{dna_seg}s can be annotated with
  \code{\link{annotation}}s, and accompagnying data can be plotted using
  \code{\link{seg_plot}}. 
  
  A phylogenetic tree (a \code{phylog} object from package \code{ade4})
  can be drawn at the left of the plot. The tree does not need to be
  ordered as the \code{dna_seg_labels}, but a permutation of the tree
  with that order should exist. If the tree is large, the number of
  permutations become too large, and the function will stop (>100000
  permutations). The solution is then to provide segments that are
  ordered in the same manner as the tree labels (or vice-versa).

  There is an (experimental) support for branch annotations. These are
  given in the Newick tree, directly after the parenthesis closing a
  node. They can be characters or integers, but so far
  \code{newick2phylog} doesn't support decimal values. Tags will be
  ignored if they start with "I", and trimmed if they start with "X".

  The format of the elements of \code{dna_segs} is previously determined
  in the object or can be globally set by \code{gene_type}. See the
  function \code{\link{gene_types}} to return the available types. Gene
  type can also be user-defined, using a function returning a
  \code{grob}. See \code{\link{gene_types}} for more details.

  \code{xlims} allow the user to plot subsegments of a
  \code{dna_seg}. \code{xlims} consists of a list composed of as many
  numeric vectors as there are segments. Each of these numeric vectors
  give pairs of left and right borders, and gives the
  direction. For example, c(1,2,6,4) will plot two subsegments, segment
  1 to 2 which is plotted left to right and segment 4 to 6, plotted
  right to left. \code{-Inf} and \code{Inf} values are
  accepted. \code{NULL} values will result in plotting the whole
  segment.

  \code{offsets} allows to user to define the placement of the
  subsegments. If a list is provided, each element of the list should
  have as many elements as there are subsegments. It will give the size
  of the gaps, including the first one from the border of the plot to
  the first subsegment.

  A main title (\code{main}) can also be added at the top of the plot,
  at the position defined by \code{main_pos}. A general scale
  can be added at the bottom right of the plot (\code{scale}).
  
  \code{dna_seg_scale} gives the ability to plot scales on one, some or
  every segment. \code{c(TRUE, FALSE, TRUE)} will add scales to the
  first and third segments.

  The four elements of \code{global_color_scheme} are (i) which column
  serves as scale to apply the color  scheme, or "auto" (default);
  (ii) if the scale is "increasing" or "decreasing" (see
  \code{apply_color_scheme} for more details), or "auto" (default);
  (iii) the color scheme to apply; (iv) the transparency to apply (0.5
  by default).
}
\value{
  %% Change this when/if would return gTree
  Nothing. A lattice graphic is plotted on the current device.
}
\author{
  Lionel Guy \email{lionel.guy@ebc.uu.se}, Jens Roat Kultima
}
\note{
  This plotting function has been tested as far as possible, but given
  its complexity and that the package is young, bugs or strange behaviors
  are possible. Please report them to the author.

  As of 10/3/2010, support for viewing exons/introns is only available
  using genbank and embl formats, not when importing ptt files.
}
\seealso{
  \code{\link{dna_seg}} and \code{\link{comparison}} for the base objects;
  \code{\link{read_dna_seg_from_tab}}, \code{\link{read_dna_seg_from_ptt}},
  \code{\link{read_comparison_from_tab}} and
  \code{\link{read_comparison_from_blast}} to read from files;
  \code{\link{annotation}} to annotate \code{dna_seg}s;
  \code{\link{seg_plot}} to draw plots next to \code{dna_seg}s;
  \code{\link{gene_types}} for \code{gene_type} argument;
  \code{\link{apply_color_scheme}} for color schemes;
}
\examples{
old.par <- par(no.readonly=TRUE)
data("three_genes")

## Segments only
plot_gene_map(dna_segs=dna_segs) 

## With comparisons
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons) 

## Tree
names <- c("A_aaa", "B_bbb", "C_ccc")
names(dna_segs) <- names
tree <- newick2phylog("(((A_aaa:4.2,B_bbb:3.9):3.1,C_ccc:7.3):1);")
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              tree=tree)
## Increasing tree width
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              tree=tree, tree_width=3)
## Annotations on the tree
tree2 <- newick2phylog("(((A_aaa:4.2,B_bbb:3.9)97:3.1,C_ccc:7.3)78:1);")
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              tree=tree2, tree_width=3)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              tree=tree2, tree_width=3, tree_branch_labels_cex=0.5)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              tree=tree2, tree_width=3, tree_branch_labels_cex=0)

## Annotation
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

## xlims
## Just returning a segment
plot_gene_map(dna_segs, comparisons,
              xlims=list(NULL, NULL, c(Inf,-Inf)),
              dna_seg_scale=TRUE)
## Removing one gene
plot_gene_map(dna_segs, comparisons,
              xlims=list(NULL, NULL, c(-Inf,2800)),
              dna_seg_scale=TRUE)

## offsets
offsets <- c(0, 0, 0)  
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, offsets=offsets)
offsets <- c(200, 400, 0)  
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, offsets=offsets)

## main
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              main="Comparison of A, B and C")
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              main="Comparison of A, B and C", main_pos="left")

## dna_seg_labels
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              dna_seg_labels=c("Huey", "Dewey", "Louie"))

## dna_seg_labels size
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              dna_seg_labels=c("Huey", "Dewey", "Louie"),
              dna_seg_label_cex=2)

## dna_seg_line
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              dna_seg_line=c("FALSE", "red", grey(0.6)))

## gene_type
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              gene_type="side_blocks")

##
## From here on, using a bigger dataset from a 4-genome comparison
##
data("barto")
## Adding a tree
tree <- newick2phylog("(BB:2.5,(BG:1.8,(BH:1,BQ:0.8):1.9):3);")
## Showing only subsegments
xlims1 <- list(c(1380000, 1445000),
               c(10000, 83000),
               c(15000, 98000),
               c(5000, 82000))
## Reducing dataset size for speed purpose
for (i in 1:length(barto$dna_segs)){
    barto$dna_segs[[i]] <- trim(barto$dna_segs[[i]], xlim=xlims1[[i]])
    if (i < length(barto$dna_segs))
        barto$comparisons[[i]] <- trim(barto$comparisons[[i]],
                                      xlim1=xlims1[[i]], xlims1[[i+1]])   
}
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree,
              xlims=xlims1,
              dna_seg_scale=TRUE)
## Showing several subsegments per genome
xlims2 <- list(c(1445000, 1415000, 1380000, 1412000),
               c(  10000,   45000,   50000,   83000, 90000, 120000),
               c(  15000,   36000,   90000,  120000, 74000,  98000),
               c(   5000,    82000))

## dna_seg_scale, global_color_scheme, size, number, color of dna_seg_scale,
## size of dna_seg_scale labels
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree, xlims=xlims2,
              dna_seg_scale=c(TRUE, FALSE, FALSE, TRUE), scale=FALSE,
              dna_seg_label_cex=1.7,
              dna_seg_label_col=c("black", "grey", "blue", "red"),
              global_color_scheme=c("e_value", "auto", "grey", "0.7"),
              n_scale_ticks=3, scale_cex=1)

## Hand-made offsets: size of all gaps
offsets2 <- list(c(10000, 10000),
                 c(2000, 2000, 2000),
                 c(10000, 5000, 2000),
                 c(10000))
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree,
              #annotations=annots,
              xlims=xlims2,
              offsets=offsets2,
              dna_seg_scale=TRUE)
##
## Exploring and modifying a previously plotted gene map plot
##
## View viewports
current.vpTree()
## Go down to one of the viewports, add an xaxis, go back up to root viewport
downViewport("dna_seg_scale.3.2")
grid.rect()
upViewport(0)
## Get all the names of the objects
grobNames <- getNames()
grobNames
## Change the color ot the scale line
grid.edit("scale.lines", gp=gpar(col="grey"))
## Remove first dna_seg_lines
grid.remove("dna_seg_line.1.1")

##
## Plot genoPlotR logo
##
col <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
         "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
cex <- 2.3
## First segment 
start1 <- c(150, 390, 570)
end1   <- c(  1, 490, 690)
genoR <- c(270, 530)
## Second segment
start2 <- c(100, 520, 550)
end2   <- c(240, 420, 650)
Plot <- c(330)
## dna_segs
ds1 <- as.dna_seg(data.frame(name=c("", "", ""),
                             start=start1, end=end1, strand=rep(1, 3),
                             col=col[c(2, 6, 1)], stringsAsFactor=FALSE))
ds_genoR <- as.dna_seg(data.frame(name=c("geno", "R"),
                             start=genoR, end=genoR, strand=rep(1, 2),
                             col=c(col[8], "black"),
                             stringsAsFactor=FALSE), cex=cex, gene_type="text")
ds2 <- as.dna_seg(data.frame(name=c("", "", ""),
                             start=start2, end=end2, strand=rep(1, 3),
                             col=col[c(5, 3, 7)],
                             stringsAsFactor=FALSE))
ds_Plot <- as.dna_seg(data.frame(name="Plot",
                             start=Plot, end=Plot, strand=1,
                             col=col[c(1)],
                             stringsAsFactor=FALSE), cex=cex, gene_type="text")
## comparison
c1 <- as.comparison(data.frame(start1=start1, end1=end1,
                               start2=start2, end2=end2,
                               col=grey(c(0.6, 0.8, 0.5))))
## Generate genoPlotR logo
\dontrun{
pdf("logo.pdf", h=0.7, w=3)
}
par(fin=c(0.7, 3))
plot_gene_map(dna_segs=list(c(ds1, ds_genoR), c(ds2, ds_Plot)),
              comparisons=list(c1), scale=FALSE, dna_seg_scale=FALSE,
              dna_seg_line=grey(0.7), offsets=c(-20,160))
\dontrun{
dev.off()
}
par(old.par)
}
\keyword{ hplot }
