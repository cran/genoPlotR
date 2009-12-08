\name{seg_plot}
\alias{seg_plot}
\alias{as.seg_plot}
\alias{is.seg_plot}
\title{
  seg_plot class and class functions
}
\description{
  An seg_plot is an object to plot data associated to a \code{dna_seg}
  object. It is a \code{list} with mandatory and optional arguments. The
  main arguments are \code{func}, which is a function returning a
  \code{\link{grob}} or a \code{\link{gList}}, and \code{args}, which
  are arguments to be passed to this function.
}
\usage{
seg_plot(func,
         args = NULL,
         xargs = c("x", "x0", "x1", "x2", "v"),
         yargs = c("y", "y0", "y1", "y2", "h"),
         ylim = NULL)
as.seg_plot(seg_plot)
is.seg_plot(seg_plot)
}
\arguments{
  \item{func}{
    Mandatory, with no defaults. A function that returns a \code{grob}
    object. See \code{\link{grid}} documentation to find ready-made
    functions. User-defined functions work too.
  }
  \item{args}{
    A list, \code{NULL} by default. The arguments that will be passed to
    the function. It is recommended that all arguments are named.
  }
  \item{xargs}{
    A vector giving the names of which of the arguments in \code{args}
    are defining the x-axis. Used, among others, by the function
    \code{\link{trim.seg_plot}}. By default, gives the most common
    x-defining arguments of the \code{grid} functions (x, x0, x1, x2, v). 
  }
  \item{yargs}{
    A vector giving the names of which of the arguments in \code{args}
    are defining the y-axis. Used when plotting the graphs to define
    a sensible \code{ylim} if not defined. By default, gives the most common
    y-defining arguments of the \code{grid} functions (y, y0, y1, y2, h). 
  }
  \item{ylim}{
    A numeric vector of length 2, defining the range of the plot when
    drawn with \code{\link{plot_gene_map}}. Derived from \code{yargs} if
    not set.
  }
  \item{seg_plot}{
    In \code{as.seg_plot}, a \code{list} object to convert to
    \code{seg_plot}. See details below.

    In \code{is.seg_plot}, an object to test.
  }
}
\value{
  \code{seg_plot} and \code{as.seg_plot} return a \code{seg_plot} object.
  \code{is.seg_plot} returns a logical.
}
\details{
  A \code{seg_plot} object is an object describing how to plot data
  associated to a \code{dna_seg}. It is a list composed of a function,
  arguments to pass to this function, two arguments to define which of
  those define x and y, and an eventual \code{ylim} to limit the
  plotting to a certain range when plotting.

  The function \code{func} should return a \code{grob} object, or a
  \code{gList} list of \code{grob}s. The predefined functions of
  \code{grid}, such as \code{linesGrob}, \code{pointsGrob},
  \code{segmentsGrob}, \code{textGrob} or \code{polygonGrob} can be
  used, or user-defined functions can be defined.

  The arguments in \code{args} should correspond to arguments passed to
  \code{func}. For example, if \code{func = pointsGrob}, \code{args}
  could contain the elements \code{x = 10:1}, \code{y = 1:10}. It will
  often also contain a \code{gp} element, the result of a call to the
  \code{\link{gpar}} function, to control graphical aspects of the plot
  such as color, fill, line width and style, fonts, etc.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{plot_gene_map}}, \code{\link{trim.seg_plot}}.
}
\examples{
## Using the existing pointsGrob
x <- 1:20
y <- rnorm(20)
sp <- seg_plot(func=pointsGrob, args=list(x=x, y=y,
                                  gp=gpar(col=1:20, cex=1:3)))
is.seg_plot(sp)
## Function seg_plot(...) is identical to as.seg_plot(list(...))
sp2 <- as.seg_plot(list(func=pointsGrob, args=list(x=x, y=y,
                                           gp=gpar(col=1:20, cex=1:3))))
identical(sp, sp2)
## For the show, plot the obtained result
grb <- do.call(sp$func, sp$args)
## Trim the seg_plot
sp_trim <- trim(sp, c(3, 10))
## Changing color and function "on the fly"
sp_trim$args$gp$col <- "blue"
sp_trim$func <- linesGrob
grb_trim <- do.call(sp_trim$func, sp_trim$args)
## Now plot
plot.new()
pushViewport(viewport(xscale=c(0,21), yscale=c(-4,4)))
grid.draw(grb)
grid.draw(grb_trim)

## Using home-made function
triangleGrob <- function(start, end, strand, col, ...) {
  x <- c(start, (start+end)/2, end)
  y1 <- 0.5 + 0.4*strand
  y <- c(y1, rep(0.5, length(y1)), y1)
  polygonGrob(x, y, gp=gpar(col=col), default.units="native",
              id=rep(1:7, 3))
}
start <- seq(1, 19, by=3)+rnorm(7)/3
end <- start + 1 + rnorm(7)
strand <- sign(rnorm(7))
sp_tr <- seg_plot(func=triangleGrob,
                  args=list(start=start, end=end, strand=strand,
                    col=1:length(start)), xargs=c("start", "end"))
grb_tr <- do.call(sp_tr$func, sp_tr$args)
plot.new()
pushViewport(viewport(xscale=c(1,22), yscale=c(-2,2)))
grid.draw(grb_tr)
## Trim
sp_tr_trim <- trim(sp_tr, xlim=c(5, 15))
str(sp_tr_trim)
## If the correct xargs are not indicated, trimming won't work
sp_tr$xargs <- c("x")
sp_tr_trim2 <- trim(sp_tr, xlim=c(5, 15))
identical(sp_tr_trim, sp_tr_trim2)

y1 <- convertY(grobY(grb_tr, "south"), "native")
y2 <- convertY(grobY(grb_tr, "north"), "native")
heightDetails(grb)
grb

## Applying it to plot_gene_maps
data(three_genes)

## Build data to plot
xs <- lapply(dna_segs, range)
colors <- c("red", "blue", "green")

seg_plots <- list()
for (i in 1:length(xs)){
  x <- seq(xs[[i]][1], xs[[i]][2], length=20)
  seg_plots[[i]] <- seg_plot(func=pointsGrob,
                             args=list(x=x, y=rnorm(20)+2*i,
                               default.units="native", pch=3,
                               gp=gpar(col=colors[i], cex=0.5)))
}
plot_gene_map(dna_segs, comparisons,
              seg_plots=seg_plots,
              seg_plot_height=0.5,
              seg_plot_height_unit="inches",
              dna_seg_scale=TRUE)

## A more complicated example
data(barto)
tree <- newick2phylog("(BB:2.5,(BG:1.8,(BH:1,BQ:0.8):1.9):3);")
## Showing several subsegments per genome
xlims2 <- list(c(1445000, 1415000, 1380000, 1412000),
               c(  10000,   45000,   50000,   83000, 90000, 120000),
               c(  15000,   36000,   90000,  120000, 74000,  98000),
               c(   5000,    82000))

## Adding fake data in 1kb windows
seg_plots <- lapply(barto$dna_segs, function(ds){
  x <- seq(1, range(ds)[2], by=1000)
  y <- jitter(seq(100, 300, length=length(x)), amount=50)
  seg_plot(func=linesGrob, args=list(x=x, y=y, gp=gpar(col=grey(0.3), lty=2)))
})
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree,
              seg_plots=seg_plots,
              seg_plot_height=0.5,
              seg_plot_height_unit="inches",
              xlims=xlims2,
              limit_to_longest_dna_seg=FALSE,
              dna_seg_scale=TRUE,
              main="Random plots for the same segment in 4 Bartonella genomes")

}
\keyword{ data }
