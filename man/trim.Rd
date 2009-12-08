\name{trim}
\alias{trim}
\alias{trim.default}
\alias{trim.dna_seg}
\alias{trim.comparison}
\alias{trim.annotation}
\alias{trim.seg_plot}
\title{
  Trimming data frames or more complex objects with >= 2 numeric columns 
}
\description{
  Trims data frames with 2 or more numeric columns using a
  xlim. \code{xlim}(s) are as used to filter rows whose numeric values are
  included in this interval.
}
\usage{
trim(x, ...)
\method{trim}{default}(x, xlim = NULL, ...)
\method{trim}{dna_seg}(x, xlim = NULL, ...)
\method{trim}{comparison}(x, xlim1 = c(-Inf, Inf), xlim2 = c(-Inf, Inf), ...)
\method{trim}{annotation}(x, xlim = NULL, ...)
\method{trim}{seg_plot}(x, xlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object to trim,. generally a data frame or a matrix, or a
    \code{seg_plot} object.
  }
  \item{xlim}{
    A numeric of length 2. In a general case, the rows whose values are
    included in this interval are returned.
  }
  \item{\dots}{
    Unused.
  }
  \item{xlim1}{
    A numeric of length 2. In the case of comparison, where the
    comparison can be filtered on two sides, the interval to filter the
    first side. 
  }
  \item{xlim2}{
    A numeric of length 2. The interval to filter the second side.
  }
}
\details{
  In the case where \code{x} is a \code{seg_plot} object, the function
  uses the \code{xargs} argument to define what are the vectors defining
  the x position (they should be the same length). Then, all the
  arguments (including those inside an eventual \code{gp} argument) that
  are the same length as the x vectors are trimmed, so that only the
  rows for which the x values are inside the \code{xlim} argument are kept.
}
\value{
  Returns the same object as input, with the rows (or subset)
  corresponding to the given interval.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{dna_seg}}, \code{\link{comparison}}, \code{\link{seg_plot}}.

}
\examples{
## Load
data(barto)
xlim_ref <- c(10000, 45000)
## Seg 2 (ref)
barto$dna_segs[[2]] <- trim(barto$dna_segs[[2]], xlim=xlim_ref)
## Seg 1
barto$comparisons[[1]] <- trim(barto$comparisons[[1]], xlim2=xlim_ref)
xlim1 <- range(barto$comparisons[[1]], overall=FALSE)$xlim1
barto$dna_segs[[1]] <- trim(barto$dna_segs[[1]], xlim=xlim1)
## Seg 3
barto$comparisons[[2]] <- trim(barto$comparisons[[2]], xlim1=xlim_ref)
xlim3 <- range(barto$comparisons[[2]], overall=FALSE)$xlim2
barto$dna_segs[[3]] <- trim(barto$dna_segs[[3]], xlim=xlim3)
## Seg 4
barto$comparisons[[3]] <- trim(barto$comparisons[[3]], xlim1=xlim3)
xlim4 <- range(barto$comparisons[[3]], overall=FALSE)$xlim2
barto$dna_segs[[4]] <- trim(barto$dna_segs[[4]], xlim=xlim4)
## Plot
plot_gene_map(barto$dna_segs, barto$comparisons)

## With seg_plot
x <- 1:20
y <- rnorm(20)
sp <- seg_plot(func=pointsGrob, args=list(x=x, y=y,
                                  gp=gpar(col=1:20, cex=1:3)))
## Trim 
sp_trim <- trim(sp, c(3, 10))
str(sp_trim)
range(sp_trim$arg$x)
}
\keyword{ data }
