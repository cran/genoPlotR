\name{apply_color_scheme}
\alias{apply_color_scheme}
\title{
  Apply a color scheme
}
\description{
  Apply a color scheme to a numeric vector, eventually taking the
  direction into account.
}
\usage{
apply_color_scheme(x, direction = NULL, color_scheme = "grey",
decreasing = FALSE, rng = NULL, transparency = 0.5)
}
\arguments{
  \item{x}{
    A numeric, that will be used to apply a gradient of colors to a
    comparison. 
  }
  \item{direction}{
    If a red-blue scheme is choosen, the vector (composed of -1 and 1
    values and of same length as \code{x}) giving the direction of the
    comparison. 
  }
  \item{color_scheme}{
    Character. One of \code{red_blue}, \code{blue_red}, \code{grey},
    \code{gray}. 
  }
  \item{decreasing}{
    Logical. Are the values of the comparisons oriented such as the
    lower the value, the closer the relationship (e.g. e-values, gaps,
    mismatches, etc)? \code{FALSE} by default. 
  }
  \item{rng}{
    Numeric of length 2. Gives the higher and lower limit to apply a
    color scheme.  
  }
  \item{transparency}{
    Numeric of length 1, between 0 and 1, or FALSE. Should the color
    scheme use transparency, and if yes how much (ratio). 0.5 by
    default.
    Not supported on all devices.
  }
}
\details{
  A color scale is calculated, with the darker color corresponding to
  the highest values of \code{x}, or the contrary is \code{decreasing}
  is \code{TRUE}. For the moment, two schemes (red-blue and grey scale)
  are used.

  For the red-blue scale (as in ACT), the direct comparisons are
  colored in red hues, and the reversed ones in blue hues.

  This is especially useful to replace comparison values (such as BLAST
  percent identity values) by color hues.
}
\value{
  A character vector of same length as \code{x}, representing colors.
}
\references{
  Artemis Comparison Tool, \url{http://www.sanger.ac.uk/Software/ACT/}
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{comparison}}
}
\examples{
## Load data
data(three_genes)

## Color schemes
## Greys
comparisons[[1]]$values <- c(70, 80, 90)
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           color_scheme="grey")
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
## Red-blue
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           direction=comparisons[[1]]$direction,
                                           color_scheme="red_blue")
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
## Decreasing
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           direction=comparisons[[1]]$direction,
                                           color_scheme="red_blue",
                                           decreasing=TRUE)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
## Range
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           direction=comparisons[[1]]$direction,
                                           color_scheme="red_blue",
                                           rng=c(30,100))
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
## Transparency
x1 <- seq(100, 600, by=50)
x2 <- seq(1100, 700, by=-50)
comparisons[[2]] <- as.comparison(data.frame(start1=c(x1, x2),
                                             end1=c(x1+250, x2+300),
                                             start2=c(x1+150, x2-300)+2000,
                                             end2=c(x1+250, x2-500)+2000
                                             ))
comparisons[[2]]$col <- apply_color_scheme(1:nrow(comparisons[[2]]),
                                           comparisons[[2]]$direction,
                                           color_scheme="blue_red")
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           color_scheme="grey",
                                           transparency=0.8)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values,
                                           color_scheme="grey",
                                           transparency=1)
comparisons[[2]]$col <- apply_color_scheme(1:nrow(comparisons[[2]]),
                                           comparisons[[2]]$direction,
                                           color_scheme="blue_red",
                                           transparency=0.2)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
}
\keyword{ data }
