\name{artemisColors}
\alias{artemisColors}
\title{
  Artemis Colors
}
\description{
  Returns a data frame with the standard artemis colors.
}
\usage{
artemisColors()
}
\value{
  A \code{data.frame} with the following columns: \code{n}, \code{names},
  \code{colors}, \code{r}, \code{g} and \code{g}. The 3 first columns
  give the Artemis color number, its name, and its equivalent in R. The
  3 last give the r, g and b values. 
}
\author{
  Lionel Guy
}
\references{
  Artemis website: http://www.sanger.ac.uk/resources/software/artemis/
}
\examples{
artCol <- artemisColors()
plot(rep(1, nrow(artCol)), artCol$n, xlim=c(1, 2), type="n")
text(rep(1, nrow(artCol)), artCol$n, labels=artCol$n, col=artCol$colors)
text(rep(1, nrow(artCol)), artCol$n, labels=artCol$names, col=artCol$colors,
     pos=4, offset=1)

}
\keyword{ data }
