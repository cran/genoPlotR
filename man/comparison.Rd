\name{comparison}
\alias{comparison}
\alias{as.comparison}
\alias{is.comparison}
\title{
  Comparison class and class functions
}
\description{
  A comparison is a collection of similarities, representing the
  comparison between two DNA segments. These functions are class
  functions to create, convert and test comparison objects.
}
\usage{
comparison(x)
as.comparison(df)
is.comparison(comparison)
}
\arguments{
  \item{x}{
    Can be a \code{list} or \code{data.frame} object. See the
    details for the columns in the \code{data.frame}.
  }
  \item{df}{
    A \code{data.frame} object. See details for the required columns.
  }
  \item{comparison}{
    Any object to test.
  }
}
\details{
  Objects (either data frames or lists) should have at least named
  elements \code{start1}, \code{end1}, \code{start2} and \code{end2}. In
  addition, it can take a \code{col}or column. Additional numeric columns can
  be used for color-coding (via \code{apply_color_scheme}.

  \code{comparison} tries to build a comparison object from either a
  data frame or a list, \code{as.comparison} accepts only data.frames.

  \code{is.comparison} returns \code{TRUE} if the object tested is a
  comparison object.
  
  Read functions such as \code{\link{read_comparison_from_tab}} and
  \code{\link{read_comparison_from_blast}} also return comparison objects.
}
\value{
  A comparison object for \code{comparison} and
  \code{as.comparison}. Comparison objects are also of class
  \code{data.frame}. They contain the columns  \code{start1},
  \code{end1}, \code{start2}, \code{end2}, \code{direction} and
  \code{col} (color).

  A logical for \code{is.comparison}.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{dna_seg}}, \code{\link{read_comparison_from_tab}},
  \code{\link{read_comparison_from_blast}},
  \code{\link{trim.comparison}}, \code{\link{reverse.comparison}}.
}
\examples{
## Get some values
starts1 <- c(2, 1000, 1050)
ends1 <- c(600, 800, 1345)
starts2 <- c(50, 800, 1200)
ends2 <- c(900, 1100, 1322)

## From a data.frame
comparison1 <- as.comparison(data.frame(start1=starts1, end1=ends1,
                                        start2=starts2, end2=ends2))
comparison1
is.comparison(comparison1)
is.data.frame(comparison1)
comparison(data.frame(start1=starts1, end1=ends1,
                      start2=starts2, end2=ends2))
## From a list
comparison(list(start1=starts1, end1=ends1,
                start2=starts2, end2=ends2))

## From a file
comparison2_file <- system.file('extdata/comparison2.tab',
                                package = 'genoPlotR')
comparison2 <- read_comparison_from_tab(comparison2_file)

}
\keyword{ data }
