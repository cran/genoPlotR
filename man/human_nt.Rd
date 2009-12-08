\name{human_nt}
\alias{human_nt}
\title{
  Human-readable nucleotide scale
}
\description{
  Return a human readable list from a nucleotide position or lenght. 
}
\usage{
human_nt(nt, signif = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nt}{
    A nucleotide position
  }
  \item{signif}{
    Either a logical or an integer. If \code{FALSE} (default), \code{nt}
    is not rounded. Else, it returns \code{signif} significant digits.
  }
}
\details{
  Return a nucleotide value in nt, kb, Mb or Gb, according to the value given.
  This is particularly useful to display nice scales without too many 
  trailing zeros.
}
\value{
  Returns a list with 4 elements
  \item{n}{A numeric value corresponding to \code{nt} divided by
    \code{mult} (see below).}
  \item{tag}{A character, giving the multiplier used in text.}
  \item{mult}{The muliplier used, in numeric value.}
  \item{text}{A character, giving the value in a human readable format.}

}
\author{
  Lionel Guy
}
\examples{
human_nt(123456)
human_nt(123456, signif=2)
human_nt(123456890, signif=2)
}
\keyword{ data }
