\name{dna_seg}
\alias{dna_seg}
\alias{as.dna_seg}
\alias{is.dna_seg}
\title{
  DNA segment (dna_seg) class and class functions
}
\description{
  A DNA segment is a collection of genes or elements along a genome, to
  be represented on a map. These functions are class functions to
  create, convert and test dna_seg objects.
}
\usage{
dna_seg(x, ...)
as.dna_seg(df, col = "blue", fill = "transparent", lty = 1, lwd = 1,
           pch = 8, cex = 1, gene_type = "arrows")
is.dna_seg(dna_seg)
}
\arguments{
  \item{x}{
    A \code{data.frame} or \code{list} that can be coerced to a data frame.
  }
  \item{\dots}{
    Arguments further passed to \code{as.dna_seg} (see below).
  }
  \item{df}{
    A data frame representing the \code{dna_seg} object. See details for
    necessary columns.
  }
  \item{col}{
    Either a color vector of the same length as \code{df} or of length
    one, to be applied to the whole object. Will be passed to
    corresponding \code{grob} object. Default to blue.
  }
  \item{fill}{
    Either a fill vector of the same length as \code{df} or of length
    one, to be applied to the whole object. Will be passed to
    corresponding \code{grob} object. Default to \code{NULL}. If
    \code{NULL}, then \code{col} is used to set the \code{fill}
    parameter.
  }
  \item{lty}{
    A vector of the same length as \code{df} or of length one, giving
    the line type around the objects.
  }
  \item{lwd}{
    Same as \code{lty}, giving the line width.
  }
  \item{pch}{
    Same as \code{lty}, giving the character representing each
    object. Goes with \code{gene_type} points.
  }
  \item{cex}{
    Same as \code{lty}, giving the character size representing each
    object. Goes with \code{gene_type} points.
  }
  \item{gene_type}{
    Vector of the same length as \code{df} or of length one, giving the
    type of representation of each object.
  }
  \item{dna_seg}{
    Object to test.
  }
}
\details{
  Objects to be converted needs to have their first 4 columns named 
  \code{name}, \code{start}, \code{end} and \code{strand}. Extra columns
  with names \code{col}, \code{lty}, \code{lwd}, \code{pch}, \code{cex},
  \code{gene_type} will be used in the plotting process. Other extra
  columns will be kept in the object, but not used.

  \code{dna_seg} tries to build a dna_seg object from a data frame or a
  list. 

  \code{as.dna_seg} tries to build a dna_seg object from a data frame.
  
  Read functions such as \code{\link{read_dna_seg_from_tab}} and
  \code{\link{read_dna_seg_from_ptt}} also return dna_seg objects.  
}
\value{
  A comparison object for \code{comparison} and
  \code{as.comparison}. DNA seg objects are also of class
  \code{data.frame}. They contain the following columns: \code{name},
  \code{start}, \code{end}, \code{strand}, \code{col}, \code{lty}, \code{lwd},
  \code{pch}, \code{cex}, \code{gene_type}.

  A logical for \code{is.comparison}.
}
\author{
  Lionel Guy
}
\seealso{
  \code{\link{read_dna_seg_from_tab}},
  \code{\link{read_dna_seg_from_ptt}}, \code{\link{gene_types}}.
}
\examples{
## generate data
names1 <- c("feat1", "feat2", "feat3")
starts1 <- c(2, 1000, 1050)
ends1 <- c(600, 800, 1345)
strands1 <- c("-", -1, 1)
cols1 <- c("blue", "grey", "red")

## create data.frame
df1 <- data.frame(name=names1, start=starts1, end=ends1,
                  strand=strands1, col=cols1)

## with dna_seg
dna_seg1 <- dna_seg(df1)
dna_seg1
as.dna_seg(df1)

## test
is.dna_seg(dna_seg1)

## directly readable with read_dna_seg_from_tab
\dontrun{
write.table(x=dna_seg1, file="dna_seg1.tab", quote=FALSE,
            row.names=FALSE, sep="\t")
}

## with only one gene and with list, or two, and merging with c.dna_seg
gene2a <- dna_seg(list(name="feat1", start=50, end=900, strand="-", col="blue"))
genes2b <- dna_seg(data.frame(name=c("feat2", "feat3"), start=c(800, 1200),
                              end=c(1100, 1322), strand=c("+", 1),
                              col=c("grey", "red"),
                              fill=c("transparent", "grey"),
                              gene_type=c("arrows", "blocks")))
dna_seg2 <- c(gene2a, genes2b)

## test
is.dna_seg(dna_seg2)

## reading from file
dna_seg3_file <- system.file('extdata/dna_seg3.tab', package = 'genoPlotR')
dna_seg3 <- read_dna_seg_from_tab(dna_seg3_file)
is.dna_seg(dna_seg3)

}
\keyword{ data }
