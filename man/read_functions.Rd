\name{read_functions}
\alias{read_functions}
\alias{read_dna_seg_from_tab}
\alias{read_dna_seg_from_file}
\alias{read_dna_seg_from_embl}
\alias{read_dna_seg_from_genbank}
\alias{read_dna_seg_from_ptt}
\alias{read_comparison_from_tab}
\alias{read_comparison_from_blast}
\alias{read_mauve_backbone}
\title{Reading functions}
\description{
  Functions to parse dna_seg objects from tab, embl, genbank, ptt files or
  from mauve backbone files, and comparison objects from tab or blast
  files.
}
\usage{
read_dna_seg_from_tab(file, header = TRUE, ...)
read_dna_seg_from_file(file, tagsToParse=c("CDS"), fileType = "detect",
                       meta_lines = 2, gene_type = "auto", header = TRUE, ...)
read_dna_seg_from_embl(file, tagsToParse=c("CDS"), ...)
read_dna_seg_from_genbank(file, tagsToParse=c("CDS"), ...)
read_dna_seg_from_ptt(file, meta_lines = 2, header = TRUE, ...)
read_comparison_from_tab(file, header = TRUE, ...)
read_comparison_from_blast(file, sort_by = "per_id",
                           filt_high_evalue = NULL,
                           filt_low_per_id = NULL,
                           filt_length = NULL,
                           color_scheme = NULL, ...)
read_mauve_backbone(file, ref = 1, gene_type = "side_blocks",
                    header = TRUE, filter_low = 0, ...)
}
\arguments{
  \item{file}{
    Path to file to load. URL are accepted.
  }
  \item{header}{
    Logical. Does the tab file has headers (column names)?
  }
  \item{tagsToParse}{
    Character vector. Tags to parse in embl or genbank files.
    Common tags are 'CDS', 'gene', 'misc_feature'.
  }
  \item{fileType}{
    Character string. Select file type, could be 'detect' for automatic
    detection, 'embl' for embl files, 'genbank' for genbank files or
    'ptt' for ptt files.
  }
  \item{meta_lines}{
    The number of lines in the ptt file that represent "meta" data, not
    counting the header lines. Standard for NCBI files is 2 (name and
    length, number of proteins. Default is also 2.
  }
  \item{gene_type}{
    Determines how genes are visualized. If 'auto' genes will appear as
    arrows in there are no introns and as blocks if there are introns.
    Can also be set to for example 'blocks' or 'arrows'. Do note, currently
    introns are not supported in the ptt file format.
    Default for mauve backbone is \code{side_blocks}.
    See \code{\link{gene_types}} page for more details, or use function
    \code{gene_types}. 
  }
  \item{sort_by}{
    In BLAST-like tabs, gives the name of the column that will be used
    to sort the comparisons. Accepted values are \code{per_id}
    (percent identity, default), \code{mism} (mismatches), \code{gaps}
    (gaps), \code{e_value} (E-value), \code{bit_score} (bit score). 
  }
  \item{filt_high_evalue}{
    A numerical, or \code{NULL} (default). Filters out all comparisons
    that have a e-value higher than this one.
  }
  \item{filt_low_per_id}{
    A numerical, or \code{NULL} (default). Filters out all comparisons
    that have a percent identity lower than this one.
  }
  \item{filt_length}{
    A numerical, or \code{NULL} (default). Filters out all comparisons
    that have alignments shorter than this value.
  }
  \item{color_scheme}{
    A color scheme to apply. See \code{apply_color_scheme} for more
    details. Possible values include \code{grey} and \code{red_blue}.
    \code{NULL} by default. Color schemes can be applied while running
    \code{plot_gene_map}.
  }
  \item{ref}{
    In mauve backbone, which of the dna segments will be the reference,
    i.e. which one will have its blocks in order.
  }
  \item{\dots}{
    Further arguments passed to generic reading functions and class
    conversion functions. See \code{\link{as.dna_seg}} and
    \code{\link{as.comparison}}.

    For \code{read_comparison*} functions, see details.
  }
  \item{filter_low}{
    A numeric. If larger than 0, all blocks smaller that this number
    will be filtered out. Defaults to 0.
  }
}
\details{
  Tab files representing DNA segements should have at least the following
  columns: name, start, end, strand (in that order. Additionally, if the
  tab file has headers, more columns will be used to define, for
  example, the color, line width and type, pch and/or cex. See
  \code{\link{dna_seg}} for more information. An example:

  \code{
    name    start   end     strand  col  
    feat1A  2       1345    1       blue
    feat1B  1399    2034    1       red
    feat1C  2101    2932    -1      grey
    feat1D  2800    3120    1       green
  }

  Embl and Genbank files are two commonly used file types. These file types
  often contain a great variety of information. To properly extract data from
  these files, the user has to choose which features to extract. Commonly 'CDS'
  features are of interest, but other feature tags such as 'gene' or 'misc_feature'
  may be of interest. Should a feature contain an inner "pseudo" tag indicating this
  CDS or gene is a pseudo gene, this will be presented as a 'CDS_pseudo' or a 'gene_pseudo'
  feature type respectively in the resulting table. Certain constraints apply to these
  file types, of which some are: embl files must contain one and only one ID tag;
  genbank files may only contain one and only one locus tag.

  Ptt (or protein table) files are a tabular format giving a bunch of
  information on each protein of a genome (or plasmid, or virus,
  etc). They are available for each published genome on the NCBI ftp
  site (ftp://ftp.ncbi.nlm.nih.gov/genomes/). As an example, look at
  ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Bartonella_henselae_Houston-1/NC_005956.ptt. 
  
  Tabular comparison files should have at least the following
  columns: start1, end1, start2, end2. If no header is specified, the
  fifth column is parsed as the color.

  \tabular{lllll}{
    start1 \tab end1   \tab start2 \tab end2  \tab  col\cr
    2      \tab 1345   \tab 10     \tab 1210  \tab  red\cr
    1399   \tab 2034   \tab 2700   \tab 1100  \tab  blue\cr
    500    \tab 800    \tab 3000   \tab 2500  \tab  blue\cr
  }

  BLAST tabular result files are produced either with blastall using -m8
  or -m9 parameter, or with any of the newer
  blastn/blastp/blastx/tblastx using -outfmt 6 or -outfmt 7.

  For both comparison functions, several further arguments can be passed
  to the generic reading functions:
  \itemize{
    \item{color_scheme}{
      Passed to \code{\link{apply_color_scheme}}. See this function for
      supported values. Default to \code{NULL}.
    }
    \item{filt_low}{
      Filter low comparison values. Only similarities whose values in
      the \code{color_from} (for BLAST files) or \code{col} columns are
      greater than \code{filt_low} will be considered.
    }
    \item{filt_len}{
      Filter comparisons by length. Comparisons whose average length is
      below \code{filt_len} are not considered. 
    }
    \item{reverse}{
      If numeric and > 0, the comparison is passed to
      \code{\link{reverse.comparison}}.
    } 
  }
  
  Mauve backbone is another tabular data file that summarizes the blocks
  that are similar between all compared genomes. Each genome gets two
  columns, one start and one end of the block. There is one row per
  block and eventually a header row. If named, columns have sequence
  numbers, not actual names, so be careful to input the same order in
  both Mauve and genoPlotR. See
  http://asap.ahabs.wisc.edu/mauve-aligner/mauve-user-guide/mauve-output-file-formats.html
  for more info on the file format. Normally, the function should be
  able to read both \code{progressiveMauve} and \code{mauveAligner}
  outputs. The function returns both the blocks as \code{dna_seg}s and
  the links between the blocks as \code{comparison}s.

}
\value{
  \code{read_dna_seg_from_tab}, \code{read_dna_seg_from_file}, \code{read_dna_seg_from_embl},
  \code{read_dna_seg_from_genbank} and \code{read_dna_seg_from_ptt} return
  \code{dna_seg} objects. \code{read_comparison_from_tab} and
  \code{read_comparison_from_blast} return \code{comparison}
  objects. \code{read_mauve_backbone} returns a list containing a list
  of \code{dna_seg}s and \code{comparison}s.
  objects. 
}
\references{
  For BLAST: http://www.ncbi.nlm.nih.gov/blast/
  For Mauve: http://asap.ahabs.wisc.edu/mauve/
}
\author{
  Lionel Guy, Jens Roat Kultima
}
\note{
  Formats are changing and it maybe that some functions are temporarily
  malfunctioning. Please report any bug to the author. Mauve examples
  were prepared with Mauve 2.3.1.
}
\seealso{
  \code{\link{comparison}}, \code{\link{dna_seg}},
  \code{\link{apply_color_scheme}}.
}
\examples{
##
## From tabs
##
## Read DNA segment from tab
dna_seg3_file <- system.file('extdata/dna_seg3.tab', package = 'genoPlotR')
dna_seg3 <- read_dna_seg_from_tab(dna_seg3_file)

## Read comparison from tab
comparison2_file <- system.file('extdata/comparison2.tab',
                                package = 'genoPlotR')
comparison2 <- read_comparison_from_tab(comparison2_file)

##
## Mauve backbone
##
## File
bbone_file <- system.file('extdata/barto.backbone', package = 'genoPlotR')
## Read backbone
bbone <- read_mauve_backbone(bbone_file)
names <- c("B_bacilliformis", "B_grahamii", "B_henselae", "B_quintana")
names(bbone$dna_segs) <- names
## Plot
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons)

## Using filter_low & changing reference sequence
bbone <- read_mauve_backbone(bbone_file, ref=2, filter_low=2000) 
names(bbone$dna_segs) <- names
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons)

## Read guide tree
tree_file <- system.file('extdata/barto.guide_tree', package = 'genoPlotR')
tree_str <- readLines(tree_file)
for (i in 1:length(names)){
  tree_str <- gsub(paste("seq", i, sep=""), names[i], tree_str)
}
tree <- newick2phylog(tree_str)
## Plot
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons,
              tree=tree)

##
## From embl file
##
bq_embl_file <- system.file('extdata/BG_plasmid.embl', package = 'genoPlotR')
bq <- read_dna_seg_from_embl(bq_embl_file)

##
## From genbank file
##
bq_genbank_file <- system.file('extdata/BG_plasmid.gbk', package = 'genoPlotR')
bq <- read_dna_seg_from_file(bq_genbank_file, fileType="detect")

##
## From ptt files
##
## From a file
bq_ptt_file <- system.file('extdata/BQ.ptt', package = 'genoPlotR')
bq <- read_dna_seg_from_ptt(bq_ptt_file)
## Read directly from NCBI ftp site:
url <- "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Bartonella_henselae_Houston-1/NC_005956.ptt"
attempt <- 0
\dontrun{
while (attempt < 5){
  attempt <- attempt + 1
  bh <- try(read_dna_seg_from_ptt(url))
  if (!inherits(bh, "try-error")) {
    attempt <- 99
  } else {
    print(paste("Tried", attempt, "times, retrying in 5s"))
    Sys.sleep(5)
  }
}
}
## If attempt to connect to internet fails
if (!exists("bh")){
  data(barto)
  bh <- barto$dna_segs[[3]]
}

##
## Read from blast
##
bh_vs_bq_file <- system.file('extdata/BH_vs_BQ.blastn.tab',
                             package = 'genoPlotR')
bh_vs_bq <- read_comparison_from_blast(bh_vs_bq_file, color_scheme="grey")

## Plot
plot_gene_map(dna_segs=list(BH=bh, BQ=bq), comparisons=list(bh_vs_bq),
              xlims=list(c(1,50000), c(1, 50000)))


}
\keyword{file}
\keyword{IO}