################################################################################
# File reading functions: read mauve backbone
################################################################################
read_mauve_backbone <- function(file, ref=1, gene_type="side_blocks",
                                header=TRUE, filter_low=0, ...){
  blocks <- read.table(file, stringsAsFactors=FALSE, h=header)
  n_orgs <- ncol(blocks)/2
  n_blocks <- nrow(blocks)
  # ordering from ref: 
  for (i in 1:n_blocks){
    if (blocks[i,ref*2] < 0) blocks[i,] <- -blocks[i,] 
  }
  if (any(blocks[,ref*2] < 0))
      stop("Not all rows in ref columns are positive. Contact author.")
  blocks <- blocks[order(blocks[,ref*2]),]
  # filter, if needed
  if (is.numeric(filter_low) && filter_low > 1){
    sizes <- matrix(NA, nrow=n_blocks, ncol=n_orgs)
    for (i in 1:n_orgs){
      sizes[,i] <- abs(blocks[,2*i]) - abs(blocks[,2*i-1])
    }
    min_sizes <- apply(sizes, 1, min, na.rm=TRUE)
    blocks <- blocks[min_sizes >= filter_low,]
    n_blocks <- nrow(blocks)
  }
  # prepare objects
  dna_segs <- list()
  comparisons <- list()
  # colors: rainbow starting from the beginning
  col <- rainbow(n=n_blocks)
  # run through organisms
  for (i in 1:n_orgs){
    # prepare dna_seg
    s <- blocks[,2*i-1]
    e <- blocks[,2*i]
    strand <- sign(s)
    df <- data.frame(name=paste("block", 1:n_blocks, sep="_"),
                     start=abs(s), end=abs(e), strand=strand,
                     stringsAsFactors=FALSE)
    df$col <- col
    dna_segs[[i]] <- as.dna_seg(df[df$strand != 0,], gene_type=gene_type, ...)
    # prepare comparison (not with i=1)
    if (i > 1){
      s0 <- blocks[,2*(i-1)-1]
      e0 <- blocks[,2*(i-1)]
      strand0 <- sign(s0)
      df <- data.frame(start1=ifelse(strand0 == 1, abs(s0), abs(e0)),
                       end1=ifelse(strand0 == 1, abs(e0), abs(s0)),
                       start2=ifelse(strand == 1, abs(s), abs(e)),
                       end2=ifelse(strand == 1, abs(e), abs(s)))
      comparison <- as.comparison(df[df$start1 != 0 & df$start2 != 0,])
      # apply red_blue color scheme
      comparison$col <- apply_color_scheme(x=NULL,
                                           direction=comparison$direction,
                                           color_scheme="red_blue")
      comparisons[[i-1]] <- comparison
    }
  }
  list(dna_segs=dna_segs, comparisons=comparisons)
}
