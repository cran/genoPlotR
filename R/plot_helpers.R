################################################################################
# Plot helpers
################################################################################
# calculate arrow coordinates from gene coordinates
arrow_coord <- function(x1, x2, y1=0.5, strand=NULL, width=1, head_len=100){
  # take care of strand, to get x1 as bottom and x2 as tip of arrow
  if (!is.null(strand) && strand == -1){
    x_temp <- x2
    x2 <- x1
    x1 <- x_temp
  }
  w2 <- width/4
  # if the head of the arrow is larger than half of the gene, reduce to half
  if (head_len > abs(x1-x2)/2){
    head_len <- abs(x1-x2)/2
  }
  # calculate xi, x "internal"
  if (x2 > x1){
    xi <- x2-head_len
  } else {
    xi <- x2+head_len
  }
  list(x=c(x1,    xi,    xi,      x2, xi,      xi,    x1),
       y=c(y1-w2, y1-w2, y1-w2*2, y1, y1+w2*2, y1+w2, y1+w2)
       )
}
# coords for a block
block_coord <- function(start, end, strand, y=0.5){
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, y + strand/2, y + strand/2, y)
  list(x=x, y=y)
}

# exon coord
exon_coord <- function(start, end, strand){
  x <- c(rep(start, 2), rep(end, 2))
  if (strand == 0 ){ y <- c(0.2, 0.8, 0.8, 0.2) }
  if (strand == 1 ){ y <- c(0.5, 0.8, 0.8, 0.5) }
  if (strand == -1 ){ y <- c(0.2, 0.5, 0.5, 0.2) }
  list(x=x, y=y)
}

# coords for a zone annotation
bracket_coord <- function(start, end, y=0, w=0.1){
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, rep(y+w, 2), y)
  list(x=x, y=y)
}
# human readable coordinates
human_nt <- function(nt, signif=FALSE){
  tag <- "nt"
  mult <- 1
  med <- median(nt)
  if (med >= 1e9){
    nt <- nt/1e9
    tag <- "Gb"
    mult <- 1e9
  } else if (med >= 1e6){
    nt <- nt/1e6
    tag <- "Mb"
    mult <- 1e6
  } else if (med >= 1e3){
    nt <- nt/1e3
    tag <- "kb"
    mult <- 1e3
  }
  if (signif) nt <- signif(nt, signif)
  list(n=nt, tag=tag, mult=mult, text=paste(nt, tag))
}
# calculate comparison coordinates
calc_comp_coor <- function(gap, xlim, comp, side){
  if (length(gap) != nrow(xlim))
    stop("gap should have the same length as xlim")
  if (side < 1 && side > 2) stop("side should be 1 or 2")
  # x is the moving cursor
  x <- 0
  start <- if (side==1) comp$start1 else comp$start2 
  end <- if (side==1) comp$end1 else comp$end2
  for (i in 1:nrow(xlim)){
    # increment by the gap length
    x <- x + gap[i]
    # select comps
    idx <- start >= xlim$x0[i] & end <= xlim$x1[i]
    # re-number by substracting the xlim and adding x
    if (xlim$strand[i] == 1){
      start[idx] <- start[idx] - xlim$x0[i] + x
      end[idx] <- end[idx] - xlim$x0[i] + x
    } else {
      start[idx] <- xlim$x1[i] - start[idx] + x
      end[idx] <- xlim$x1[i] - end[idx] + x
    }
    # increment x by the length of the segment
    x <- x + xlim$length[i]
  }
  # reattribute start and stop
  if (side==1) comp$start1 <- start else comp$start2 <- start
  if (side==1) comp$end1 <- end else comp$end2 <- end
  # return the modified comp
  comp
}
middle <- function(dna_seg){
  if (!is.dna_seg(dna_seg)) stop("argument should be a dna_seg object")
  apply(dna_seg[,c("start", "end")], 1, mean)
}
