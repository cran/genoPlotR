################################################################################
# Minimize comparison sizes
################################################################################
# from a set of comparisons and lengths, determine the best possible
# arrangement of maps
minimize_comps <- function(comparisons, xlims, lengths, prel_offsets){
  # function to minimize. Calculates the mean of the absolute differences
  # between starts and ends for direct hits only
  ## mean_w_offset <- function(offset, s_ref, e_ref, s, e){
  ##   direction <- sign(e_ref-s_ref)*sign(e-s)
  ##   dists <- c(abs(s_ref-(s+offset)),
  ##            abs(e_ref-(e+offset)))[direction > 0]
  ##   mean(dists)
  ## }
  mean_w_gaps <- function(offsets, offsets_ref, xlim, xlim_ref, comp, side_ref){
    side_test <- if (side_ref == 2) 1 else 2
    # recalc for ref side
    comp <- calc_comp_coor(gap=offsets_ref, xlim=xlim_ref, comp=comp,
                           side=side_ref)
    # test
    comp <- calc_comp_coor(gap=offsets, xlim=xlim, comp=comp, side=side_test)
    direction <- sign(comp$end1-comp$start1)*sign(comp$end2-comp$start2)
    dists <- c(abs(comp$start1-comp$start2),
             abs(comp$end1-comp$end2))#[direction > 0]
    lengths <- abs(comp$end1-comp$start1) + abs(comp$end2-comp$start2)
    weighted.mean(dists, c(lengths, lengths))
  }
  n_org <- length(lengths)
  offsets <- prel_offsets
  if (length(comparisons) < 1) return(offsets)
  idx_ref <- which.max(lengths)
  max_len <- max(lengths)

  # go up from ref
  if (idx_ref > 1){
    # comp i is between org i and i+1
    for (i in (idx_ref-1):1){
      # optimise
      opt <- optim(par=offsets[[i]], fn=mean_w_gaps, method="L-BFGS-B",
                   lower=offsets[[i]],
                   ## upper=rep(max_len-lengths[i], length(offsets[[i]])),
                   offsets_ref=offsets[[i+1]],
                   xlim=xlims[[i]], xlim_ref=xlims[[i+1]],
                   comp=comparisons[[i]], side_ref=2)
      offsets[[i]] <- opt$par
    }
  }
  # go down
  if (idx_ref < n_org){
    for (i in idx_ref:(n_org-1)){
      # optimise
      opt <- optim(par=offsets[[i+1]], fn=mean_w_gaps, method="L-BFGS-B",
                   lower=offsets[[i+1]],
                   ## upper=rep(max_len-lengths[i+1], length(offsets[[i+1]])),
                   offsets_ref=offsets[[i]],
                   xlim=xlims[[i+1]], xlim_ref=xlims[[i]],
                   comp=comparisons[[i]], side_ref=1)
      offsets[[i+1]] <- opt$par
    }
  }
  offsets
}
