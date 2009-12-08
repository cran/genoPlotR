################################################################################
# Phylo tree
################################################################################
# a function to plot phylogenetic trees ported to grid from ade4.
# to be rewritten...
phylog_grob <-
  function (x, y=NULL, f.phylog = 1, cleaves = 0, cnodes = 0, 
            labels.leaves = names(x$leaves), clabel.leaves = 1,
            labels.nodes = names(x$nodes), 
            clabel.nodes = 0, sub = "", csub = 1.25, possub = "bottomleft", 
            draw.box = FALSE, ...) 
{
  if (!inherits(x, "phylog")) 
    stop("Non convenient data")
  leaves.number <- length(x$leaves)
  leaves.names <- names(x$leaves)
  nodes.number <- length(x$nodes)
  nodes.names <- names(x$nodes)
  if (length(labels.leaves) != leaves.number) 
    labels.leaves <- names(x$leaves)
  if (length(labels.nodes) != nodes.number) 
    labels.nodes <- names(x$nodes)
  leaves.car <- gsub("[_]", " ", labels.leaves)
  nodes.car <- gsub("[_]", " ", labels.nodes)
  if (f.phylog < 0.05) 
    f.phylog <- 0.05
  if (f.phylog > 0.95) 
    f.phylog <- 0.95
  maxx <- max(x$droot)
# plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
#     yaxt = "n", xlim = c(-maxx * 0.15, maxx/f.phylog), ylim = c(-0.05, 
#         1), xaxs = "i", yaxs = "i", frame.plot = FALSE)
  
  x.leaves <- x$droot[leaves.names]
  x.nodes <- x$droot[nodes.names]
  if (!is.null(y)){
    # check that it is constrained between 0 and 1
    if (!all(range(y) %in% c(0,1)))
      stop("y should be constrained between 0 and 1")
  } else {
    y <- seq(1, 0, len=leaves.number)
  }
  names(y) <- leaves.names
  xcar <- maxx * 1.05
  xx <- c(x.leaves, x.nodes)
  # prepare grobs
  labelGrobs <- gList()
  labelSegGrobs <- gList()
  branchesGrobs <- gList()
  # leaves labels and segments leading to it
  if (clabel.leaves > 0) {
    for (i in 1:leaves.number) {
      labelGrobs[[i]] <-
        textGrob(x=0, y=y[i], label=leaves.car[i], just="left",
                 gp=gpar(cex=clabel.leaves), default.units="native",
                 name=paste("tree.leave.label", i, sep="."))
      labelSegGrobs[[i]] <-
        segmentsGrob(xcar, y[i], xx[i], y[i],
                     gp=gpar(col=grey(0.7)), default.units="native",
                     name=paste("tree.leave.segment", i, sep="."))
    }
  }
  yleaves <- y[1:leaves.number]
  xleaves <- xx[1:leaves.number]
# if (cleaves > 0) {
#     for (i in 1:leaves.number) {
#         points(xx[i], y[i], pch = 21, bg = 1, cex = par("cex") * 
#             cleaves)
#     }
# }
  yn <- rep(0, nodes.number)
  names(yn) <- nodes.names
  y <- c(y, yn)
  # plot tree
  for (i in 1:length(x$parts)) {
    w <- x$parts[[i]]
    but <- names(x$parts)[i]
    y[but] <- mean(y[w])
    b <- range(y[w])
    # vertical branches
    branchesGrobs[[i]] <-
      segmentsGrob(xx[but], b[1], xx[but], b[2], default.units="native",
                   name=paste("tree.branch.vert.segment", i, sep="."))
    x1 <- xx[w]
    y1 <- y[w]
    x2 <- rep(xx[but], length(w))
    # horizontal branches
    branchesGrobs[[i+length(x$parts)]] <-
      segmentsGrob(x1, y1, x2, y1, default.units="native",
                   name=paste("tree.branch.horiz.segment", i, sep="."))
}
# if (cnodes > 0) {
#     for (i in nodes.names) {
#         points(xx[i], y[i], pch = 21, bg = "white", cex = cnodes)
#     }
# }
# if (clabel.nodes > 0) {
#     scatterutil.eti(xx[names(x.nodes)], y[names(x.nodes)], 
#         nodes.car, clabel.nodes)
# }
  x <- x.leaves
  y <- y[leaves.names]
  xbase <- xcar
# if (csub > 0) 
#     scatterutil.sub(sub, csub = csub, possub = possub)
# if (draw.box) 
#     box()
# if (cleaves > 0) 
#     points(xleaves, yleaves, pch = 21, bg = 1, cex = par("cex") * 
#         cleaves)
  # creating gTree for branches
  branchesTree <- gTree(children=gList(branchesGrobs, labelSegGrobs),
                        vp=viewport(xscale=c(0, xcar), yscale=c(0, 1),
                          name="tree.branches"),
                        name="tree.branchesTree")
  labelTree <- gTree(children=labelGrobs,
                     vp=viewport(xscale=c(0, 1), yscale=c(0, 1),
                       name="tree.labels"),
                     name="tree.labelsTree")
  # finally plotting
  label_width <- unit(1, "grobwidth",
                      labelGrobs[[which.max(nchar(leaves.car))]])
  layout <- grid.layout(1, 2, widths=unit.c(unit(1, "null"), label_width))
  fg <- frameGrob(layout=layout, name="treeFrameGrob",
                  vp=viewport(name="treeFrame"))
  fg <- placeGrob(fg, branchesTree, col=1)
  fg <- placeGrob(fg, labelTree, col=2)
  return(invisible(list(xy = data.frame(x = x, y = y), xbase = xbase, 
                        cleaves = cleaves, grob=fg, width=label_width)))
}


# permute tree leaves to match labels
permute_tree <- function(tree, labels){
  ref <- names(tree$leaves)
  n <- length(ref)
  wanted_permut <- rep(NA, n)
  for (i in 1:n){
    idx <- which(ref[i] == labels)
    if(length(idx) != 1) stop("Non-unique or non-matching labels")
    wanted_permut[i] <- idx
  }
  permuts <- enum.phylog(tree)
  equals <- apply(permuts, 1, function(x)
                  identical(as.numeric(x), as.numeric(wanted_permut)))
  if (!any(equals))
    stop("No tree permutation compatible with label order. Change input order")
  if (!sum(equals))
    stop("Several permutations matching. Something went wrong, contact the author")
  return(permuts[equals,])
}
