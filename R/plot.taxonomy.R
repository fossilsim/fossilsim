#' Plot simulated taxonomy
#'
#' @description
#' This function is adapted from the \emph{ape} function \code{plot.phylo} used to plot phylogenetic trees.
#' The function can be used to plot simulated taxonomy along with the corresponding tree.
#'
#' @param x Taxonomy object.
#' @param tree Phylo object.
#' @param show.mode Indicate speciation mode.
#' @param show.legend Add a legend for the symbols indicating different speciation modes.
#' @param legend.position Position of the legend. Options include \code{"bottomleft"} (default), \code{"topleft"}, \code{"bottomright"}, \code{"topright"} or a vector of length 2 with x, y coordinates.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @param hide.edge If TRUE hide the root edge but still incorporate it into the automatic timescale (default = FALSE).
#' @param edge.width A numeric vector giving the width of the branches of the plotted phylogeny. These are taken to be in the same order as the component edge of \code{tree}. If fewer widths are given than the number of edges, then the values are recycled.
#' @param show.tip.label Whether to show the tip labels on the phylogeny (defaults to FALSE).
#' @param align.tip.label A logical value or an integer. If TRUE, the tips are aligned and dotted lines are drawn between the tips of the tree and the labels. If an integer, the tips are aligned and this gives the type of the lines (lty).
#' @param taxa.palette Colour palette used for taxa. Colours are assigned to taxa using the function \code{grDevices::hcl.colors()} and the default palette is "viridis". Other colour blind friendly palettes include \code{"Blue-Red 3"} and \code{"Green-Brown"}.
#' @param cex Numeric value giving the factor used to scale the points representing the fossils when \code{show.fossils = TRUE}.
#' @param ... Additional parameters to be passed to \code{plot.default}.
#'
#' @examples
#' set.seed(123)
#'
#' ## simulate tree
#' t = TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]
#'
#' ## simulate taxonomy
#' s = sim.taxonomy(t, 0.5, 1)
#'
#' ## plot the output
#' plot(s, t)
#'
#' @export
#' @importFrom graphics legend plot
plot.taxonomy = function(x, tree, show.mode = TRUE, show.legend = TRUE, legend.position = "bottomleft", root.edge = TRUE, hide.edge = FALSE, edge.width = 1, show.tip.label = FALSE, align.tip.label = FALSE, taxa.palette = "viridis", cex = 1.2,...){

  taxonomy = x

  # hard coded options for tree appearance
  edge.color = "black"
  edge.lty = 1
  font = 3 # italic
  tip.color = "black"
  label.offset = 0.02
  align.tip.label.lty = 3
  underscore = FALSE
  plot = TRUE
  node.depth = 1
  no.margin = FALSE
  adj = NULL
  srt = 0

  if(!"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  # tolerance for extant tips and interval/ fossil age comparisons #NEEDED
  tol = min((min(tree$edge.length)/100), 1e-8)

  if(is.null(tree$root.edge))
    root.edge = FALSE

  max.age = NULL
  if(is.null(max.age))
    ba = tree.max(tree, root.edge = root.edge)
  else ba = max.age

  # check the tree
  Ntip <- length(tree$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }

  if (any(tabulate(tree$edge[, 1]) == 1))
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")

  if(!ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(tree$edge.length))
    stop("tree must have edge lengths")

  Nnode <- tree$Nnode
  if (any(tree$edge < 1) || any(tree$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1

  if(!all(tree$edge %in% taxonomy$edge))
    stop("Mismatch between tree and taxonomy objects")

  # collect data for plotting the tree
  type = "phylogram"
  direction = "rightwards"
  horizontal = TRUE # = "rightwards"

  xe = tree$edge # used in the last part of the fxn
  yy = numeric(Ntip + Nnode)
  TIPS = tree$edge[tree$edge[, 2] <= Ntip, 2]
  yy[TIPS] = 1:Ntip

  yy = ape::node.height(tree)
  xx = ape::node.depth.edgelength(tree)

  if(root.edge)
    xx = xx + tree$root.edge

  x.lim = NULL
  y.lim = NULL

  # x.lim is defined by max(ba, tree age)
  if(ba > max(xx))
    xx = xx + (ba - max(xx))

  if (is.null(x.lim)) {
    x.lim = c(0, NA)
    pin1 = par("pin")[1]
    strWi = graphics::strwidth(tree$tip.label, "inches", cex = par("cex"))
    xx.tips = xx[1:Ntip] * 1.04
    alp = try(stats::uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
    if (is.character(alp)) {
      tmp = max(xx.tips)
      if (show.tip.label)
        tmp = tmp * 1.5
    } else {
      tmp = if (show.tip.label)
        max(xx.tips + strWi/alp)
      else max(xx.tips)
    }
    if (show.tip.label)
      tmp = tmp + label.offset
    x.lim[2] = tmp
  }
  else if (length(x.lim) == 1) {
    x.lim = c(0, x.lim)
  }

  if (is.null(y.lim)) {
    y.lim = c(1, Ntip)
  }
  else if (length(y.lim) == 1) {
    y.lim = c(0, y.lim)
    y.lim[1] = 1
  }

  # the order in which you use plot and par here is very important
  old.par = list(xpd = par("xpd"))
  par(xpd=NA)
  # open a new plot window
  graphics::plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", axes = FALSE, asp = NA, ...)

  taxonomy$col = "red"

  # taxonomy colours
  sps = unique(taxonomy$sp)
  col = grDevices::hcl.colors(length(sps), palette = taxa.palette)
  j = 0
  for(i in sps){
    j = j + 1
    taxonomy$col[which(taxonomy$sp == i)] = col[j]
  }

  if (plot) {

    # plot the tree
    ape::phylogram.plot(tree$edge, Ntip, Nnode, xx, yy, horizontal, edge.color, edge.width, edge.lty)

    # format the root edge
    if (root.edge && !hide.edge) {
      rootcol = if(length(edge.color) == 1)
        edge.color
      else "black"
      rootw = if(length(edge.width) == 1)
        edge.width
      else 1
      rootlty = if(length(edge.lty) == 1)
        edge.lty
      else 1

      # plot the root edge
      segments(xx[ROOT] - tree$root.edge, yy[ROOT], xx[ROOT], yy[ROOT], col = rootcol, lwd = rootw, lty = rootlty)
    }

    # format tip labels
    if (show.tip.label) {

      # x and y co-ordinates
      adj = 0
      MAXSTRING <- max(graphics::strwidth(tree$tip.label, cex = cex))
      loy <- 0
      lox <- label.offset + MAXSTRING * 1.05 * adj

      if (is.expression(tree$tip.label))
        underscore <- TRUE
      if (!underscore)
        tree$tip.label <- gsub("_", " ", tree$tip.label)

      if (align.tip.label) {
        xx.tmp <- max(xx[1:Ntip])
        yy.tmp <- yy[1:Ntip]
        segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, lty = align.tip.label.lty)
      }
      else {
        xx.tmp <- xx[1:Ntip]
        yy.tmp <- yy[1:Ntip]
      }
      text(xx.tmp + lox, yy.tmp + loy, tree$tip.label,adj = adj, font = font, srt = srt, cex = cex, col = tip.color)
    }

    # for show.taxonomy = TRUE
    # fetch oldest and youngest edge associated with a species
    # deal with the oldest part of the ranges
    # deal with the youngest part of the ranges
    # everyting inbetween

    for(i in unique(taxonomy$sp)){

      mn = min(taxonomy$end[which(taxonomy$sp == i)])
      mx = max(taxonomy$start[which(taxonomy$sp == i)])
      edge.mn = taxonomy$edge[which(taxonomy$sp == i & taxonomy$end == mn)]
      edge.mx = taxonomy$edge[which(taxonomy$sp == i & taxonomy$start == mx)][1]

      edges = find.edges.inbetween(edge.mn, edge.mx, tree)

      col = taxonomy$col[which(taxonomy$sp == i)]

      for(j in edges){

        # single edge
        range = max(xx) - c(taxonomy$start[which(taxonomy$edge == j & taxonomy$sp == i)],
                            taxonomy$end[which(taxonomy$edge == j
                                               & taxonomy$sp == i)])
        # plot ranges
        sp = yy[j]

        w = 0.1

        # plot ranges
        rect(min(range), sp+w, max(range), sp-w, col=adjustcolor(col, alpha.f = 0.5))

        # for speciation mode
        if(j == edge.mx && show.mode) {
          mode =  taxonomy$mode[which(taxonomy$edge == edge.mx & taxonomy$sp == i)]

          if(mode == "r" || mode == "o") next
          # plot mode
          if(mode == "a") pch = 17
          else if(mode == "b") pch = 15
          else if(mode == "s") pch = 16

          points(min(range), sp, cex = 1.2, pch = pch, col = 1, bg = 1)
        }
      }
    }
    # legend
    if(show.legend) {
      x = legend.position[1]
      y = if(length(legend.position) > 1) legend.position[2] else NULL
      legend(x = x, y = y, legend=c("Budding", "Bifurcation", "Anagenesis"), pch = c(15, 16, 17))
    }
  }
  par(old.par)
  L <- list(type = type, use.edge.length = TRUE,
            node.pos = NULL, node.depth = node.depth, show.tip.label = show.tip.label,
            show.node.label = FALSE, font = font, cex = cex,
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
            x.lim = x.lim, y.lim = y.lim, direction = direction,
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = tree$root.time,
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
         envir = ape::.PlotPhyloEnv)
  invisible(L)
}
