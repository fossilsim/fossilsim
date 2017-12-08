#' Plot simulated fossils
#'
#' @description
#' This function is adapted from the \emph{ape} function \code{plot.phylo} used to plot phylogenetic trees.
#' The function can be used to plot simulated fossils (\code{show.fossils = TRUE}), with or without the corresponding tree (\code{show.tree = TRUE}),
#' stratigraphic intervals (\code{show.strata = TRUE}), stratigraphic ranges (\code{show.ranges = TRUE}) and sampling proxy data (\code{show.proxy = TRUE}).
#' Interval ages can be specified as a vector (\code{interval.ages}) or a uniform set of interval ages can be specified using the
#' number of intervals (\code{strata}) and maximum interval age (\code{max}), where interval length \eqn{= basin.age/strata}.
#' If no maximum age is specified, the function calculates a maximum interval age slightly older than the root edge (or root age if \code{root.edge = FALSE}),
#' using the function \code{basin.age}.
#'
#' @param fossils Fossils object.
#' @param tree Phylo object.
#' @param show.fossils If TRUE plot fossils (default = TRUE).
#' @param show.tree If TRUE plot the tree  (default = TRUE).
#' @param show.ranges If TRUE plot stratigraphic ranges (default = FALSE).
#' @param show.strata If TRUE plot strata  (default = FALSE).
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param strata Number of stratigraphic intervals.
#' @param max Maximum age of a set of equal length intervals. If no value is specified (max = NULL), the function uses a maximum age based on tree height.
#' @param show.axis If TRUE plot x-axis (default = TRUE).
#' @param binned If TRUE fossils are plotted at the mid point of each interval.
#' @param show.proxy If TRUE add water depth profile (default = FALSE).
#' @param proxy.data Vector of sampling proxy data (default = NULL).
#' @param show.preferred.environ If TRUE add species prefferred environmental value (e.g. water depth) (default = FALSE).
#' @param preferred.environ Prefferred environmental value (e.g. water depth).
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @param hide.edge If TRUE hide the root edge but still incorporate it into the automatic timescale (default = FALSE).
#' @param edge.width A numeric vector giving the width of the branches of the plotted phylogeny. These are taken to be in the same order as the component edge of \code{tree}. If fewer widths are given than the number of edges, then the values are recycled.
#' @param show.tip.label Whether to show the tip labels on the phylogeny (defaults to FALSE).
#' @param align.tip.label A logical value or an integer. If TRUE, the tips are aligned and dotted lines are drawn between the tips of the tree and the labels. If an integer, the tips are aligned and this gives the type of the lines (lty).
#' @param fcex Numeric value giving the factor used to scale the points representing the fossils. Only used if \code{show.fossils = TRUE}.
#' @param fcol Color of fossil occurrences or ranges.
#' @param ecol Color of extant samples. This only works if binned = FALSE.
#' @param ... Additional parameters to be passed to \code{plot.default}.
#'
#' @examples
#' set.seed(123)
#'
#' ## simulate tree
#' t <- TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]
#'
#' ## simulate fossils under a Poisson sampling process
#' f <- sim.fossils.poisson(rate = 3, tree = t)
#' plot(f, t)
#' # add a set of equal length strata
#' plot(f, t, show.strata = TRUE, strata = 4)
#'
#' ## simulate fossils under a non-uniform model of preservation
#' # assign a max interval based on tree height
#' max = basin.age(t)
#' times = c(0, 0.3, 1, max)
#' rates = c(4, 1, 0.1)
#' f <- sim.fossils.intervals(t, interval.ages = times, rates = rates)
#' plot(f, t, show.strata = TRUE, interval.ages = times)
#' # add proxy data
#' plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = TRUE, proxy.data = rates)
#'
#' @export
#' @importFrom graphics par points lines text axis mtext segments
plot.fossils<-function(fossils, tree, show.fossils = TRUE, show.tree = TRUE, show.ranges = FALSE,
                       # age info/options
                       show.strata = FALSE, strata = 1, max = NULL, interval.ages = NULL, binned = FALSE, show.axis = TRUE,
                       # proxy stuff
                       show.proxy = FALSE, proxy.data = NULL,
                       show.preferred.environ = FALSE, preferred.environ = NULL,
                       # tree appearance
                       root.edge = TRUE, hide.edge = FALSE, edge.width = 1, show.tip.label = FALSE, align.tip.label = FALSE,
                       # fossil appearance
                       fcex = 1.2, fcol = "darkorange", ecol = NULL, ...) {

  # TODO probably not appropriate all the time
  fossils$h = (fossils$hmin + fossils$hmax)/2

  x<-tree  # tree
  ba<-max

  if(!show.tree)
    align.tip.label = TRUE

  # other possible options
  show.node.label = FALSE # this doesn't do anything
  edge.color = "black"
  edge.lty = 1
  font = 3
  cex = par("cex")
  tip.color = "black"
  label.offset = 0
  underscore = FALSE
  plot = TRUE
  node.depth = 1
  no.margin = FALSE
  x.lim = NULL
  y.lim = NULL
  adj = NULL
  srt = 0

  if(is.null(ecol))
    ecol = fcol

  if(!(is.fossils(fossils)))
    stop("fossils must be an object of class \"fossils\"")

  # check the tree
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  if (any(tabulate(x$edge[, 1]) == 1))
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")

  # check edge lengths != null
  # check tree is rooted
  # check tree is binary

  if(is.null(x$root.edge))
    root.edge = FALSE

  if(show.strata || show.proxy){
    if( (is.null(interval.ages)) && (is.null(strata)) )
      stop("To plot interval info specify interval.ages OR number of strata, else use show.strata = FALSE")
  }

  if(show.proxy && is.null(proxy.data))
    stop("Specify sampling profile")

  # is there a more efficient way of doing this?
  if(show.proxy){
    if(!is.null(interval.ages)){
      if( (length(interval.ages) - 1) != length(proxy.data) )
        stop("Make sure number of sampling proxy data points matches the number of intervals")
    } else {
      if(strata != length(proxy.data))
        stop("Make sure number of sampling proxy data points matches the number of intervals")
    }
    if(any(is.na(proxy.data)))
      stop("I can't handle NA proxy values right now, please use 0 for the time being")
  }

  # required ape C fxns
  .nodeHeight <- function(edge, Nedge, yy) .C(ape::node_height, as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(ape::node_depth, as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, edge.length) .C(ape::node_depth_edgelength, as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), double(Ntip + Nnode))[[5]]

  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1

  type <- "phylogram"
  direction <- "rightwards"

  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  } else {
    if (align.tip.label)
      align.tip.label.lty <- 3
  }

  phyloORclado <- TRUE # = "phylogram"
  horizontal <- TRUE # = "rightwards"

  xe <- x$edge # used in the last part of the fxn
  yy <- numeric(Ntip + Nnode)
  TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
  yy[TIPS] <- 1:Ntip
  z <- stats::reorder(x, order = "postorder")

  #yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
  yy <- .nodeHeight(z$edge, Nedge, yy)
  xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, z$edge.length)

  if (root.edge) {
    xx <- xx + x$root.edge
  }

  #if (no.margin)
  # par(mai = rep(0, 4))
  if (show.tip.label)
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  if (is.null(x.lim)) {
    x.lim <- c(0, NA)
    pin1 <- par("pin")[1]
    strWi <- graphics::strwidth(x$tip.label, "inches", cex = cex)
    xx.tips <- xx[1:Ntip] * 1.04
    alp <- try(stats::uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
    if (is.character(alp)) {
      tmp <- max(xx.tips)
      if (show.tip.label)
        tmp <- tmp * 1.5
    }
    else {
      tmp <- if (show.tip.label)
        max(xx.tips + strWi/alp)
      else max(xx.tips)
    }
    if (show.tip.label)
      tmp <- tmp + label.offset
    x.lim[2] <- tmp
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
  }
  if (is.null(y.lim)) {
    y.lim <- c(1, Ntip)
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    y.lim[1] <- 1
  }

  # the order in which you use plot and par here is very important
  if(show.proxy){
    old.par = par("usr", "mar", "oma", "xpd", "mgp","fig")
    par(mar=c(1, 3.5, 0, 0.5)) # to reduce the margins around each plot - bottom, left, top, right -- this is harder to manipulate
    par(oma=c(2, 0, 2, 0)) # to add an outer margin to the top and bottom of the graph -- bottom, left, top, right
    par(xpd=NA) # allow content to protrude into outer margin (and beyond)
    par(mgp=c(1.5, .5, 0)) # to reduce the spacing between the figure plotting region and the axis labels -- axis label at 1.5 rows distance, tick labels at .5 row
    par(fig=c(0,1,0,0.4))
    graphics::plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", axes = FALSE, asp = NA, ...)
    par(fig=c(0,1,0.4,1))
  }
  else{
    old.par = par("xpd")
    par(xpd=NA)
    # open a new plot window
    graphics::plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", axes = FALSE, asp = NA, ...)
  }

  if (plot) {

    if(show.strata || show.axis){
      if( (is.null(interval.ages)) ){
        if(is.null(ba))
          ba = basin.age(x, root.edge = root.edge)
        s1 = ba / strata # horizon length (= max age of youngest horizon)
        horizons.max = seq(s1, ba, length = strata)
        horizons.min = horizons.max - s1
        s1 = rev(horizons.max - horizons.min) # rev is unneccessary here
      } else {
        horizons.min = utils::head(interval.ages, -1)
        horizons.max = interval.ages[-1]
        ba = max(horizons.max)
        s1 = rev(horizons.max - horizons.min)
        strata = length(horizons.max)
      }
    }

    # add colored strata
    # rect(xleft, ybottom, xright, ytop)
    if(show.strata || show.axis){

      # y-axis:
      #y.bottom = 0
      #y.top = max(yy)+1
      y.bottom = 0.5
      y.top = max(yy) + 0.5

      # x-axis:
      #s1 = ba / strata
      x.left = 0 - (ba - max(xx))
      x.right = x.left + s1[1]
      cc = 1 # color switch

      axis.strata = x.left

      for(i in 1:strata){
        if(cc %% 2 == 0)
          col="grey90"
        else
          col="grey95"
        if(show.strata)
          graphics::rect(xleft = x.left, xright = x.right, ybottom = y.bottom, ytop = y.top, col=col, border=NA)
        x.left = x.right
        x.right = x.left + s1[i+1]
        cc = cc + 1
        axis.strata = c(axis.strata, x.left)
      }

      x.labs = rev(round(c(0, horizons.max),2))

      if(show.proxy)
        labs = FALSE
      else
        labs = x.labs

      if(show.axis)
        axis(1, col = 'grey75', at = axis.strata, labels = labs, lwd = 2, line = 0.5, col.axis = 'grey75', cex.axis = .7)
    }

    # plot the tree
    if(show.tree)
      ape::phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, edge.color, edge.width, edge.lty)

    # format the root edge
    if (root.edge && show.tree && !hide.edge) {
      rootcol <- if (length(edge.color) == 1)
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1)
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1)
        edge.lty
      else 1

      # plot the root edge
      segments(0, yy[ROOT], x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, lty = rootlty)

    }

    # format tip labels
    if (show.tip.label) {

      # x and y co-ordinates
      adj = 0
      MAXSTRING <- max(graphics::strwidth(x$tip.label, cex = cex))
      loy <- 0
      lox <- label.offset + MAXSTRING * 1.05 * adj

      if (is.expression(x$tip.label))
        underscore <- TRUE
      if (!underscore)
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- max(xx[1:Ntip])
          yy.tmp <- yy[1:Ntip]
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label,adj = adj, font = font, srt = srt, cex = cex, col = tip.color)
      }
    }

    # add node labels
    if (show.node.label){
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
           x$node.label, adj = adj, font = font, srt = srt,
           cex = cex)
    }

    # fossils
    if(show.fossils){
      if(binned){
        # TODO check correction
        #if(attr(fossils, "ages") == "interval.max"){
        #y = sapply(fossils$h, function(x) max(which(horizons.max == x)) )
        #} else {
        # treat fossil data as continuous
        y = sapply(fossils$hmin, function(x) max(which(horizons.min <= x)) )
        #}
        points(max(xx) - horizons.max[y] + (rev(s1)[y]/2), yy[fossils$edge] , col = fcol, pch = 19, cex = fcex)
      }
      else{
        if(fcol == ecol)
          points(max(xx) - fossils$h, yy[fossils$edge], col = fcol, pch = 19, cex = fcex)
        else{
          extant = fossils[which(fossils$h == 0), ]
          extinct = fossils[which(fossils$h != 0), ]
          points(max(xx) - extant$h, yy[extant$edge], col = ecol, pch = 19, cex = fcex)
          points(max(xx) - extinct$h, yy[extinct$edge], col = fcol, pch = 19, cex = fcex)
        }
      }
    }

    # ranges
    if(show.ranges){
      buffer = 0.01 * max(xx) # buffer for singletons

      if(binned){
        if(attr(fossils, "ages") == "interval.max"){
          y = sapply(fossils$h, function(x) max(which(horizons.max == x)) )
        } else {
          # treat fossil data as continuous
          y = sapply(fossils$h, function(x) max(which(horizons.min <= x)) )
        }
        fossils$r = max(xx) - horizons.max[y] + (rev(s1)[y]/2)
      } else {
        fossils$r = max(xx) - fossils$h
      }

      for(i in unique(fossils$edge)) {

        range = fossils$r[which(fossils$edge==i)]

        if(length(range) == 1)
          range = c(range-buffer,range+buffer)

        species = rep(yy[i], length(range))

        lines(y = species, x = range, lwd = 6, col=fcol)

      }
      fossils$r = NULL
    }

    # water depth profile
    add.depth.axis = TRUE
    if(show.proxy){
      add.depth.profile(proxy.data, axis.strata, strata, show.axis, add.depth.axis,
                        show.preferred.depth = show.preferred.environ, PD = preferred.environ, x.labs = x.labs)
    }
  }

  par(old.par)
  L <- list(type = type, use.edge.length = TRUE,
            node.pos = NULL, node.depth = node.depth, show.tip.label = show.tip.label,
            show.node.label = show.node.label, font = font, cex = cex,
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
            x.lim = x.lim, y.lim = y.lim, direction = direction,
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time,
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
         envir = ape::.PlotPhyloEnv)
  invisible(L)
}

add.depth.profile<-function(depth.profile,axis.strata,strata,show.axis,add.depth.axis,show.preferred.depth=TRUE,PD=NULL,x.labs=FALSE){
  par(fig=c(0,1,0,0.4), new = T)
  # change the y-axis scale for depth
  u = par("usr") # current scale
  depth.profile = rev(depth.profile)
  depth = depth.profile # there is some redundancy here
  tol = max(depth) * 0.1
  par(usr = c(u[1], u[2], min(depth) - tol, max(depth) + tol))
  # define x-axis values (time)
  #time = axis.strata[1:strata] + ((axis.strata[2] - axis.strata[1])/2)
  time = c()
  depth = c()
  for(i in 1:length(depth.profile)){
    if(i == 1)
      time = c(time, axis.strata[i]+0.01, axis.strata[i+1])
    else
      time = c(time, axis.strata[i], axis.strata[i+1])
    depth = c(depth, depth.profile[i], depth.profile[i])
  }
  # plot proxy
  lines(time,depth)
  if(show.preferred.depth)
    lines(x = axis.strata, y = rep(PD, length(axis.strata)), col = "grey75", lwd = 2, lty = 3)
  if(show.axis){
    axis(1, col = 'grey75', at = axis.strata, labels = x.labs, lwd = 2, col.axis = 'grey75', cex.axis = .7)
    mtext(1, col = 'grey75', text="Time before present", line = 1.5)
  }
  if(add.depth.axis){
    axis(2, col = 'grey75', labels = TRUE, lwd = 2, las = 2, col.axis = 'grey75', line = 0.5, cex.axis = .7)
    mtext(2, col = 'grey75', text="Sampling proxy", line = 2)
  }
}
