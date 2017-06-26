#' Plot simulated fossils
#'
#' @param fossils Dataframe of sampled fossils (sp = edge labels, h = ages).
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
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @param hide.edge If TRUE hide the root edge but still incorporate it into the automatic timescale (default = FALSE).
#' @param fcol Color of fossil occurrences or ranges.
#'
#' @examples
#' set.seed(123)
#' # simulate tree
#' t<-TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]
#'
#' # simulate fossils under a Poisson sampling process
#' f<-sim.fossils.poisson(t,3)
#' plot(f, t)
#' # add a set of equal length strata
#' plot(f, t, show.strata = TRUE, strata = 4)
#'
#' # simulate fossils under an alternative non-uniform model of preservation
#' # assign a max interval based on tree height
#' max = basin.age(t)
#' times = c(0, 0.3, 1, max)
#' rates = c(4, 1, 0.1)
#' f<-sim.fossils.non.unif(t, times, rates)
#' plot(f, t, show.strata = TRUE, interval.ages = times)
#' # add proxy data
#' plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = T, proxy.data = rates)
#'
#' @export
plot.fossils<-function(fossils, tree, show.fossils = TRUE, show.tree = TRUE, show.ranges = FALSE,
                       show.strata = FALSE, strata = 1, max = NULL, interval.ages = NULL,
                       show.axis = TRUE, show.proxy = FALSE, proxy.data = NULL, binned = FALSE,
                       root.edge = TRUE, hide.edge = FALSE, show.tip.label = FALSE, align.tip.label = FALSE, fcol = "darkorange",
                       fcex = 1.2, edge.width = 1, ...) {
  fossils<-fossils
  x<-tree  # tree
  show.fossils<-show.fossils
  show.tree<-show.tree
  show.ranges<-show.ranges
  # age info/options
  show.strata<-show.strata
  interval.ages<-interval.ages
  strata<-strata
  ba<-max
  show.axis<-show.axis
  binned<-binned
  # proxy stuff
  show.profile<-show.proxy
  depth.profile<-proxy.data
  # tree appearance
  root.edge<-root.edge
  hide.edge<-hide.edge
  show.tip.label<-show.tip.label
  align.tip.label<-align.tip.label
  edge.width<-edge.width
  # fossil appearance
  fcol<-fcol
  fcex<-fcex

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

  if(show.strata){
    if( (is.null(interval.ages)) && (is.null(strata)) )
      stop("To plot intervals specify interval.ages OR number of strata, else use show.strata = FALSE")
  }

  # check edge lengths != null
  # check tree is rooted
  # check tree is binary

  if(is.null(x$root.edge))
    root.edge = FALSE

  if(show.profile && is.null(depth.profile))
    stop("Specify sampling profile")

  # required C fxns
  .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C(ape::node_height,as.integer(Ntip), as.integer(Nnode), as.integer(edge[,1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy))[[6]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(ape::node_depth,as.integer(Ntip), as.integer(Nnode), as.integer(edge[,1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[6]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, edge.length) .C(ape::node_depth_edgelength, as.integer(Ntip),as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[,2]), as.integer(Nedge), as.double(edge.length), double(Ntip + Nnode))[[7]]

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
  z <- reorder(x, order = "postorder")

  yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
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
    strWi <- strwidth(x$tip.label, "inches", cex = cex)
    xx.tips <- xx[1:Ntip] * 1.04
    alp <- try(uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
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
  if(show.profile){
    old.par = par("usr", "mar", "oma", "xpd", "mgp","fig")
    par(mar=c(1, 3.5, 0, 0.5)) # to reduce the margins around each plot - bottom, left, top, right -- this is harder to manipulate
    par(oma=c(2, 0, 2, 0)) # to add an outer margin to the top and bottom of the graph -- bottom, left, top, right
    par(xpd=NA) # allow content to protrude into outer margin (and beyond)
    par(mgp=c(1.5, .5, 0)) # to reduce the spacing between the figure plotting region and the axis labels -- axis label at 1.5 rows distance, tick labels at .5 row
    par(fig=c(0,1,0,0.4))
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", axes = FALSE, asp = NA, ...)
    par(fig=c(0,1,0.4,1))
  }
  else{
    old.par = par("xpd")
    par(xpd=NA)
    # open a new plot window
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", axes = FALSE, asp = NA, ...)
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
        horizons.min = head(interval.ages, -1)
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
          rect(xleft = x.left, xright = x.right, ybottom = y.bottom, ytop = y.top, col=col, border=NA)
        x.left = x.right
        x.right = x.left + s1[i+1]
        cc = cc + 1
        axis.strata = c(axis.strata, x.left)
      }

      if(show.axis)
        axis(1, col = 'grey75', at = axis.strata, labels = FALSE, lwd = 2, line = 0.5)
      #axis(1, col = 'grey75', at = axis.strata, labels = seq(ba, 0, by=-s1) )

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
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
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
      # use attr(f, "ages"); confirm the age assignment below is also correct
      if(binned){
        y = sapply(f$h, function(x) max(which(horizons.min <= x)) )
        points(max(xx) - horizons.max[y] + (rev(s1)[y]/2), yy[fossils$sp] , col = fcol, pch = 19, cex = fcex)
      }
      else
        points(max(xx) - fossils$h, yy[fossils$sp], col = fcol, pch = 19, cex = fcex)
    }

    # ranges
    if(show.ranges){
      buffer = 0.01 * max(xx) # buffer for singletons
      s2 = 0
      if(binned)
        s2 = (ba / strata)/2

      for(i in unique(fossils$sp)) {

        range = fossils$h[which(fossils$sp==i)]
        if(length(range) == 1)
          range = c(range-buffer,range+buffer)

        species = rep(yy[i], length(range))

        lines(y = species, x = max(xx)-range+s2, lwd = 6, col=fcol)
      }
    }

    # water depth profile
    add.depth.axis = TRUE
    if(show.profile){
      add.depth.profile(depth.profile,axis.strata,strata,show.axis,add.depth.axis)
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

add.depth.profile<-function(depth.profile,axis.strata,strata,show.axis,add.depth.axis,show.preferred.depth=TRUE,PD=1){
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
    lines(x = axis.strata, y = rep(100, length(axis.strata)), col = "grey75", lwd = 2, lty = 3)
  if(show.axis){
    axis(1, col = 'grey75', at = axis.strata, labels = FALSE, lwd = 2)
    mtext(1, col = 'grey75', text="Time before present", line = 1)
  }
  if(add.depth.axis){
    axis(2, col = 'grey75', labels = TRUE, lwd = 2, las = 2, col.axis = 'grey75', line = 0.5, cex.axis = .7)
    mtext(2, col = 'grey75', text="Sampling proxy", line= 2)
  }
}


