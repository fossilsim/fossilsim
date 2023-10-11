#' Plot oriented tree with stratigraphic ranges
#'
#' @param x object of type \code{fbdrange} containing orientation and range data
#' @param smart.labels whether to label the ranges (default \code{FALSE}, requires package [ggrepel] to place labels)
#' @param ... other arguments to be passed to the plot labels (does nothing if \code{smart.labels = FALSE})
#'
#' @return a ggtree plot which can be combined with any other commands from [ggplot2] or [ggtree]
#' @export
#'
#' @examples
#' @importFrom ggplot2 aes
#' @importFrom utils modifyList
plot.fbdrange <- function(x, smart.labels = FALSE, ...) {
  p1 <- ggplot2::ggplot(x, aes(color = range)) + ggtree::geom_tree(linewidth=2) + ggplot2::geom_point(aes(color = range), size = 1.3) 
  
  ## smartly spaced out tip/range labels
  if(smart.labels) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Package ggrepel is needed for smart labels. Please install it.", call. = FALSE)
    }
    
    labels <- x@phylo$tip.label
    species <- remove_last_substring(labels)
    species_count <- table(species)
    
    ## to avoid labeling the range twice, we make amty labels for first occurence
    new_labels <- sapply(1:length(labels), function(l) {
      if (species_count[species[l]] == 1 | grepl("last", labels[l])) gsub("_", " ", species[l])
      else ""
    })
    
    ## create a dataframe where each occurence label maps to newly created label
    labeldf <- data.frame(label = labels, new_labels = new_labels)
    
    default_args = list(size = 3, xlim = c(NA, Inf), min.segment.length = 0, force = 0.3, nudge_x = 0.5, direction = "both", hjust = 0, segment.size = 0.2)
    label_args = modifyList(default_args, list(...))
    p1 <- p1 + do.call(ggrepel::geom_label_repel, c(list(labeldf, aes(label = new_labels, alpha = 0.7, fontface = 4), label_args)))
  }
  
  p1
}

#' @export
#' @importFrom ggplot2 fortify
#' @importFrom rlang .data
fortify.fbdrange <- function (model, data, layout = "rectangular", ladderize = TRUE, 
                              right = FALSE, branch.length = "branch.length", mrsd = NULL, 
                              as.Date = FALSE, yscale = "none", root.position = 0, ...) {
  x <- ape::as.phylo(model)
  label <- x$tip.label[x$edge[,2][x$edge[,2]<=length(x$tip.label)]]
  if (ladderize) {
    x <- ape::ladderize(x, right = right)
  }
  
  if (!is.null(x$edge.length)) {
    if (anyNA(x$edge.length)) {
      warning("'edge.length' contains NA values...\n## setting 'edge.length' to NULL automatically when plotting the tree...")
      x$edge.length <- NULL
    }
  }
  
  if (layout %in% c("equal_angle", "daylight", "ape")) {
    res <- ggtree:::layout.unrooted(model, layout.method = layout, branch.length = branch.length, ...)
  }
  else {
    dd <- tidytree::arrange(tidytree::get.data(model), .data$node)
    
    ypos <- getYcoord.range(x, node.orientation=dd$orientation, match="ancestor",
                            tip.order=label)
    N <- ape::Nnode(x, internal.only = FALSE)
    if (is.null(x$edge.length) || branch.length == "none") {
      if (layout == "slanted") {
        sbp <- .convert_tips2ancestors_sbp(x, include.root = TRUE)
        xpos <- getXcoord_no_length_slanted(sbp)
        ypos <- getYcoord_no_length_slanted(sbp)
      }
      else {
        xpos <- ggtree:::getXcoord_no_length(x)
      }
    }
    else {
      xpos <- ggtree:::getXcoord(x)
    }
    xypos <- tidytree::tibble(node = 1:N, x = xpos + root.position, y = ypos)
    df <- tidytree::mutate(tidytree::as_tibble(model), isTip = !.data$node %in% .data$parent)
    res <- tidytree::full_join(df, xypos, by = "node")
  }
  res <- ggtree:::calculate_branch_mid(res, layout = layout)
  if (!is.null(mrsd)) {
    res <- ggtree:::scaleX_by_time_from_mrsd(res, mrsd, as.Date)
  }
  if (layout == "slanted") {
    res <- ggtree:::add_angle_slanted(res)
  }
  else {
    res <- ggtree:::calculate_angle(res)
  }
  res <- ggtree:::scaleY(ape::as.phylo(model), res, yscale, layout, ...)
  res <- adjust_hclust_tip.edge.len(res, x)
  class(res) <- c("tbl_tree", class(res))
  attr(res, "layout") <- layout
  return(res)
}

getYcoord.range <- function(tr, step=5, tip.order = NULL, node.orientation = NULL, match = NULL) {
  Ntip <- length(tr[["tip.label"]])
  N <- ape::Nnode(tr, internal.only = FALSE)
  
  edge <- tr[["edge"]]
  edge_length <- tr[["edge.length"]]
  parent <- edge[,1]
  child <- edge[,2]
  
  label<-tr$tip.label[edge[,2][edge[,2]<=Ntip]]
  
  cl <- split(child, parent)
  child_list <- list()
  child_list[as.numeric(names(cl))] <- cl
  
  y <- numeric(N)
  if (is.null(tip.order)) {
    tip.idx <- child[child <= Ntip]
    y[tip.idx] <- 1:Ntip * step
  } else {
    tip.idx <- 1:Ntip
    y[tip.idx] <- match(tr$tip.label, tip.order) * step
  }
  
  for (t in 1:Ntip){
    t_ <- tip.order[1]
    i <-1
    s = 0
    while(tr$tip.label[t]!=t_){
      tip_n = which(tr$tip.label==t_)
      row_id = which(child==tip_n)
      if (edge_length[row_id]==0){
        s=s+1
      }
      i <- i+1
      t_<-tip.order[i]
    }
    y[t] <- y[t]-step*s
  }
  
  y[-tip.idx] <- NA
  pvec <- edge2vec(tr)
  
  currentNode <- 1:Ntip
  while(anyNA(y)) {
    pNode <- unique(pvec[currentNode])
    
    idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
    newNode <- pNode[idx]
    
    if (!is.null(node.orientation)){
      for (n in newNode){
        ch <- child_list[[n]]
        if (node.orientation[ch[1]]==node.orientation[ch[2]]){
          # choose the one with longer branch length
          l1 <- edge_length[which(edge[,1]==n & edge[,2]==ch[1])]
          if (l1!=0){
            y[n] <- y[ch[1]]
          } else{
            y[n] <- y[ch[2]]
          }
        } else{
          l1 <- edge_length[which(edge[,1]==n & edge[,2]==ch[1])]
          l2 <- edge_length[which(edge[,1]==n & edge[,2]==ch[2])]
          if ((node.orientation[ch[1]]==match & l1!=0) | l2==0){
            y[n] <- y[ch[1]]
          } else {
            y[n] <- y[ch[2]]
          }
        }
      }
    } else {
      y[newNode] <- sapply(newNode, function(i) {
        mean(y[child_list[[i]]], na.rm=TRUE)
      })
    }
    
    for (t in 1:Ntip){
      row_id = which(child==t)
      if (edge_length[row_id]==0){
        p = parent[row_id]
        ch <- child_list[[p]]
        ch <- ch[which(ch!=t)]
        y[t] <- y[ch]
      }
    }
    
    currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
  }
  
  return(y)
}

remove_last_substring <- function(str) {
  gsub("_[^_]*$", "", str)
}

### All functions below copied from ggtree
### Credit: Guangchuang Yu

edge2vec <- function(tr) {
  parent <- tr$edge[,1]
  child <- tr$edge[,2]
  
  ## use lookup table
  pvec <- integer(max(tr$edge))
  pvec[child] <- parent
  return(pvec)
}

adjust_hclust_tip.edge.len <- function(df, phylo){
  if (inherits(phylo, 'treedata')){
    tip.edge.len <- attr(phylo@phylo, 'tip.edge.len')
  }else{
    tip.edge.len <- attr(phylo, 'tip.edge.len')
  }
  if (!is.null(tip.edge.len)){
    mx <- max(df$x, na.rm=TRUE)
    df$x <- df$x - mx
    df$branch <- df$branch - mx
    df[df$isTip, "x", drop=TRUE] <- tip.edge.len
    attr(df, 'revts.done') = TRUE
  }                       
  return(df)
}

getXcoord_no_length_slanted <- function(x){
  x <- -colSums(x)
  x <- unname(x[order(as.numeric(names(x)))])
  x <- x + max(abs(x))
  return(x)
}

getYcoord_no_length_slanted <- function(y){
  y <- seq_len(nrow(y)) * y
  y[y==0] <- NA
  y <- colMeans(y, na.rm = TRUE)
  y <- unname(y[order(as.numeric(names(y)))])
  return(y)
}

.convert_tips2ancestors_sbp <- function (tree, include.root = FALSE){
  all.nodes <- unique(as.vector(tree@phylo$edge))
  if (!include.root) {
    all.nodes <- setdiff(all.nodes, tidytree::rootnode(tree))
  }
  
  tip.nodes <- .nodeId(tree, type = "tips")
  tn <- lapply(tip.nodes, 
               .internal_ancestor, 
               .data = tree, 
               all.nodes = all.nodes)
  sbp <- stats::setNames(tn, tip.nodes) 
  sbp <- do.call(rbind, sbp) 
  colnames(sbp) <- all.nodes
  return(sbp)
}

.nodeId <- function (tree, type = "all"){
  type <- match.arg(type, c("all", "tips", "internal"))
  if (inherits(tree, "treedata")) {
    tree <- tree@phylo
  }
  nodes <- unique(as.vector(tree$edge))
  if (type == "all") {
    return(nodes)
  }
  edge <- tree$edge
  tips <- edge[!edge[, 2] %in% edge[, 1], 2]
  if (type == "tips"){
    return(tips)
  }
  else if (type == "internal") {
    return(setdiff(nodes, tips))
  }
}

.internal_ancestor <- function(.data, .node, all.nodes){
  x <- tidytree::ancestor(.data=.data, .node=.node)
  x <- c(x, .node)
  x <- all.nodes %in% x
  return (x)
}