## adapted from ape::drop.tip to account for the taxonomy updates
## returns a list of tree and taxonomy
drop.tip.with.taxonomy = function(phy, remove_tips, taxonomy = NULL) {
  Ntip <- length(phy$tip.label)
  if(is.character(remove_tips)) tip <- which(phy$tip.label %in% remove_tips)
  else tip <- remove_tips
  
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) return(list(tree = phy, taxonomy = taxonomy))
  if (length(tip) == Ntip) {
    warning("drop all tips of the tree: returning NULL")
    return(NULL)
  }
  
  if (length(tip) == Ntip - 1) {
    i <- which(phy$edge[, 2] == (1:Ntip)[-tip])
    res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i, 2]], Nnode = 1L)
    class(res) <- "phylo"
    res$edge.length <- phy$edge.length[i]
    if (!is.null(phy$node.label)) res$node.label <- phy$node.label[phy$edge[i, 1] - Ntip]
    
    if(!is.null(taxonomy)) {
      taxonomy = taxonomy[taxonomy$edge == phy$edge[i, 2],]
      taxonomy$parent[which.max(taxonomy$start)] = 0
      taxonomy = rbind(taxonomy, list(sp = taxonomy$sp[which.max(taxonomy$start)], edge = 2, mode = "r", parent = 0, 
                                      start = max(taxonomy$start), end = max(taxonomy$start),
                                      cryptic = F, cryptic.id = taxonomy$sp[which.max(taxonomy$start)]))
    }
  }
  
  phy <- ape::reorder.phylo(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  
  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) break
    keep[sel] <- FALSE
  }
  
  if(!is.null(taxonomy)) taxonomy = taxonomy[taxonomy$edge %in% phy$edge[keep,],]
  
  phy$edge <- phy$edge[keep, ]
  phy$edge.length <- phy$edge.length[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  if (length(tip)) phy$tip.label <- phy$tip.label[-tip]
  
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label)) phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
  
  newNb[oldNo.ofNewTips] = rank(phy$edge[TERMS, 2])
  if(!is.null(taxonomy)) taxonomy$edge = newNb[taxonomy$edge]
  
  collapse.singles.taxonomy(phy, taxonomy)
}

## adapted from ape::collapse.singles to account for the taxonomy updates
## returns a list of tree and taxonomy
collapse.singles.taxonomy = function (tree, taxonomy = NULL) {
  n <- length(tree$tip.label)
  
  tree <- ape::reorder.phylo(tree)
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  tab <- tabulate(e1)
  if (all(tab[-c(1:n)] > 1)) return(list(tree = tree, taxonomy = taxonomy))
  
  el <- tree$edge.length
  ROOTEDGE <- if(!is.null(tree$root.edge)) tree$root.edge else 0
  ROOT <- n + 1L
  while (tab[ROOT] == 1) {
    i <- which(e1 == ROOT)
    NEWROOT <- e2[i]
    ROOTEDGE <- ROOTEDGE + el[i]
    el <- el[-i]
    e1 <- e1[-i]
    e2 <- e2[-i]
    
    if(!is.null(taxonomy)) {
      tax_id = find_last_id_on_edge(taxonomy, ROOT)
      new_tax_id = find_first_id_on_edge(taxonomy, NEWROOT)
      if(taxonomy$sp[tax_id] == taxonomy$sp[new_tax_id]) {
        taxonomy$end[tax_id] = taxonomy$end[new_tax_id]
        taxonomy = taxonomy[-new_tax_id,]
      } else {
        taxonomy$mode[new_tax_id] = "a"
      }
      taxonomy$edge[taxonomy$edge == ROOT] = NEWROOT
    }
    
    ROOT = NEWROOT
  }
  
  singles <- which(tabulate(e1) == 1)
  if (length(singles) > 0) {
    ii <- sort(match(singles, e1), decreasing = TRUE)
    jj <- match(e1[ii], e2)
    for (i in 1:length(singles)) {
      
      if(!is.null(taxonomy)) {
        tax_id = find_last_id_on_edge(taxonomy, e2[jj[i]])
        new_tax_id = find_first_id_on_edge(taxonomy, e2[ii[i]])
        if(taxonomy$sp[tax_id] == taxonomy$sp[new_tax_id]) {
          taxonomy$end[tax_id] = taxonomy$end[new_tax_id]
          taxonomy = taxonomy[-new_tax_id,]
        } else {
          if(taxonomy$end[tax_id] != taxonomy$start[new_tax_id]) stop(paste(e2[jj[i]], e2[ii[i]]))
          taxonomy$mode[new_tax_id] = "a"
        }
        taxonomy$edge[taxonomy$edge == e2[jj[i]]] = e2[ii[i]]
      }
      
      e2[jj[i]] <- e2[ii[i]]
      el[jj[i]] <- el[jj[i]] + el[ii[i]]
    }
    e1 <- e1[-ii]
    e2 <- e2[-ii]
    el <- el[-ii]
  }
  
  Nnode <- length(e1) - n + 1L
  oldnodes <- unique(e1)
  if (!is.null(tree$node.label)) tree$node.label <- tree$node.label[oldnodes - n]
  
  newNb <- integer(max(oldnodes))
  newNb[ROOT] <- n + 1L
  sndcol <- e2 > n
  e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
  e1 <- newNb[e1]
  
  if(!is.null(taxonomy)) {
    upd_tax = taxonomy$edge > n
    taxonomy$edge[upd_tax] = newNb[taxonomy$edge[upd_tax]]
  }
  
  tree$edge <- cbind(e1, e2, deparse.level = 0)
  tree$Nnode <- Nnode
  tree$root.edge <- ROOTEDGE
  tree$edge.length <- el
  
  list(tree = tree, taxonomy = taxonomy)
}

## find what is the last (i.e. youngest) species on a specific edge (identified by the bottom node)
## returns the corresponding row idx in the taxonomy object
find_last_id_on_edge = function(taxonomy, edge) {
  idxs = which(taxonomy$edge == edge)
  last = which.min(taxonomy$end[idxs])
  idxs[last]
}

## find what is the first (i.e. oldest) species on a specific edge (identified by the bottom node)
## returns the corresponding row idx in the taxonomy object
find_first_id_on_edge = function(taxonomy, edge) {
  idxs = which(taxonomy$edge == edge)
  first = which.max(taxonomy$start[idxs])
  idxs[first]
}