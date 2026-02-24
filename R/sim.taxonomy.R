#' Simulate taxonomy
#'
#' Simulate a taxonomy object relating species identity to a phylo object under a mixed model of speciation.
#' Anagenetic and cryptic species can also be added later using the \code{sim.anagenetic.species} and \code{sim.cryptic.species} functions.
#'
#' @param tree Phylo object.
#' @param beta Probability of bifurcating speciation. Default = 0.
#' @param lambda.a Rate of anagenetic speciation. Default = 0.
#' @param kappa Probability that speciation event is cryptic. Default = 0.
#' @param root.edge If TRUE include root edge. Default = TRUE.
#' @return An object of class taxonomy.
#'
#' @examples
#' t = ape::rtree(10)
#' sim.taxonomy(t, 0.5, 0.1, 0.5)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
sim.taxonomy = function(tree, beta = 0, lambda.a = 0, kappa = 0, root.edge = TRUE){
  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!(beta >= 0 && beta <= 1))
    stop("beta must be a probability between 0 and 1")
  if(lambda.a < 0)
    stop("lambda.a must be zero or positive")
  if(!(kappa >= 0 && kappa <= 1))
    stop("kappa must be a probability between 0 and 1")
  
  if(is.null(tree$edge.length))
    stop("tree must have edge lengths")
  
  if(!ape::is.rooted(tree))
    stop("tree must be rooted")
  
  # assign symmetric and asymmetric species
  node.ages = n.ages(tree)
  
  species<-data.frame(sp = integer(),
                      edge = integer(),
                      parent = integer(),
                      start = numeric(),
                      end = numeric(),
                      mode = character(),
                      cryptic = logical(),
                      cryptic.id = integer()
  )
  
  # identify the root
  root = length(tree$tip.label) + 1
  
  if(root.edge && exists("root.edge",tree)){
    start = node.ages[root] + tree$root.edge
    mode = "o"
  } else {
    start = node.ages[root]
    mode = "r"
  }
  
  species <- rbind(species, data.frame(sp = root, edge = root, parent = 0, start = start, end = node.ages[root],
                                       mode = mode, cryptic = 0, cryptic.id = root))
  
  aux = function(node, p) {
    # fetch the two descendants
    descendants = tree$edge[which(tree$edge[,1] == node), 2]
    if(length(descendants) == 0) {
      return(p)
    }
    d1 <- descendants[1]
    d2 <- descendants[2]
    
    if(beta == 1 || (beta > 0 && runif(1) > (1 - beta))){
      # speciation event is symmetric
      a <- p$sp[which(p$edge == node)]
      
      p <- rbind(p, data.frame(sp = d1, edge = d1, parent = a, start = node.ages[a], end = node.ages[d1],
                               mode="s", cryptic = 0, cryptic.id = d1))
      p <- rbind(p, data.frame(sp = d2, edge = d2, parent = a, start = node.ages[a], end = node.ages[d2],
                               mode="s", cryptic = 0, cryptic.id = d2))
    } else{
      # speciation event is asymmetric/budding
      
      # update all labels & ages associated with d1 since (d1 is younger now)
      p$sp[which(p$sp == node)] = d1
      p$parent[which(p$parent == node)] = d1
      p$cryptic.id[which(p$cryptic.id == node)] = d1
      
      a <- p$parent[which(p$edge == node)]
      m <- p$mode[which(p$sp == d1)][1]
      
      p <- rbind(p, data.frame(sp = d1, edge = d1, parent = a, start = node.ages[ancestor(d1, tree)], end = node.ages[d1],
                               mode = m, cryptic = 0, cryptic.id = d1))
      # new species
      start = node.ages[ancestor(d2, tree)]
      p <- rbind(p, data.frame(sp = d2, edge = d2, parent = d1, start = start, end = node.ages[d2],
                               mode="b", cryptic = 0, cryptic.id = d2))
    }
    
    p = aux(d1, p)
    p = aux(d2, p)
    p
  }
  species = aux(root, species)
  
  species = species[order(species$sp),]
  
  species = taxonomy(species)
  
  # simulate anagenetic species
  if(lambda.a > 0) species = sim.anagenetic.species(tree, species, lambda.a)
  
  # simulate cryptic species
  if(kappa > 0) species = sim.cryptic.species(species, kappa)
  
  rownames(species) = NULL
  return(species)
}

#' Simulate anagenetic species on a taxonomy object
#'
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param lambda.a Rate of anagenetic speciation. Default = 0.
#' @return Object of class taxonomy.
#'
#' @examples
#' t = ape::rtree(10)
#' sp = sim.taxonomy(t, 1)
#' sim.anagenetic.species(t, sp, 0.1)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
sim.anagenetic.species = function(tree, species, lambda.a){
  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")
  if(any(species$mode=="a"))
    stop("taxonomy object already contains anagenetic species")
  if(lambda.a < 0)
    stop("lambda.a must be zero or positive")
  
  node.ages = n.ages(tree)
  
  if(lambda.a == 0) return(species)
  
  # identify the root
  root = root(tree)
  
  # identify the start of the new species id count
  # this is set to greater than the greatest edge label id to avoid any confusion with species labels
  # since species labels can change with the addition of anagenetic species
  species.counter = max(as.numeric(names(node.ages))) + 1
  
  # calculate edge start and end times
  edges = data.frame(edge = numeric(), start = numeric(), end = numeric())
  
  for(i in unique(species$edge)){
    edges = rbind(edges, data.frame(edge = i, start = species$start[which(species$edge == i)][1],
                                    end = species$end[which(species$edge == i)][1]))
  }
  
  for(i in unique(species$sp)){
    sp.start = max(species$start[which(species$sp == i)])
    sp.end = min(species$end[which(species$sp == i)])
    
    # sample n.eventsom number from a poisson distribution with rate = branch duration x lambda.a
    n.events = rpois(1, (sp.start - sp.end)  * lambda.a)
    
    if (n.events > 0) {
      
      # associated edges
      potential.edges = species$edge[which(species$sp == i)]
      p = subset(edges, edges$edge %in% potential.edges)
      p = p[order(p$end, decreasing = TRUE),]
      
      # subset edges associated with i from species df
      s = subset(species, species$edge %in% potential.edges)
      # remove edges associated with i from species df
      species = subset(species, !species$edge %in% potential.edges)
      
      # sample anagenetic speciation times from the branch duration
      h = runif(n.events, min = sp.end, max = sp.start)
      
      # order the new speciation times
      h = sort(h, decreasing = T)
      
      for(j in 1:(n.events+1)){
        
        # assign new identifier using species.counter - last label doesn't change
        sp = if(j == n.events + 1) i else species.counter
        start = if(j == 1) sp.start else h[j-1]
        end = if(j == n.events + 1) sp.end else h[j]
        
        # this is the edge on which the speciation event occurred
        edge = if(j == 1) p$edge[which(p$start == start)] else p$edge[which(p$start > start & p$end < start)]
        
        # additional younger edges associated with sp
        edge = c(edge, p$edge[which(p$start < start & p$start > end)])
        
        edge.start = p$start[which(p$edge %in% edge)]
        edge.end = p$end[which(p$edge %in% edge)]
        
        # replace first edge.start and last edge.end with species j start and end
        edge.start[which.max(edge.start)] = start
        edge.end[which.min(edge.end)] = end
        
        mode = if(j == 1) s$mode[1] else "a"
        parent = if(j == 1) s$parent[1] else species.counter - 1
        
        species <- rbind(species, data.frame(sp = sp, edge = edge, parent = parent, start = edge.start, end = edge.end,
                                             mode = mode, cryptic = 0, cryptic.id = sp))
        
        # deal with potential asym descendants - only if species label changed
        if(j != n.events + 1) {
          # find species with parent == i and beginning of species within range of new species
          for(desc_sp in unique(species$sp[species$parent == i])) {
            sp_start = max(species$start[species$sp == desc_sp])
            if(sp_start <= start && sp_start > end) species$parent[species$sp == desc_sp] = sp
          }
        }
        
        if(j != n.events + 1) species.counter = species.counter + 1
      }
    }
  }
  
  rownames(species) = NULL
  return(species)
}

#' Simulate cryptic species on a taxonomy object
#'
#' @param species Taxonomy object.
#' @param kappa Probability that speciation event is cryptic.
#' @return An object of class taxonomy. Note the origin or root can not be cryptic.
#'
#' @examples
#' t = ape::rtree(10)
#' sp = sim.taxonomy(t, 1)
#' sim.cryptic.species(sp, 0.5)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
sim.cryptic.species = function(species, kappa){
  if(!"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")
  if(any(species$cryptic == 1))
    stop("taxonomy object already contains cryptic species")
  if(!(kappa >= 0 && kappa <= 1))
    stop("kappa must be a probability between 0 and 1")
  
  if(kappa > 0){
    for(i in unique(species$sp)){
      
      # origin or root can not be cryptic
      if(species$mode[which(species$sp == i)][1] == "o") next
      if(species$mode[which(species$sp == i)][1] == "r") next
      
      if(runif(1) < kappa){
        
        # speciation event is cryptic
        species$cryptic[which(species$sp == i)] = 1
        
        # cryptic id based on crytptic parent id
        parent = species$parent[which(species$sp == i)][1]
        c = species$cryptic.id[which(species$sp == parent)][1]
        
        # assign cryptic species id using parent identity
        species$cryptic.id[which(species$sp == i)] = c
        
        # reassign cryptic ids that depend on the current species
        species$cryptic.id[which(species$cryptic.id == i)] = c
        
      }
    }
  }
  return(species)
}
