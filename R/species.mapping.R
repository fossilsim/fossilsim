#' Create taxonomy
#'
#' @param tree Phylo object.
#' @param beta Probability of bifurcating speciation. Default = 0.
#' @param lambda.a Rate of anagenic speciation. Default = 0.
#' @param kappa Probability that speciation event is cryptic. Default = 0.
#' @param root.edge If TRUE include root edge. Default = TRUE.
#' @return An object of class taxonomy.
#'
#' @examples
#' t = ape::rtree(10)
#' create.taxonomy(t, 0.5, 0.1, 0.5)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
create.taxonomy<-function(tree, beta = 0, lambda.a = 0, kappa = 0, root.edge = TRUE){
  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!(beta >= 0 && beta <= 1))
    stop("beta must be a probability between 0 and 1")
  if(lambda.a < 0)
    stop("lambda.a must be zero or positive")
  if(!(kappa >= 0 && kappa <= 1))
    stop("kappa must be a probability between 0 and 1")

  # assign symmetric and asymmetric species
  node.ages = n.ages(tree)

  species<-data.frame(sp = integer(),
                      edge = integer(),
                      parent = integer(),
                      start = numeric(),
                      end = numeric(),
                      mode = character(),
                      origin = integer(),
                      cryptic = integer(),
                      cryptic.id = integer(),
                      edge.start = numeric(),
                      edge.end = numeric()
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
                                       mode = mode, origin = root, cryptic = 0, cryptic.id = root,
                                       edge.start = start, edge.end = node.ages[root]))

  aux = function(node, p) {
    # fetch the two descendants
    descendants = tree$edge[which(tree$edge[,1] == node), 2]
    if(length(descendants) == 0) {
      return(p)
    }
    d1 <- descendants[1]
    d2 <- descendants[2]

    if(runif(1) > (1 - beta)){
      # speciation event is symmetric
      a <- p$sp[which(p$edge == node)]
      p <- rbind(p, data.frame(sp = d1, edge = d1, parent = a, start = node.ages[a], end = node.ages[d1],
                               mode="s", origin = d1, cryptic = 0, cryptic.id = d1,
                               edge.start = node.ages[a], edge.end = node.ages[d1]))
      p <- rbind(p, data.frame(sp = d2, edge = d2, parent = a, start = node.ages[a], end = node.ages[d2],
                               mode="s", origin = d2, cryptic = 0, cryptic.id = d2,
                               edge.start = node.ages[a], edge.end = node.ages[d2]))
    } else{
      # speciation event is asymmetric/budding

      # update all labels & ages associated with d1 since (d1 is younger now)
      p$sp[which(p$sp == node)] = d1
      p$end[which(p$sp == d1)] = node.ages[d1]
      p$parent[which(p$parent == node)] = d1
      p$cryptic.id[which(p$cryptic.id == node)] = d1

      a <- p$parent[which(p$edge == node)]
      s <- p$start[which(p$sp == d1)][1]
      m <- p$mode[which(p$sp == d1)][1]
      o <- p$origin[which(p$sp == d1)][1]
      p <- rbind(p, data.frame(sp = d1, edge = d1, parent = a, start = s, end = node.ages[d1],
                               mode = m, origin = o, cryptic = 0, cryptic.id = d1,
                               edge.start = node.ages[ancestor(d1, tree)], edge.end = node.ages[d1]))
      # new species
      start = node.ages[ancestor(d2, tree)]
      p <- rbind(p, data.frame(sp = d2, edge = d2, parent = d1, start = start, end = node.ages[d2],
                               mode="b", origin = d2, cryptic = 0, cryptic.id = d2,
                               edge.start = start, edge.end = node.ages[d2]))
    }

    p = aux(d1, p)
    p = aux(d2, p)
    p
  }
  species = aux(root, species)

  species = species[order(species$sp),]

  species = taxonomy(species)

  # simulate anagenic species
  if(lambda.a > 0) species = add.anagenic.species(tree, species, lambda.a)

  # simulate cryptic species
  if(kappa > 0) species = add.cryptic.species(species, kappa)

  rownames(species) = NULL
  return(species)
  # eof
}

#' Add anagenic species to a taxonomy object
#'
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param lambda.a Rate of anagenic speciation. Default = 0.
#' @return Object of class taxonomy.
#'
#' @examples
#' t = ape::rtree(10)
#' sp = create.taxonomy(t, 1)
#' add.anagenic.species(t, sp, 0.1)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
add.anagenic.species<-function(tree, species, lambda.a){
  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")
  if(any(species$mode=="a"))
    warning("taxonomy object already contains anagenic species")
  if(lambda.a < 0)
    stop("lambda.a must be zero or positive")

  node.ages = n.ages(tree)

  if(lambda.a > 0){

    # identify the root
    root = root(tree)

    # identify the start of the new species id count
    # this is set to greater than the greatest edge label id to avoid any confusion with species labels
    # since species labels can change with the addition of anagenic species
    species.counter = max(as.numeric(names(node.ages))) + 1

    # calculate edge start and end times (I think this is probably not necessary now)
    edges = data.frame(edge = numeric(), start = numeric(), end = numeric())

    for(i in unique(tree$edge[,2])){
      start = node.ages[ancestor(i, tree)]
      end = node.ages[i]
      edges = rbind(edges, data.frame(edge = i, start = start, end = end))
    }
    # deal with the root / root edge
    root.age = node.ages[root]
    if(exists("root.edge",tree)){
      edges = rbind(edges, data.frame(edge = root, start = root.age + tree$root.edge, end = root.age))
    } else edges = rbind(edges, data.frame(edge = root, start = root.age, end = root.age))

    # order by edge label so that root edge label appears before descendant labels
    edges = edges[order(edges$edge),]

    for(i in unique(species$sp)){

      sp.start = species$start[which(species$sp == i)][1]
      sp.end = species$end[which(species$sp == i)][1]

      # sample random number from a poisson distribution with rate = branch duriation x lambda.a
      rand = rpois(1, (sp.start - sp.end)  * lambda.a)

      if (rand > 0) {

        # associated edges
        potential.edges = species$edge[which(species$sp == i)]
        p = subset(edges, edges$edge %in% potential.edges)

        # subset edges associated with i from species df
        s = subset(species, species$edge %in% potential.edges)
        # remove edges associated with i from species df
        species = subset(species, !species$edge %in% potential.edges)

        # sample anagenic speciation times from the branch duration
        h = runif(rand, min = sp.end, max = sp.start)

        # order the new speciation times
        h = sort(h, decreasing = T)

        for(j in 1:(length(h)+1)){

          # original species beginning the branch
          if(j == 1){

            # assign new identifier using species.counter
            sp = species.counter
            start = sp.start
            end = h[j]

            # this is the edge on which the speciation event occurred
            edge = p$edge[which(p$start == start)]
            # additional younger edges associated with sp
            edge = c(edge, p$edge[which(p$start < start & p$start > end)])

            edge.start = edges$start[which(edges$edge %in% edge)]
            edge.end = edges$end[which(edges$edge %in% edge)]

            #mode = as.character(s$mode[which(s$mode != "NA")])
            mode = s$mode[1]

            parent = s$parent[1]

            species <- rbind(species, data.frame(sp = sp, edge = edge, parent = parent, start = start, end = end,
                                                 mode = mode, origin = edge[1], cryptic = 0, cryptic.id = sp,
                                                 edge.start = edge.start, edge.end = edge.end))

            # deal with potential asym descendants
            if(any(species$parent == i & species$start <= start & species$start > end)) {

              species$parent[which(species$parent == i & species$start <= start
                                   & species$start > end)] = sp

            }
            species.counter = species.counter + 1
          } else if (j == (length(h)+1)){ # species ending the branch

            sp = i
            start = h[j-1]
            end = sp.end

            # this is the edge on which the speciation event occurred
            edge = p$edge[which(p$start > start & p$end < start)]
            # additional younger edges associated with sp
            edge = c(edge, p$edge[which(p$start < start & p$start > end)])

            edge.start = edges$start[which(edges$edge %in% edge)]
            edge.end = edges$end[which(edges$edge %in% edge)]

            mode = "a"

            parent = species.counter - 1

            species <- rbind(species, data.frame(sp = sp, edge = edge, parent = parent, start = start, end = end,
                                                 mode = mode, origin = edge[1], cryptic = 0, cryptic.id = sp,
                                                 edge.start = edge.start, edge.end = edge.end))

            # any descendant parent labels associated with this species shouldn't need to change

          } else { # intermediate anagenic species

            sp = species.counter
            start = h[j-1]
            end = h[j]

            # this is the edge on which the speciation event occurred
            edge = p$edge[which(p$start > start & p$end < start)]
            # additional younger edges associated with sp
            edge = c(edge, p$edge[which(p$start < start & p$start > end)])

            edge.start = edges$start[which(edges$edge %in% edge)]
            edge.end = edges$end[which(edges$edge %in% edge)]

            mode = "a"

            parent = species.counter - 1

            species <- rbind(species, data.frame(sp = sp, edge = edge, parent = parent, start = start, end = end,
                                                 mode = mode, origin = edge[1], cryptic = 0, cryptic.id = sp,
                                                 edge.start = edge.start, edge.end = edge.end))

            # deal with potential asym descendants
            if(any(species$parent == i & species$start < start & species$start > end)) {

              species$parent[which(species$parent == i & species$start < start
                                   & species$start > end)] = sp

            }
            species.counter = species.counter + 1
          }
        }

      }
    }
  }
  rownames(species) = NULL
  return(species)
}

#' Add cryptic species to a taxonomy object
#'
#' @param species Taxonomy object.
#' @param kappa Probability that speciation event is cryptic.
#' @return An object of class taxonomy. Note the origin or root can not be cryptic.
#'
#' @examples
#' t = ape::rtree(10)
#' sp = create.taxonomy(t, 1)
#' add.cryptic.species(sp, 0.5)
#'
#' @seealso \code{\link{taxonomy}}
#'
#' @export
add.cryptic.species<-function(species, kappa){
  if(!"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")
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
