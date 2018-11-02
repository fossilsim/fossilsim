#' Create a set of BEAST2 constraints to construct a DPPDIV style fixed extant
#' topology FBD analysis
#'
#' If complete = FALSE, only the extant taxa are used to construct the taxon
#' constraints, resulting in a DPPDIV style analysis in which the extant
#' topology is fixed and fossils can float in the tree. The resulting output uses the stronglyMonophyletic
#' taxon constraint on the root, this means that all fossil taxa will be sampled in the crown group, and never
#' in a position below the root.
#'
#' @param tree an object of class "phylo", representing the tree upon which the
#'   fossil occurrences were simulated.
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences for the "tree" argument.
#' @param file the name of the file to which the constraints will be written,
#'   defaults to "BEASTconstraints.xml".
#' @param complete logical, if TRUE then taxon constraints are built for the
#'   complete tree, if FALSE then constraints are built for the crown clades
#'   only. Default value is FALSE.
#' @param tree.name the name of the tree as used in the BEAST2 xml format.
#' @return NULL.
#' @export
fossils.to.BEAST.constraints <- function(fossils, tree, file = "BEASTconstraints.xml", complete = FALSE, tree.name = "beastTree"){

  if (missing(tree)) {
    stop("tree must be supplied")
  }

  if (missing(fossils)) {
    stop("fossils must be supplied")
  }

  if (any(fossils$sp == root(tree))) {
    stop("Can't handle fossil samples on the root.edge")
  }

  # if the complete tree is required
  if (complete) {
    anc <-
      place.fossils(tree = tree,
                    fossils = fossils)
  } else {
    # if the extant only tree is required
    ext.tree <- prune.fossil.tips(tree)
    anc <-
      place.fossils(tree = tree,
                    fossils = fossils,
                    ext.tree = ext.tree)
    tree <- ext.tree
  }

  if (!all(anc %in% tree$edge[, 1])) {
    stop("fossil mismatch, fossils belong to non existent nodes in the input tree")
  }

  fossils$anc_ext <- anc # assign the ancestral nodes and move on

  # Create some arbitrary names for the fossils
  fossils$names <- paste0("fossil_", c(1:length(fossils$sp)))

  # Writing Functions, these are based on the write.nexus.data function in ape
  filename <- file
  zz <- file(filename, "w")
  "fcat" <-
    function (..., file = zz)
    {
      cat(...,
          file = filename,
          sep = "",
          append = TRUE)
    }

  "ccat" <-
    function (..., file = zz)
    {
      cat(...,
          file = filename,
          sep = " ",
          append = TRUE)
    }

  # First, construct L/R calibs for the entire tree and store them in a list
  L <- list()
  R <- list()
  dec_store <- c()
  mult_dec <-c()
  nodes <- sort(unique(tree$edge[, 1]))
  j = 0
  for (node in nodes) {
    j <- j + 1
    # Get the direct descendent nodes
    #dec <- phangorn::Descendants(tree, node, "children")
    dec <- get.dec.nodes(tree, node)
    dec_store <- rbind(dec_store, dec)
    #dec_t <- phangorn::Descendants(tree, dec, "tips")
    dec_t <- list()
    dec_t[[1]] <-
      which(tree$tip.label %in% fetch.descendants(dec[1], tree = tree))
    dec_t[[2]] <-
      which(tree$tip.label %in% fetch.descendants(dec[2], tree = tree))

    # Do the L and R sides have multiple taxa?
    mult_dec <- rbind(mult_dec, sapply(dec_t, length) > 1)

    # Construct the L/R constraint
    R[[j]] <- tree$tip.label[dec_t[[1]]]
    L[[j]] <- tree$tip.label[dec_t[[2]]]
  }
  names(L) <- nodes
  names(R) <- nodes
  row.names(mult_dec) <- nodes
  colnames(mult_dec) <- c("Right", "Left")
  row.names(dec_store) <- nodes

  for (j in 1:length(fossils$anc_ext)) {
    #ancs <-
    #  c(phangorn::Ancestors(tree, fossils$anc_ext[j]),
    #    fossils$anc_ext[j]) # get the node number of all clades the fossil is in
    ancs <-
      c(find.edges.inbetween(
        j = min(tree$edge[, 1]),
        i = fossils$anc_ext[j],
        tree = tree
      )[-1],
      fossils$anc_ext[j]) # get the node number of all clades the fossil is in

    if (any(ancs == root(tree))) {
      ancs <- ancs[-which(ancs == root(tree))]
    }

    if (length(ancs) > 0) {
      # if we still have ancestors after removing the root
      for (k in 1:length(ancs)) {
        # add the fossil on the correct side
        # side_tmp <- phangorn::Ancestors(tree, ancs[k], "parent")
        side_tmp <-
          find.edges.inbetween(j = min(tree$edge[, 1]),
                               i = ancs[k],
                               tree = tree)[-1][1]
        side_tmp <- which(row.names(dec_store) == side_tmp)
        side <- which(dec_store[side_tmp, ] == ancs[k])

        base_clade <-
          row.names(dec_store)[side_tmp]

        if (side == 1) {
          # if the side is the right
          id <- which(names(R) == base_clade)
          R[[id]] <- c(fossils$names[j], R[[id]])
        }
        if (side == 2) {
          # if the side is the left
          id <- which(names(L) == base_clade)
          L[[id]] <- c(fossils$names[j], L[[id]])
        }
      }
    }
  }

  root_monophyly <-
    c(tree$tip.label, fossils$names)

  # Whole group Monophyly
  fcat(
    "            <distribution id=\"whole_group.prior\" stronglyMonophyletic=\"true\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:",
    tree.name,
    "\">\n"
  )
  fcat("                <taxonsetIn id=\"whole_group\" spec=\"TaxonSet\">\n")
  for (i in 1:length(root_monophyly)) {
    fcat("                    <taxon id=\"")
    fcat(root_monophyly[i])
    fcat("\" spec=\"Taxon\"/>\n")
  }
  fcat("                </taxonsetIn>\n")
  fcat("            </distribution>\n")


  # individual clades
  # For each clade, do we need both L AND R
  for (i in 1:length(L)) {
    lengths <- c()
    lengths[1] <- length(R[[i]])
    lengths[2] <- length(L[[i]])
    # dec_node_id <- phangorn::Descendants(tree, nodes[i], "children")
    dec_node_id <- get.dec.nodes(tree, nodes[i])
    if (all(lengths > 1)) {
      # if both sides of the bifurcation have multiple taxa
      # L in, R out: R in ,L out
      fcat("            <distribution id=\"Const_")
      fcat(
        dec_node_id[2],
        ".prior\" stronglyMonophyletic=\"true\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:",
        tree.name,
        "\">\n"
      )
      # Define the In constraint
      fcat("                <taxonsetIn id=\"Const_")
      fcat(dec_node_id[2], "\" spec=\"TaxonSet\">\n")
      for (j in 1:length(L[[i]])) {
        fcat("                    <taxon idref=\"")
        fcat(L[[i]][j], "\" spec=\"Taxon\"/>\n")
      }
      fcat("                </taxonsetIn>\n")
      # Define the Out constraint
      fcat("                <taxonsetOut ")
      fcat(" spec=\"TaxonSet\">\n")
      for (j in 1:length(R[[i]])) {
        fcat("                    <taxon idref=\"")
        fcat(R[[i]][j], "\" spec=\"Taxon\"/>\n")
      }
      fcat("                </taxonsetOut>\n")
      fcat("            </distribution>\n")

      # Make the constraint symmetrical
      fcat("            <distribution id=\"Const_")
      fcat(
        dec_node_id[1],
        ".prior\" stronglyMonophyletic=\"true\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:",
        tree.name,
        "\">\n"
      )
      # Define the In constraint
      fcat("                <taxonsetIn id=\"Const_")
      fcat(dec_node_id[1], "\" spec=\"TaxonSet\">\n")
      for (j in 1:length(R[[i]])) {
        fcat("                    <taxon idref=\"")
        fcat(R[[i]][j], "\" spec=\"Taxon\"/>\n")
      }
      fcat("                </taxonsetIn>\n")
      # Define the Out constraint
      fcat("                <taxonsetOut ")
      fcat(" spec=\"TaxonSet\">\n")
      for (j in 1:length(L[[i]])) {
        fcat("                    <taxon idref=\"")
        fcat(L[[i]][j], "\" spec=\"Taxon\"/>\n")
      }
      fcat("                </taxonsetOut>\n")
      fcat("            </distribution>\n")
    }

    # if only one side has multiple taxa, make this the taxonsetIn side
    if (sum(lengths == 1) == 1) {
      in_side <- which(lengths > 1)

      if (in_side == 1) {
        #if the right side has multiple taxa and the left has one taxon
        fcat("            <distribution id=\"Const_")
        fcat(
          dec_node_id[1],
          ".prior\" stronglyMonophyletic=\"true\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:",
          tree.name,
          "\">\n"
        )
        # Define the In constraint
        fcat("                <taxonsetIn id=\"Const_")
        fcat(dec_node_id[1], "\" spec=\"TaxonSet\">\n")
        for (j in 1:length(R[[i]])) {
          fcat("                    <taxon idref=\"")
          fcat(R[[i]][j], "\" spec=\"Taxon\"/>\n")
        }
        fcat("                </taxonsetIn>\n")
        # Define the Out constraint
        fcat("                <taxonsetOut ")
        fcat(" spec=\"TaxonSet\">\n")
        for (j in 1:length(L[[i]])) {
          fcat("                    <taxon idref=\"")
          fcat(L[[i]][j], "\" spec=\"Taxon\"/>\n")
        }
        fcat("                </taxonsetOut>\n")
        fcat("            </distribution>\n")
      }
      if (in_side == 2) {
        # if the left side has multiple taxa and the right has a single taxon
        fcat("            <distribution id=\"Const_")
        fcat(
          dec_node_id[2],
          ".prior\" stronglyMonophyletic=\"true\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:",
          tree.name,
          "\">\n"
        )
        # Define the In constraint
        fcat("                <taxonsetIn id=\"Const_")
        fcat(dec_node_id[2], "\" spec=\"TaxonSet\">\n")
        for (j in 1:length(L[[i]])) {
          fcat("                    <taxon idref=\"")
          fcat(L[[i]][j], "\" spec=\"Taxon\"/>\n")
        }
        fcat("                </taxonsetIn>\n")
        # Define the Out constraint
        fcat("                <taxonsetOut ")
        fcat(" spec=\"TaxonSet\">\n")
        for (j in 1:length(R[[i]])) {
          fcat("                    <taxon idref=\"")
          fcat(R[[i]][j], "\" spec=\"Taxon\"/>\n")
        }
        fcat("                </taxonsetOut>\n")
        fcat("            </distribution>\n")
      }
    }
  }
  close(zz)

}

#' Create a suitable starting tree for a DPPDIV style FBD analysis in BEAST2
#'
#' If complete = FALSE, only the extant taxa are used to construct the tree,
#' resulting in a DPPDIV style analysis in which the extant topology is fixed
#' and fossils can float in the tree.
#'
#' @param tree an object of class "phylo", representing the tree upon which the
#'   fossil occurrences were simulated.
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences for the "tree" argument.
#' @param complete logical, if TRUE then the tree are built for the complete
#'   tree, if FALSE then the tree is built for the crown clades only.
#' @return a string representing the starting tree in newick format.
#' @export
fossils.to.BEAST.start.tree <- function(tree, fossils, complete = FALSE){

  if (any(fossils$sp == min(tree$edge[,1]))) {
    stop("Can't handle fossil samples on the root.edge")
  }

  if("root.edge" %in% names(tree)){ # if a root edge exists, remove its length
    tree$root.edge <- NULL
  }

  # test for the branch lengths
  brlen <- tree$edge.length
  if(length(brlen) == 0){
    stop("Tree must have branch lengths")
  }

  # if the output from the complete tree is required
  if (complete) {
    anc <- place.fossils(tree = tree, fossils = fossils)
    aug_tree <- tree
  } else {
    # if the extant only tree is required
    ext.tree <- prune.fossil.tips(tree)
    fossils <- remove.stem.fossils(fossils, tree)
    anc <-
      place.fossils(tree = tree,
                    fossils = fossils,
                    ext.tree = ext.tree)
    aug_tree <- ext.tree
    tree <- ext.tree
  }
  fossils$anc_ext <- anc # assign the ancestral nodes and move on

  for (i in 1:length(fossils$anc_ext)) {
    if (fossils$anc_ext[i] != 0) {
      # if the subtending node is not the origin
      # Get the descendent tips in the extant only tree
      # dec <-
      #  phangorn::Descendants(node = fossils$anc_ext[i], x = tree)[[1]]
      dec <-
        which(tree$tip.label %in% fetch.descendants(edge = fossils$anc_ext[i],tree = tree))
      dec <- tree$tip.label[dec] # get the names of descendent taxa
      dec <-
        ape::getMRCA(aug_tree, dec)
      # dec <-
      #  phangorn::Descendants(node = dec, x = aug_tree)[[1]]
      dec <-
        which(aug_tree$tip.label %in% fetch.descendants(edge = dec, tree = aug_tree))
      dec <-
        sample(x = dec, size = 1) #randomly choose a sister taxon in the descendents
      # Suppress warnings as bind.tip complains about binding a single tip
      #aug_tree <-
      #  suppressWarnings(phytools::bind.tip(
      #    aug_tree,
      #    tip.label = paste0("fossil_", i),
      #    where = dec
      #  ))
      aug_tree <- bind.to.tip(aug_tree, dec, paste0("fossil_", i))

      aug_tree$edge.length <-
        rep(1, length(aug_tree$edge.length)) # remove brLens
    } else {
      # if it is the origin, the fossil needs to be added outside the crown
      # this is not possible
      stop("fossil is not a member of the crown group")

      # crown <- prune.fossil.tips(tree) # get taxa that are outside of the crown
      # crown <- crown$tip.label
      # crownNode <- ape::getMRCA(tree, crown)
      # crownTax <- fetch.descendants(crownNode, tree)
      # stemTax <- setdiff(fetch.descendants(min(tree$edge[,1]), tree), crownTax)
      # placement <- sample(x = stemTax, size = 1)
      # aug_tree <- bind.to.tip(aug_tree, placement, paste0("fossil_", i))
      #
      # aug_tree$edge.length <-
      #   rep(1, length(aug_tree$edge.length)) # remove brLens
    }
  }
  aug_tree$edge.length <-
    rep(NA, length(aug_tree$edge.length)) # remove brLens
  aug_newick <-
    gsub(pattern = ":NA|;",
         ape::write.tree(aug_tree),
         replacement = "") # Convert to newick

  return(aug_newick)
}
