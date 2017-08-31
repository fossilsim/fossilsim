#' Add extant and extinct tip samples
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param rho Tip species sampling.
#' @return An object of class fossils.
#' sp = edge labels. h = ages.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # add tip samples
#' f<-add.tip.samples(t, f, rho = 0.5)
#' plot(f, t)
#' @export
add.tip.samples<-function(tree, fossils, rho = 1) {

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")
  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")
  node.ages <- n.ages(tree)
  for (i in tree$edge[, 2]) {
    if(is.tip(i,tree)){
      row = which(tree$edge[, 2] == i)
      ancestor = tree$edge[, 1][row]
      a = which(names(node.ages) == ancestor)
      lineage.start = node.ages[[a]]
      b = tree$edge.length[row]
      lineage.end = lineage.start - b
      lineage.end = round(lineage.end, 7)
      if (runif(1) < rho)
        fossils <- rbind(fossils, data.frame(h = lineage.end, sp = i))
    }
  }
  return(fossils)
}

#' Add extant occurrence samples
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param rho Extant species sampling.
#' @return An object of class fossils.
#' sp = edge labels. h = ages.
#' @examples
#' # simulate tree
#' lambda = 0.1
#' mu = 0.05
#' tips = 8
#' t<-TreeSim::sim.bd.taxa(tips, 1, lambda, mu)
#' t<-t[[1]]
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 0.5, root.edge = FALSE)
#' # add extant samples
#' f<-add.extant.occ(t, f, rho = 0.5)
#' plot(f, t)
#' @export
add.extant.occ<-function(tree, fossils, rho = 1){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")
  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")

  # store speciation mode
  speciation = attr(fossils, "speciation")

  # work out node ages
  node.ages<-n.ages(tree)

  for (i in tree$edge[,2]){ # internal nodes + tips

    # work out the max age of the lineage (e.g. when that lineage became extant)
    # & get ancestor
    row=which(tree$edge[,2]==i)
    ancestor=tree$edge[,1][row]

    # get the age of the ancestor
    a=which(names(node.ages)==ancestor)
    lineage.start=node.ages[[a]]

    # work out the min age of the lineage (e.g. when that lineage became extinct)
    # & get the branch length
    b=tree$edge.length[row]
    lineage.end=lineage.start-b # branch length

    lineage.end=round(lineage.end,7)

    if(lineage.end==0){
      if(runif(1) < rho)
        fossils<-rbind(fossils,data.frame(h=0,sp=i))
    }
  }
  if(!is.fossils(fossils))
    fossils = as.fossils(fossils, speciation.mode = speciation)
  return(fossils)
  #eof
}

#' Map asymmetric fossil lineages
#'
#' Map fossils assuming asymmetric (budding) speciation.
#'
#' @param tree Phylo object.
#' @param fossil Fossils object.
#' @return An object of class fossils.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # add extant samples
#' f<-add.extant.occ(t, f, rho = 0.5)
#' # asymmetric mapping
#' f<-asymmetric.fossil.mapping(t, f)
#' @export
asymmetric.fossil.mapping<-function(tree,fossils){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(attr(fossils, "speciation") == "asymmetric")
    stop("Species have already been asymmetrically mapped")

  p<-asymmetric.identities(tree)

  af<-data.frame(h=numeric(),sp=numeric()) # h=horizon, sp=node

  for (i in 1:length(fossils[,2])) {

    f=fossils$sp[i]
    h=fossils$h[i]
    ai=p[f] # asym identity
    af<-rbind(af,data.frame(h=h,sp=ai))

  }
  af<-fossils(af, speciation.mode = "asymmetric", ages = attr(fossils, "ages"))
  return(af)
  # EOF
}

#' Identify attachment ages and extinction times in an incompletely sampled tree
#'
#' The age at which a species attaches to a tree may not be equivalent to the time of origin of a species
#' if sampling is incomplete.
#' This function takes an object of class fossils and the corresponding phylo object and calculates
#' the speciation (= attachment) times taking into account incomplete sampling.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param asymmetric.mapping If TRUE fossil sampling is assymmetric.
#' If speciation is "symmetric" and asymmetric.mapping = TRUE, the function calls asymmetric.fossil.mapping.
#' If speciation is "asymmetric" and asymmetric.mapping = FALSE, the function returns an error.
#' @return Dataframe containing the speciation & extinction times in an incompletely sampled tree.
#' @examples
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # add extant samples
#' f<-add.extant.occ(t, f, rho = 0.5)
#' # asymmetric mapping
#' f<-asymmetric.fossil.mapping(t, f)
#' # calculate attachment times
#' attachment.times(t, f)
#' @export
attachment.times<-function(tree,fossils,asymmetric.mapping=TRUE){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  # attachment identities & asymmetric / symmetric ages
  if(asymmetric.mapping){
    if(attr(fossils, "speciation") == "symmetric")
      fossils<-asymmetric.fossil.mapping(tree, fossils)
    attach.ident<-attachment.identities(tree,fossils)
    ages<-asymmetric.ages(tree)
  } else{
    if(attr(fossils, "speciation") == "asymmetric")
      stop("asymmetric.mapping = FALSE but speciation = \"asymmetric\"")
    warning("generating symmetric attachment times - this is experimental!\n")
    attach.ident<-symmetric.attachment.identities(tree,fossils)
    ages<-symmetric.ages(tree)
    nages<-n.ages(tree)
  }

  a<-data.frame(sp=numeric(),lineage.starts=numeric(),lineage.ends=numeric(),first.appearance=numeric())

  for(i in unique(fossils$sp)){

    # find the extinction time of i
    lineage.end=ages$end[which(ages$sp==i)]

    # find the attachment identity of i
    attaches=attach.ident$attaches[which(attach.ident$sp==i)]

    # find the speciation time of a
    if(asymmetric.mapping)
      lineage.start=ages$start[which(ages$sp==attaches)]
    else
      lineage.start=nages[which(names(nages)==attaches)]

    # Find the first appearance
    fa=max(fossils$h[which(fossils$sp==i)])

    a<-rbind(a,data.frame(sp=i,lineage.starts=lineage.start,lineage.ends=lineage.end,first.appearance=fa))
  }
  return(a)
  #EOF
}

### mixed speciation

#' Map asymmetric & symmetric lineages
#'
#' Asymmetric (or budding) speciation occurs with probability \eqn{\beta}. Asymmetric speciation gives rise to one new (morpho)species, while symmetric speciation gives rise to two new species and results in the extinction of the ancestor.
#' Note that if \eqn{\beta = 0} all speciation events will be asymmetric and if \eqn{\beta = 1} all speciation events will be symmetric.
#'
#' @param tree Phylo object.
#' @param beta Probability of symmetric speciation.
#' @examples
#' t<-ape::rtree(6)
#' mixed.speciation(t, 0.5)
#' @return
#' Dataframe of asymmetric or symmetric edge labels, parent edge labels and the mode of speciation.
#' "s" = symmetric speciation, "b" = budding (asymmetric) speciation, "o" = origin.
#' @export
mixed.speciation<-function(tree, beta){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!(beta >= 0 && beta <= 1))
    stop("beta must be a probability between 0 and 1")

  # sp = species; p = ancestor
  p<-data.frame(sp=numeric(),p=numeric(),equivalent.to=numeric(),mode=character())

  done = c()

  # identify the root
  root=length(tree$tip.label)+1

  p<-rbind(p,data.frame(sp=root,p=root,equivalent.to=root,mode="o"))

  ancestor = root

  process.complete = 0

  while(process.complete==0) {

    # fetch the two descendants
    row=which(tree$edge[,1]==ancestor)
    descendants=tree$edge[,2][row]
    d1<-descendants[1]
    d2<-descendants[2]

    if(runif(1) > (1 - beta)){
      # speciation event is symmetric

      if(!d2 %in% p[[1]]) {
        ai=p$equivalent.to[which(p$sp==ancestor)] # asymmetric identity
        p<-rbind(p,data.frame(sp=d1,p=ai,equivalent.to=d1,mode="s"))
        p<-rbind(p,data.frame(sp=d2,p=ai,equivalent.to=d2,mode="s"))
      }

    }
    else{
      # speciation event is asymmetric

      # unless d2 is not in the table (d1 also won't be in the table)
      if(!d2 %in% p[[1]]) {

        # find out how the parent of d1 is defined (e.g. as itself, or by an older ancestor)
        row=which(p$sp==ancestor) # one step ai=p$equivalent.to[which(p$sp==ancestor)]
        ai=p$equivalent.to[row] # asymmetric identity
        ap=p$p[which(p$sp==ai)]
        p<-rbind(p,data.frame(sp=d1,p=ap,equivalent.to=ai,mode="NA"))

        # define d2 its itself (e.g the budding champion)
        p<-rbind(p,data.frame(sp=d2,p=ai,equivalent.to=d2,mode="b"))

      }
    }

    # the following is simply a way of tranversing the tree in a particular order
    if(!d1 %in% done) {
      if ((is.tip(d1,tree)) == 1) {
        done<-c(done, d1)
      }
      else {
        ancestor=d1
      }
    }
    else if (!d2 %in% done) {
      if ((is.tip(d2,tree)) == 1) {
        done<-c(done, d2)
      }
      else {
        ancestor=d2
      }
    }
    else {
      if(ancestor==root){
        process.complete=1
      }
      else {
        done<-c(done, ancestor)
        row=which(tree$edge[,2]==ancestor)
        ancestor=tree$edge[,1][row]
      }
    }
  }
  return(p)
  # eof
}

#' Print out a data.frame with mixed speciation (asymmetric and symmetric) ages
#'
#' Function uses mixed.speciation() to assign asymmetric & symmetric lineages.
#'
#' @param tree Phylo object.
#' @param beta Probability of symmetric speciation.
#' @param root.edge If TRUE include root edge.
#' @return Dataframe with internal asymmetric or symmetric branch ages (min and max). Note that if root edge = T the oldest lineage incorporates the origin.
#' "s" = symmetric speciation, "b" = budding (asymmetric) speciation, "o" = origin.
#' @examples
#' t<-ape::rtree(6)
#' mixed.ages(t, 0.5)
#' @export
mixed.ages<-function(tree,beta,root.edge=T){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!(beta >= 0 && beta <= 1))
    stop("f must be a probability between 0 and 1")

  sym.ages<-symmetric.ages(tree, root.edge = root.edge)
  mixed.ident<-mixed.speciation(tree, beta = beta)
  origin=root(tree) # identify the root

  ages<-data.frame(sp=numeric(),p=numeric(),start=numeric(),end=numeric(),mode=character())

  for(i in unique(mixed.ident$equivalent.to)){

    # record the mode of speciation
    mode=mixed.ident$mode[which(mixed.ident$equivalent.to==i)][1]

    # record the parent
    parent=mixed.ident$p[which(mixed.ident$equivalent.to==i)][1]

    # identify all equivalent asymmetric lineages
    lineages=mixed.ident$sp[which(mixed.ident$equivalent.to==i)]

    # if i is the root and has no asymmetric descendants
    if ((length(lineages) < 2 ) && (!is.element(lineages, sym.ages$sp))){

      lineage.start=max(sym.ages$start)

      lineage.end=lineage.start

    }
    else {

      # find the oldest start time
      lineage.start=max(subset(sym.ages,sp %in% lineages)$start)

      # find the youngest end time
      lineage.end=min(subset(sym.ages,sp %in% lineages)$end)

    }

    ages<-rbind(ages,data.frame(sp=i,p=parent,start=lineage.start,end=lineage.end,mode=mode))

  }

  return(ages)
  #eof
}

#' Simulate anagenic species along branches
#'
#' @param ages Dataframe of branching times (for asymmetric and / or symmetric species).
#' @param lambda.a Rate of anagenic speciation.
#' @param parent.labels If TRUE print out parent edge labels.
#' @return Dataframe containing species labels, (morpho)species speciation & extinction times, and mode of speciation.
#' "a" = anagenic speciation, s" = symmetric speciation, "b" = budding (asymmetric) speciation, "o" = origin.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # assign symmetric and asymmetric species
#' sp1<-mixed.ages(t, 0.5)
#' # simulate anagenic species
#' sp2<-anagenic.species(sp1, 0.1)
#' @export
anagenic.species<-function(ages,lambda.a=0.1,parent.labels=FALSE){

  # identify the start of the new species id count
  species.count = max(ages$sp) + 1

  # sp = species id, p = parent, b = branch
  all.species<-data.frame(sp=numeric(),p=numeric(),b=numeric(),start=numeric(),end=numeric(),mode=character())

  for(i in 1:length(ages$sp)){

    sp = ages$sp[i]
    branch = ages$sp[i]
    parent = ages$p[i]
    mode = ages$mode[i]
    lineage.start = ages$start[i]
    lineage.end = ages$end[i]

    # work out total branch duration
    b = lineage.start - lineage.end

    # sample random number from a poisson distribution with rate = branch duriation x lambda.a
    rand = rpois(1, b * lambda.a)

    if (rand > 0) {
      # sample anageic speciation times from the branch duration
      h = runif(rand, min = lineage.end, max = lineage.start)

      # order and label the new species
      h = sort(h,decreasing = T)

      # deal with the original lineage
      all.species<-rbind(all.species,data.frame(sp=sp,p=parent,b=branch,start=lineage.start,end=h[1],mode=mode))

      for(i in 1:length(h)){

        start = h[i]
        if(i == length(h))
          end = lineage.end
        else
          end = h[i+1]

        if(i == 1)
          parent = sp
        else
          parent = species.count-1

        all.species<-rbind(all.species,data.frame(sp=species.count,p=parent,b=branch,start=start,end=end,mode="a"))

        species.count = species.count + 1
      }
    }
    else{
      all.species<-rbind(all.species,data.frame(sp=sp,p=parent,b=branch,start=lineage.start,end=lineage.end,mode=mode))
    }
  }

  # note at this stage the parent labels for sym and asym speciation events may be incorrect
  # the following loop assigns the correct labels if parent.labels = T
  # I made it optional because this might be slow for large trees
  if(parent.labels){
    for(i in 1:length(all.species$sp)){
      sp = all.species$sp[i]
      parent = all.species$p[i]
      mode = all.species$mode[i]
      if(mode == "s"){
        all.species$p[i] = max(all.species$sp[which(all.species$b==parent)])
      }
      else if(mode == "b") {
        birth.time = all.species$start[i]

        # identify parent branch members
        bm = all.species$sp[which(all.species$b==parent)]

        for(j in bm){
          bm.birth.time = all.species$start[which(all.species$sp==j)]
          bm.death.time = all.species$end[which(all.species$sp==j)]

          if((bm.birth.time > birth.time) && (bm.death.time < birth.time)){
            all.species$p[i] = j
          }
        }

      }
    }
  }
  else{
    all.species$p<-NULL
  }
  return(all.species)
  #eof
}

#' Generate datasets of cryptic species
#'
#' @param ages Dataframe of branching times for mixed species (asymmetric, symmetric, anagenic).
#' @param kappa Probability that a speciation event generates a cryptic species.
#'
#' @return Dataframe containing species labels, (morpho)species speciation & extinction times, mode of speciation, cryptic indicator, and corresponding cryptic speciation lables.
#' "a" = anagenic speciation, s" = symmetric speciation, "b" = budding (asymmetric) speciation, "o" = origin.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # assign symmetric and asymmetric species
#' sp1<-mixed.ages(t, 0.5)
#' # simulate anagenic species
#' sp2<-anagenic.species(sp1, 0.1, parent.labels = T)
#' assign cryptic speciation events
#' sp3<-cryptic.speciation(sp2, 0.5)
#' @export
cryptic.speciation<-function(ages, kappa){

  if(is.null(ages$p))
    stop("For cryptic speciation you must include parent labels")

  ages$cryptic = 0
  ages$cryptic.id = 0

  for(i in 1:length(ages$sp)){
    sp = ages$sp[i]
    if(runif(1) < kappa){
      # origin id is always = sp
      if(ages$mode[i] == "o")
        ages$cryptic.id[i] = sp
      else{
        # speciation event is cryptic
        ages$cryptic[i] = 1
        # identify parent
        parent = ages$p[i]
        # identify parent cryptic label
        parent.c = ages$cryptic.id[which(ages$sp == parent)]
        # assign cryptic species id
        ages$cryptic.id[i] = parent.c
        # and all other species that share the same cryptic id
        ages$cryptic.id[which(ages$cryptic.id == sp)] = parent.c
      }
    }else{
      ages$cryptic.id[i] = sp
    }
  }
  return(ages)
}

### tree functions

#' Create a data.frame with symmetric species durations (not node ages)
#
#' @param tree Phylo object.
#' @param root.edge If TRUE include root edge. Root edge takes the root label.
#' @return Dataframe with internal branch ages (min and max).
#' @examples
#' t<-ape::rtree(6)
#' symmetric.ages(t)
#' @export
# Function required by asymmetric.ages
symmetric.ages<-function(tree, root.edge=TRUE){

  node.ages=n.ages(tree)
  origin=root(tree) # identify the root

  ages<-data.frame(sp=numeric(),start=numeric(),end=numeric())

  for (i in tree$edge[,2]){ # internal nodes + tips

    # work out the max age of the lineage (e.g. when that lineage became extant)
    # & get ancestor
    row=which(tree$edge[,2]==i)
    ancestor=tree$edge[,1][row]

    # get the age of the ancestor
    a=which(names(node.ages)==ancestor)
    lineage.start=round(node.ages[[a]],7)

    # work out the min age of the lineage (e.g. when that lineage became extinct)
    # & get the branch length
    b=tree$edge.length[row]
    lineage.end=round(lineage.start-b,7) # branch length

    ages<-rbind(ages,data.frame(sp=i,start=lineage.start,end=lineage.end))
  }

  if(root.edge && exists("root.edge",tree) ){
    end = max(node.ages)
    start = end + tree$root.edge
    ages<-rbind(ages,data.frame(sp=origin,start=start,end=end))
  }

  return(ages)
  #eof
}

#' Create a data.frame with asymmetric species durations
#
#' @param tree Phylo object.
#' @param root.edge If TRUE include root edge.
#' @return Dataframe with internal asymmetric branch ages (min and max). Note that if root edge = TRUE the oldest lineage incorporates the origin.
#' @examples
#' t<-ape::rtree(6)
#' asymmetric.ages(t)
#' @export
asymmetric.ages<-function(tree,root.edge=TRUE){

  sym.ages<-symmetric.ages(tree, root.edge = FALSE) # note if root.edge = TRUE, the origin is handled below
  asym.ident<-asymmetric.identities(tree)
  origin=root(tree) # identify the root

  ages<-data.frame(sp=numeric(),start=numeric(),end=numeric())

  for(i in unique(asym.ident)){

    # identify all equivalent asymmetric lineages
    lineages=which(asym.ident==i)

    # find the oldest start time
    lineage.start=max(subset(sym.ages,sp %in% lineages)$start)

    # find the youngest end time
    lineage.end=min(subset(sym.ages,sp %in% lineages)$end)

    ages<-rbind(ages,data.frame(sp=i,start=lineage.start,end=lineage.end))
  }

  if(root.edge && exists("root.edge",tree) ){
    ages$start[which(ages$sp==origin)] = ages$start[which(ages$sp==origin)] + tree$root.edge
  }

  return(ages)
  #eof
}

#' Fetch descendant lineages in a symmetric tree
#
#' @param edge Edge label.
#' @param tree Phylo object.
#' @param return.edge.labels If TRUE return all descendant edge labels instead of tips.
#' @examples
#' t<-ape::rtree(6)
#' fetch.descendants(7,t)
#' fetch.descendants(7,t,return.edge.labels=TRUE)
#' @return
#' Vector of symmetric descendants
#' @export
fetch.descendants<-function(edge,tree,return.edge.labels=F){
  ancestor<-edge

  if(is.tip(edge, tree))
    return(NULL)

  # create a vector for nodes, tips & tracking descendent
  tips<-c()
  done<-data.frame(a=numeric()) # this data frame contains descendants (nodes+tips)

  coi=ancestor # clade of interest
  process.complete=0
  # count=0 # debugging code

  if(is.tip(ancestor,tree)){
    tip.label=tree$tip[ancestor]
    tips<-c(tips,tip.label)
  }
  else{

    while(process.complete==0) {

      # fetch the two descendants
      row=which(tree$edge[,1]==ancestor)
      descendants=tree$edge[,2][row]
      d1<-descendants[1]
      d2<-descendants[2]

      if(!d1 %in% done[[1]]) {
        if ((is.tip(d1,tree)) == 1) {
          done<-rbind(done,data.frame(a=d1))
          tip.label=tree$tip[d1]
          tips<-c(tips,tip.label)
        }
        else {
          ancestor=d1
        }
      }
      else if (!d2 %in% done[[1]]) {
        if ((is.tip(d2,tree)) == 1) {
          done<-rbind(done,data.frame(a=d2))
          tip.label=tree$tip[d2]
          tips<-c(tips,tip.label)
        }
        else {
          ancestor=d2
        }
      }
      else {
        if(ancestor==coi){
          process.complete=1
        }
        else {
          done<-rbind(done,data.frame(a=ancestor))
          row=which(tree$edge[,2]==ancestor)
          ancestor=tree$edge[,1][row]
        }
      }
      #	if (count==100) {
      #		process.complete=1
      #	}

      # count=count+1
    }
  }
  if(return.edge.labels)
    return(done$a)
  else
    return(tips)
  # EOF
}

# Fetch descendent lineages in an asymmetric tree
#
# @param tree Phylo object.
# @examples
# t<-ape::rtree(6)
# fetch.asymmetric.descendants(7,t)
# @return
# Vector of asymmetric descendants
# Function required by attachment.identities
fetch.asymmetric.descendants<-function(edge,tree){

  api<-asymmetric.parent.identities(tree)

  p<-data.frame(dec=numeric(),done=numeric())

  rows=which(api$parent==edge)
  if(length(rows)==0)
    return(NA) # no descendants
  decs=unique(api$equivalent.to[rows])

  p<-rbind(p,data.frame(dec=decs,done=0))

  process=0
  while(process==0){

    rows=which(p$done==0)
    if(length(rows)==0)
      process=1

    ancs=p$dec[rows]
    for(a in ancs){
      rows=which(api$parent==a)
      if(length(rows) > 0) {
        decs=unique(api$equivalent.to[rows])
        p<-rbind(p,data.frame(dec=decs,done=0))
      }
      p$done[which(p$dec==a)]=1
    }

  }
  return(p$dec)
  #eof
}

### hidden functions

### tree + fossils

# Identify attachment lineages in an incomplete asymmetrically mapped tree
#
# @param tree Phylo object.
# @param fossils Fossils object.
# @return Dataframe containing edge labels and attachment IDs.
# @examples
# # simulate tree
# t<-ape::rtree(6)
# # simulate fossils
# f<-sim.fossils.poisson(t, 2)
# # asymmetric fossil mapping
# f<-asymmetric.fossil.mapping(t, f)
# # asymmetric attachment identities
# attachment.identities(t, f)
# Function required by attachment.times
attachment.identities<-function(tree,fossils) {

  if(attr(fossils, "speciation") != "asymmetric")
    stop("Fossils speciation mode != asymmetric")

  api<-asymmetric.parent.identities(tree)

  # identify the root
  root=length(tree$tip.label)+1

  p<-data.frame(sp=numeric(),attaches=numeric())

  for(i in unique(fossils$sp)){

    # if am I the root
    if(i == root){
      p<-rbind(p,data.frame(sp=i,attaches=i))
    }
    else {
      # identify the asymmetric parent
      parent=api$parent[which(api$species==i)]
      j=i

      process=0
      while(process==0){

        if(parent == root){
          # if the root has been sampled
          if(root %in% fossils$sp){
            p<-rbind(p,data.frame(sp=i,attaches=j))
          }
          else{

            # fetch the root descendants
            root.decs=fetch.asymmetric.descendants(root,tree)
            # if i is not the only other sampled descendant
            if(length(which(root.decs[!root.decs==i] %in% fossils$sp)) > 0){

              # identify decs also sampled
              decs=root.decs[which(root.decs %in% fossils$sp)]

              # identify the right most sample
              row=min(which(tree$edge[,2] %in% decs))
              right=tree$edge[,2][row]

              # if i/j is not the right most sample in the sym tree
              if(i != right){
                p<-rbind(p,data.frame(sp=i,attaches=j))
              }
              else{
                p<-rbind(p,data.frame(sp=i,attaches=root))
              }
            }
            else{
              p<-rbind(p,data.frame(sp=i,attaches=root))
            }
          }
          process=1
        }

        # if parent is sampled
        # -> attachment identity = self (i) or nearest ancestor (j)
        else if(parent %in% fossils$sp){
          p<-rbind(p,data.frame(sp=i,attaches=j))
          process=1
        }

        else {
          decs=fetch.asymmetric.descendants(parent,tree)
          # if i is not the only other sampled descendant
          if ( length(which(decs[!decs==i] %in% fossils$sp)) > 0 ){

            # identify decs also sampled
            decs=decs[which(decs %in% fossils$sp)]

            # identify the right most sample
            # it will appear highest in the table
            row=min(which(tree$edge[,2] %in% decs))
            right=tree$edge[,2][row]

            # if i/j is not the right most sample in the sym tree
            # -> attachment identity = self i/j
            #if(i == left){
            if(i != right){
              p<-rbind(p,data.frame(sp=i,attaches=j))
              process=1
            }
            else{
              j=parent
              parent=api$parent[which(api$species==j)]
            }
          }
          else{
            j=parent
            parent=api$parent[which(api$species==j)]
          }
        }
      }
    }
  }
  return(p)
  #eof
}

# Identify attachment lineages in an incomplete symmetrically mapped tree
#
# @param tree Phylo object.
# @param fossils Fossil dataframe.
#
symmetric.attachment.identities<-function(tree, fossils){

  p<-data.frame(sp=numeric(),attaches=numeric())

  for(i in unique(fossils$sp)){
    if(i == root(tree)){
      p<-rbind(p,data.frame(sp=i,attaches=i))
      next
    }
    a = ancestor(i, tree)
    if(a %in% fossils$sp)
      p<-rbind(p,data.frame(sp=i,attaches=a))
    else{
      j = i
      process = 0
      while(process == 0){
        # fetch immediate descendants
        d = descendants(a, tree)
        # identify the other descendant
        d2 = d[which(d != j)]
        # if the other immediate descendant is sampled
        if(d2 %in% fossils$sp){
          p<-rbind(p,data.frame(sp=i,attaches=a))
          process = 1
          next
        }
        # fetch all d2 descendants
        dd = fetch.descendants(d2, tree, return.edge.labels = T)
        # if any dd species are sampled, node j exists in the tree
        df = which(dd %in% fossils$sp)
        #if(any(dd == fossils$sp)){
        if(length(df) > 0){
          p<-rbind(p,data.frame(sp=i,attaches=a))
          process = 1
        }
        else{
          j = a
          a = ancestor(j, tree)
        }
      }
    }
  }
  return(p)
  #eof
}

### tree functions

# Sample asymmetric lineages
#
# @param tree Phylo object.
# @return
# Vector of asymmetric edge labels.
# @examples
# t<-ape::rtree(6)
# asymmetric.identities(t)
asymmetric.identities<-function(tree){

  parents = rep(0, length(tree$tip.label)+tree$Nnode)
  # identify the root
  root=length(tree$tip.label)+1
  parents[root] = root

  aux = function(node, par) {
    desc = tree$edge[which(tree$edge[,1] == node), 2]
    if(length(desc) == 0) return(par)
    par[desc[1]] = par[node]
    par = aux(desc[1], par)
    par[desc[2]] = desc[2]
    par = aux(desc[2], par)
    par
  }

  parents = aux(root, parents)
  parents
}

# Identify parent lineages in an asymmetric tree
#
# @param tree Phylo object.
# @return
# @examples
# t<-ape::rtree(6)
# asymmetric.parent.identities(t)
# Dataframe of asymmetric parent labels.
# Function required by attachment.identities
asymmetric.parent.identities<-function(tree){

  p<-data.frame(species=numeric(),equivalent.to=numeric(),parent=numeric())

  done<-data.frame(a=numeric())

  # identify the root
  root=length(tree$tip.label)+1
  parent=NA

  p<-rbind(p,data.frame(species=root,equivalent.to=root,parent=parent))

  ancestor=root

  process.complete=0

  while(process.complete==0) {

    # fetch the two descendants
    row=which(tree$edge[,1]==ancestor)
    descendants=tree$edge[,2][row]
    d1<-descendants[1]
    d2<-descendants[2]

    # unless d2 is not in the table (d1 also won't be in the table)
    if(!d2 %in% p[[1]]) {

      # find out how the parent of d1 is defined (e.g. as itself, or by an older ancestor)
      row=which(p$species==ancestor) # one step ai=p$equivalent.to[which(p$species==ancestor)]
      ai=p$equivalent.to[row] # asymmetric identity
      if(ai==root)
        pi=NA
      else {
        pi=ancestor(ai,tree)
        pi=p$equivalent.to[which(p$species==pi)]
      }
      p<-rbind(p,data.frame(species=d1,equivalent.to=ai,parent=pi)) # this is incorrect

      # define d2 as its itself (e.g the budding champion = new species)
      p<-rbind(p,data.frame(species=d2,equivalent.to=d2,parent=ai)) # this is correct

    }

    # the following is simply a way of tranversing the tree in a particular order
    if(!d1 %in% done[[1]]) {
      if ((is.tip(d1,tree)) == 1) {
        done<-rbind(done,data.frame(a=d1))
      }
      else {
        ancestor=d1
      }
    }
    else if (!d2 %in% done[[1]]) {
      if ((is.tip(d2,tree)) == 1) {
        done<-rbind(done,data.frame(a=d2))
      }
      else {
        ancestor=d2
      }
    }
    else {
      if(ancestor==root){
        process.complete=1
      }
      else {
        done<-rbind(done,data.frame(a=ancestor))
        row=which(tree$edge[,2]==ancestor)
        ancestor=tree$edge[,1][row]
      }

    }
    #	if (count==100) {
    #		process.complete=1
    #	}

    # count=count+1
  }

  #eof  # the number of asym lineages should be 2n
  return(p)
}


# Assign species labels to asymmetric lineages
#
# @param tree Phylo object.
# @return
# Vector of asymmetric species labels.
# @examples
# t<-ape::rtree(6)
# asymmetric.species.identities(t)
asymmetric.species.identities<-function(tree){

  children = rep(0, length(tree$tip.label)+tree$Nnode)
  # identify the root
  root=length(tree$tip.label)+1

  aux = function(node, par) {
    desc = tree$edge[which(tree$edge[,1] == node), 2]
    if(length(desc) == 0) {
      par[node] = node
      return(par)
    }
    par = aux(desc[1], par)
    par = aux(desc[2], par)
    par[node] = par[desc[1]]
    par
  }

  children = aux(root, children)
  children
}
