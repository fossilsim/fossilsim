% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsampling.R
\name{prune.fossil.tips}
\alias{prune.fossil.tips}
\title{Remove fossil lineages from a tree}
\usage{
prune.fossil.tips(tree)
}
\arguments{
\item{tree}{an object of class "Phylo".}
}
\value{
an object of class "Phylo". If fossil lineages were found in the tree
  these will be pruned, if not then the original tree is returned.
}
\description{
Remove fossil lineages from a tree
}