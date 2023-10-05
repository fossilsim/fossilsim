setClass("fbdrange", contains = "treedata")

#TODO conditional dependency on treeio
get_fbdrange_from_file = function(input_file) {
  tree <- treeio::read.beast(input.treefile)
  as.fbdrange(tree)
}

as.fbdrange.treedata = function(data, ...) {
  if(is.null(data@data$range) || is.null(data@data$orientation)) stop("Fbdrange objects require range and orientation data")
  fbd = new("fbdrange", data)
}

as.fbdrange.default = function(data, ...) {
  fbdrange(data, ...)
}

#TODO checks that fbd range object contains elements "range" and "orientation" as part of data