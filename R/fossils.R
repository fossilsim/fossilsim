
# define class fossils
fossils<-function(f, speciation.mode = "asymmetric"){
  # do some checking: are sp and h available?
  me <- f
  class(me) <- c("fossils", class(me))
  attr(me, "speciation") <- speciation.mode
  return(me)
}

print.fossils<-function(f, max.length = 10, ...){
  print(as.data.frame(f)[1:max.length,])
  if(length(f$sp) > max.length)
    cat("...\n")
  cat("Fossil record with", length(f$h), "occurrences repesenting", length(unique(f$sp)), "species\n")
  if(!is.null(attr(f,"speciation")))
    cat("Speciation mode:", attr(f, "speciation"))
}



# as.fossil.record<-function(f,...){
#   if(identical(class(f),"fossil.record"))
#     return(f)
#   else
#     UseMethod("as.fossil.record")
# }
