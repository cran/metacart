#' Print function for FEmrt
#'
#' Print the results of a FEmrt object
#'
#' @param x fitted tree of class \code{FEmrt}.
#' @param \dots additional arguments to be passed.
#' @return Printed output of a FE meta-tree
#' @details
#' The function returns the objects concerning the analysis results.
#' @export
print.FEmrt<- function(x, ...){
  if (length(x$n) < 2) {
    cat("\n")
    cat("Fixed Effects Meta-Tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )
    cat("\n")
    cat("use summary() to see the meta-analysis results")


  } else {
    cat("\n")
    cat("Fixed Effects Meta-tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("A tree with ", length(x$n), " terminal nodes was detected", sep="" )
    cat("\n")
    cat("The moderators are ", paste(as.character(x$moderators), collapse = ", "), sep = "")
    cat("\n")
    cat("Use summary() and plot() to inspect the moderator analysis results and the tree structure.")
    cat("\n")

    #Rounding change:
    x$tree$splits[,4L]<-round(x$tree$splits[,4L],2)
      
    print(x$tree)
    #print(x$tree)
  }
}


#' Print function for REmrt
#'
#' Print the results of a REmrt object
#'
#' @param x fitted tree of class \code{FEmrt}.
#' @param \dots additional arguments to be passed.
#' @return Printed output of a RE meta-tree
#' @details
#' The function returns the results (e.g., the value of the Q-between) after each split of the tree.
#' @export
print.REmrt<- function(x, ...){
  
  if(!is.null(x$pruned)){
    if(x$pruned==FALSE){
      cat("\n")
      cat("The following tree is the initial tree and it has not been pruned yet.")
      cat("\n")
      cat("Random Effects Meta-tree (K = ", sum(x$initial), " studies); ",
          sep = "")
      cat("\n")
      print(x$call)
      cat("\n")
      cat("A tree with ", length(x$initial), " terminal nodes was detected", sep="" )
      cat("\n")
      cat("The moderators are ", paste(as.character(x$mod[!is.na(x$mod)]), collapse = ", "), sep = "")
      cat("\n")
      cat("use summary() and plot() to see the moderator analysis results and the tree structure")
      cat("\n")
      cat("\n")
      #Rounding change:
      rows<-grep("<", x$tree$split)
      for (i in rows) {
        x$tree$split[i]<-paste0(substr(x$tree$split[i], 1, unlist(gregexpr('<', x$tree$split[i]))+1), " " ,round(as.numeric(substr(x$tree$split[i], unlist(gregexpr('<', x$tree$split[i]))+2, nchar(x$tree$split[i]))),2))
      }
      
      print(x$tree)
      
      
      
    }
  }else{
  
  if (length(x$n) < 2) {
    cat("\n")
    cat("Random Effects Meta-Tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )
    cat("\n")
    cat("Use summary() to inspect the meta-analysis results")


  } else {
    cat("\n")
    cat("Random Effects Meta-tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("A tree with ", length(x$n), " terminal nodes was detected", sep="" )
    cat("\n")
    cat("The moderators are ", paste(as.character(x$moderators), collapse = ", "), sep = "")
    cat("\n")
    cat("use summary() and plot() to see the moderator analysis results and the tree structure")
    cat("\n")

    #Rounding change:
    rows<-grep("<", x$tree$split)
    for (i in rows) {
      x$tree$split[i]<-paste0(substr(x$tree$split[i], 1, unlist(gregexpr('<', x$tree$split[i]))+1), " " ,round(as.numeric(substr(x$tree$split[i], unlist(gregexpr('<', x$tree$split[i]))+2, nchar(x$tree$split[i]))),2))
    }
    print(x$tree)
    #print(x$tree)
  }
  }
}
