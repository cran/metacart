#' Summary of the results of a FE meta-tree object
#'
#' @param object fitted tree of class \code{FEmrt}.
#' @param digits specified number of decimals in the printed results.
#' @param \dots additional arguments to be passed.
#' @return Summary of a FE meta-tree
#' @details If no moderator effect is detected,
#' the summary function will show the standard meta-analysis results.
#' Otherwise, the summary function will show the subgroup meta-analysis results,
#' with the significance test results for moderator effects, the split points of the moderators,
#' and the estimated subgroup summary effect sizes.
#' @importFrom stats symnum
#' @export
summary.FEmrt <- function(object, digits = 3, ...){
  if (!is.element("FEmrt", class(object))) {
    stop("Argument 'object' must be an object of class \"FEtree\".")
  } else {
    if (length(object$n) == 1) {
      cat("\n")
      cat("Fixed Effects meta-tree (K = ", sum(object$n), " studies); ",
          sep = "")
      cat("\n")
      print(object$call)
      cat("\n")
      cat("No moderator effect was detected" )
      cat("\n\n")
      cat("Test for Heterogeneity under FE assumption:")
      cat("\n")
      Q <- formatC(object$Q, digits=digits, format="f")
      pval <- format.pval(object$pval, eps = 10^(-digits-1))
      cat("Q = ", object$Q," (df = ", object$df, "), ", 
          ifelse(is.null(object$pval.perm),"p-value ", "Permutation p-value"), ifelse( is.null(object$pval.perm) ,pval, format.pval(object$pval.perm, eps = 10^(-digits-1)))
          , ";", sep = "")
      cat("\n\n")
      cat("Fixed Effect Meta-analysis Results:")
      cat("\n")
      sig <- symnum(object$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
      res.table <- c(object$n,
                     formatC(c(object$g, object$se, object$zval, object$pval, object$ci.lb, object$ci.ub), digits, format = "f"),
                     sig)
      names(res.table) <- c("no.", "Est.", "se", "zval", "pval", "ci.lb", "ci.ub", " ")
      print(res.table, quote = FALSE, right = TRUE, ...)
      cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")
      
    } else {
      cat("\n")
      cat("Fixed Effects meta-tree (K = ", sum(object$n), " studies); ",
          sep = "")
      cat("\n")
      print(object$call)
      cat("\n")
      cat("A tree with ", length(object$n), " terminal nodes was detected", sep="" )
      cat("\n")
      cat("Moderators were detected as: ", paste(as.character(object$moderators), collapse = ", "), sep="")
      cat("\n\n")
      cat("Test for Between-Subgroups Heterogeneity under FE assumption:")
      cat("\n")
      Qb <- formatC(object$Qb, digits=digits, format="f")
      pval.Qb <- format.pval(object$pval.Qb, eps = 10^(-digits-1))
      cat("Qb = ", Qb," (df = ", object$df, "), ",
          ifelse(is.null(object$pval.perm),"p-value ", "Permutation p-value"), ifelse( is.null(object$pval.perm) ,pval.Qb, format.pval(object$pval.perm, eps = 10^(-digits-1)))
          , ";", sep = "")
      cat("\n\n")
      cat("Subgroup Meta-analysis Results:")
      cat("\n")
      sig <- symnum(object$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
      res.table <- cbind(object$n,
                         formatC(cbind(object$Qw, object$g, object$se, object$zval, object$pval, object$ci.lb, object$ci.ub), digits, format = "f"),
                         sig)
      colnames(res.table) <- c( "K", "Qw", "Est.", "se", "zval",  "pval",
                                "ci.lb", "ci.ub", " ")
      rownames(res.table) <- rownames(object$tree$frame[object$tree$frame$var == "<leaf>", ])
      print(res.table, quote = FALSE, right = TRUE, ...)
      cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")
    }
    
  }
}


#' Summary of the results of a RE meta-tree object
#'
#' @param object fitted tree of class \code{REmrt}.
#' @param digits specified number of decimals in the printed results.
#' @param \dots additional arguments to be passed.
#' @return Summary of a RE meta-tree
#' @details If no moderator effect is detected,
#' the summary function will show the standard meta-analysis results.
#' Otherwise, the summary function will show the subgroup meta-analysis results,
#' with the significance test results for moderator effects, the split points of the moderators,
#' and the estimated subgroup summary effect sizes.
#' @importFrom stats symnum
#' @export
summary.REmrt <- function(object, digits = 3, ...){
  if (!is.element("REmrt", class(object))) {
    stop("Argument 'object' must be an object of class \"REmrt\".")
  } else {
    
    if(!is.null(object$pruned)){
      if(object$pruned==FALSE){
        
      cat("\n")
      cat("Random Effects meta-tree (K = ", sum(object$initial), " studies); ",
          sep = "")
      cat("\n")
      print(object$call)
      cat("\n")
      cat("A tree with ", length(object$initial), " terminal nodes was detected", sep="" )
      cat("\n")
      cat("Moderators were detected as: ", paste(as.character(object$mod[!is.na(object$mod)]), collapse = ", "), sep="")
      cat("\n\n")
      cat("Between-Subgroups Heterogeneity under RE assumption:")
      cat("\n")
      Qb <- formatC(object$Qb, digits=digits, format="f")
      cat("Qb = ", Qb[length(Qb)])
      cat("\n")
      cat("Test for Between-Subgroups Heterogeneity under RE assumption is not performed.")
      cat("\n\n")
      
      
      if(!is.null(object$boot_tau2)){
        tau2 <- formatC(object$boot_tau2, digits=digits, format="f")
        cat("The estimate using bootstrap for the residual heterogeneity tau2 = ", tau2[length(tau2)], ";", sep = "" )
        
      }else{
        tau2 <- formatC(object$tau2, digits=digits, format="f")
        cat("The estimate without bootstrap for the residual heterogeneity tau2 = ", tau2[length(tau2)], ";", sep = "" )
      }
      
      cat("\n\n")
        
      }
    }else{
    
    if (length(object$n) == 1) {
      cat("\n")
      cat("Random Effects meta-tree (K = ", sum(object$n), " studies); ",
          sep = "")
      cat("\n")
      print(object$call)
      cat("\n")
      cat("No moderator effect was detected" )
      cat("\n\n")
      cat("Test for Heterogeneity")
      cat("\n")
      Q <- formatC(object$Q, digits=digits, format="f")
      pval <- format.pval(object$pval, eps = 10^(-digits-1))
      cat("Q = ", object$Q," (df = ", object$df, "), ",
          ifelse(is.null(object$pval.perm),"p-value ", "Permutation p-value"), ifelse( is.null(object$pval.perm) ,pval, format.pval(object$pval.perm, eps = 10^(-digits-1)))
          , ";", sep = "")
      cat("\n")
      if(!is.null(object$boot_tau2)){
        tau2 <- formatC(object$boot_tau2, digits=digits, format="f")
        cat("The estimate using bootstrap for the residual heterogeneity tau2 = ", tau2, ";", sep = "" )
        
      }else{
        tau2 <- formatC(object$tau2, digits=digits, format="f")
        cat("The estimate without bootstrap for the residual heterogeneity tau2 = ", tau2, ";", sep = "" )
      }
      
      cat("\n\n")
      if(!is.null(object$boot_res)){
        cat("Random Effects Meta-analysis Results with Bootstrap:")
        cat("\n")
        res.table <- c(K=object$n,
                       formatC(object$boot_res, digits, format = "f"))
        print(res.table, quote = FALSE, right = TRUE, ...)

      }else{
        cat("Random Effects Meta-analysis Results without Bootstrap:")
        cat("\n")
        sig <- symnum(object$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " "))
        res.table <- c(object$n,
                       formatC(c(object$g, object$se, object$zval, object$pval, object$ci.lb, object$ci.ub), digits, format = "f"),
                       sig)
        names(res.table) <- c("K", "est", "se", "zval", "pval", "ci.lb", "ci.ub", " ")
        print(res.table, quote = FALSE, right = TRUE, ...)
        cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")
      }

      
      
    } else {
      cat("\n")
      cat("Random Effects meta-tree (K = ", sum(object$n), " studies); ",
          sep = "")
      cat("\n")
      print(object$call)
      cat("\n")
      cat("A tree with ", length(object$n), " terminal nodes was detected", sep="" )
      cat("\n")
      cat("Moderators were detected as: ", paste(as.character(object$moderators), collapse = ", "), sep="")
      cat("\n\n")
      cat("Test for Between-Subgroups Heterogeneity under RE assumption:")
      cat("\n")
      Qb <- formatC(object$Qb, digits=digits, format="f")
      pval.Qb <- format.pval(object$pval.Qb, eps = 10^(-digits-1))
      cat("Qb = ", Qb," (df = ", object$df, "), ",
          ifelse(is.null(object$pval.perm),"p-value ", "Permutation p-value "), ifelse( is.null(object$pval.perm) ,pval.Qb, format.pval(object$pval.perm, eps = 10^(-digits-1)))
          , ";", sep = "")
      cat("\n")
      
      if(!is.null(object$boot_tau2)){
        tau2 <- formatC(object$boot_tau2, digits=digits, format="f")
        cat("The estimate using bootstrap for the residual heterogeneity tau2 = ", tau2, ";", sep = "" )
        
      }else{
        tau2 <- formatC(object$tau2, digits=digits, format="f")
        cat("The estimate without bootstrap for the residual heterogeneity tau2 = ", tau2, ";", sep = "" )
      }
      
      cat("\n\n")
      
      
      if(!is.null(object$boot_res)){
        cat("Random Effects Meta-analysis Results with Bootstrap:")
        cat("\n")
        res.table <- cbind(K=object$n, formatC(object$boot_res, digits, format = "f"))
        print(res.table, quote = FALSE, right = TRUE, ...)
        
      }else{
        cat("Random Effects Meta-analysis Results without Bootstrap:")
        cat("\n")
        sig <- symnum(object$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " "))
        res.table <- cbind(object$n,
                           formatC(cbind(object$g, object$se, object$zval, object$pval, object$ci.lb, object$ci.ub), digits, format = "f"),
                           sig)
        colnames(res.table) <- c( "K", "est", "se", "zval",  "pval",
                                  "ci.lb", "ci.ub", " ")
        print(res.table, quote = FALSE, right = TRUE, ...)
        cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")
      }
      
      
    }
    }  
  }
}
