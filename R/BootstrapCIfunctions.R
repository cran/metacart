#Functions needed in metacart
# =========================================
# SEND down MC 1
# =========================================
send.down.MC <- function(fit) { #, newdata){ 
  #pred.new <- predict(fit, newdata) #works
  #newdata$Tnodes <- as.factor(pred.new$TNodes) #works
  #res<- metafor::rma(vi = vi, y = g, mods = ~ Tnodes-1, method = "DL",
  # data = newdata) #works
  train <- data.frame(node = names(fit$n), no.train = fit$n, 
                      est.train = fit$g, se.train = fit$se)#works
  #test <- data.frame(node =  as.character(sort(unique(newdata$Tnodes))), #works
  #  no.test = tapply(newdata$g, newdata$Tnodes, length), 
  #  est.test = res$b,
  #  se.test = res$se)
  output <- train
  #output <- merge(train, test, by="node", all.x = FALSE) #works
  tau2 <- fit$tau2
  #tau2 <- c(fit$tau2, res$tau2) #works
  #names(tau2) <- c("tau2.train", "tau2.test") #works
  list(output = output, tau2 = tau2) #works
  
}

# =========================================
# SEND down MC 2
# =========================================
send.down.MC2 <- function(fit){ # ,newdata){
  #pred.new <- predict(fit, newdata)
  #newdata$Tnodes <- as.factor(pred.new$TNodes)
  # (take the g (or estimate) ^ 2, devided by sampling variance per number of 
  # terminal nodes) minus
  #testQ <- tapply(newdata$g^2/newdata$vi, newdata$Tnodes, sum) -  
  # (Take the g (or estimate), devided by sampling per terminal nodes together) 
  # ^ 2 devided by
  #tapply(newdata$g/newdata$vi, newdata$Tnodes, function(x)(sum(x))^2)/ 
  # 1 / sampling variance (weights) per terminal node together
  #tapply(1/newdata$vi, newdata$Tnodes, sum) 
  olddata <- fit$data #old data
  #make sure that output data always have the outcome in the first column
  olddata$g <- olddata[ ,1]
  #same but with train data
  trainQ <- tapply(olddata$g^2/olddata$`(vi)`, olddata$term.node, sum) - 
    tapply(olddata$g/olddata$`(vi)`, olddata$term.node, function(x)(sum(x))^2)/
    tapply(1/olddata$`(vi)`, olddata$term.node, sum)
  #number of tree with number of terminal node, number of train, Q statistic of
  #this number of terminal nodes
  train <- data.frame(node = names(fit$n), no.train = fit$n, 
                      Q.train = trainQ)
  
  train$Q.by.n.train <- train$Q.train/train$no.train # Q statistic of trees with 
  # different terminal nodes/ devided by the number of trees with certain terminal 
  # nodes (mean)
  #test <- data.frame(node =  names(testQ), #same for test
  # no.test = tapply(newdata$g, newdata$Tnodes, length),
  # Q.test = testQ)
  #test$Q.by.n.test <- test$Q.test/test$no.test #same for test
  output <- train
  #output <- merge(train, test, by="node", all.x = FALSE) #together in one dataset
  output
  
}

###############################################################################
########################## function for bootstrapping ##########################
Bootstrap_bias_correction <- function(fit = fit, data = data, B = 50) { 
  
  colnames(data)[which(colnames(data)=="(vi)")]<-"vi" #Change the name from (`(vi)` to vi). Required for the REmrt function below
  
  #add fitted (original) tree, test data, original dataset & number of bootstraps
  out.Q0 <- send.down.MC2(fit) #, test) 
  K <- sum(fit$n)
  node0 <- fit$data$term.node #take every terminal node of study of fitted tree
  bias <- 0 #set bias to zero
  formula <- fit$formula
  vi <- fit$call$vi
  c.pruning <- fit$call$c.pruning

  maxL <-ifelse(is.null(fit$call$maxL),5, fit$call$maxL)
  minbucket <- ifelse(is.null(fit$call$minbucket),3, fit$call$minbucket)
  minsplit <- ifelse(is.null(fit$call$minsplit),6, fit$call$minsplit)
  
  skip.c <- 0
  
  
  for (b in 1:B) { #start bootstrapping
    #take random studies of dataset with replacement 
    id.b <- sample(1:K, size=K, replace=TRUE)  
    dat.b <- data[id.b,] #select the selected random 120 studies of the data file
    # fit new RE tree with same specifics as fit
    fit.b <- REmrt(formula, data = dat.b, vi = vi, c.pruning = c.pruning,
                   maxL = maxL, minbucket = minbucket, minsplit = minsplit)
    #JVM when bootstrapped tree has no moderators
    if (length(fit.b$moderators) == 0) {
      skip.c <- skip.c + 1
      next
    }
    #get new dataset with specifics of train/test with trees with terminal nodes
    info.b <- send.down.MC2(fit.b)
    #get Q statistic of trees with different terminal nodes & test - train
    #bias.b <- info.b$Q.by.n.test - info.b$Q.by.n.train 
    bias.b <- info.b$Q.by.n.train 
    #every value below zero becomes zero (ED: why?)
    positive.bias = TRUE
    if (positive.bias) bias.b <- pmax(bias.b, 0) 
    #predict number of terminal nodes based on new fitted value in REtree

    predobj<- predict(fit.b, data)
    node.b <- predobj$TNodes #does not work yet in version 1_3.3
    #generate cross table with "original" number of trees with certain
    #terminal nodes vs test number of trees with certain terminal nodes
    tab <- table(node0, node.b) 
    M.prop <- round(prop.table(tab, 1),4) #add round JVM, proportional table
    #multiplication of every row with bias of "original" terminal tree
    bias.b <- M.prop%*%bias.b 
    bias <- bias + bias.b #add bias of every bootstrap
  }
  
  newB = B - skip.c #recalculation of the performed bootstrapped trees
  meanbias <- bias/newB #divide by the number of bootstraps
  out.Q0$Q.c <- (out.Q0$Q.by.n.train + meanbias[,1])*out.Q0$no.train
  output <- list(table = out.Q0, 
                 total_b =  bias,
                 mean_bias  =  meanbias,
                 perf_b = newB)
  
  ##########recalculate average mean squares and add to output
  output$table$out.after.Qc <- output$table$Q.c/output$table$no.train
  return(output)
}


#------------------------------------------------------------------------------

################################Tau function####################################

#add fitted tree, test data, bootstrap data & original data
Recalc_tau_function <- function(fit = fit, Bootstrapoutput, data = data ) { 
  out.info0 <- send.down.MC(fit)  
  node0 <- fit$data$term.node #take every terminal node of study of fitted tree 
  ### equation 8 in manuscript & equation 4 in own
  C <- tapply(1/data$`(vi)`, node0, sum) -  #1/ sampling variance of list of trees 
    #with number of terminal nodes & sumup, minus
    tapply(1/data$`(vi)`, node0, function(x) sum(x^2))/ # 1/ sampling variance of 
    #list of trees with no term.nodes& sum up and square (positive numbers?)
    tapply(1/data$`(vi)`, node0, sum) #devided by 1/ sampling variance of list of
  #trees with number of terminal nodes (weights?)
  df <- Bootstrapoutput$table$no.train-1 #degrees of free per leaf
  ### equation 7 in manuscript
  tau2.co <- (sum(Bootstrapoutput$table$Q.c) - sum(df))/sum(C) #(sum every 
  #corrected Q (calculated with bias)) minus ((the sum of of degrees of freedom)) 
  #/ sum of components
  output <-  list(tau.co = tau2.co, table = out.info0)
  output
}

################################################################################


### SE Function ################################################################

#add fitted data and calculated tau output
Recalc_SE_function <- function(fit = fit, tau_out = tau_out, data = data) { 
  out.info0 <- tau_out$table
  node0 <- fit$data$term.node #take every terminal node of study of fitted tree 
  #orginal Se per leaf
  
  #se.test0 <- out.info0$output$se.test*sqrt(out.info0$output$no.test) 
  # standerd error of test * square root of number of test tree with certain nodes
  se.train0 <- out.info0$output$se.train*sqrt(out.info0$output$no.train) 
  # standerd error of train * square root of number of test tree with certain nodes
  #corrected SE per leaf, equation 5 in manuscript

  train.co <- data.frame(g.co = tapply(fit$data[,1]/(data$`(vi)` + tau_out$tau.co), node0, sum)/ 
                           #sum of g/vi per tree with number of nodes
                           tapply(1/(data$`(vi)` + tau_out$tau.co), node0, sum), 
                         #divided by weights per tree with number of nodes
                         se.co = tapply(data$`(vi)` + tau_out$tau.co, node0, 
                                        function(x) sqrt(1/sum(1/x)))) 
  #generate new dataset with corrected effect size & standerd error per tree with certain nodes
  #train.co$se.co <- train.co$se.co*sqrt(out.info0$output$no.train)
  output <- list(se.train = se.train0, train.co = train.co)
  #output <- list(se.test = se.test0, se.train = se.train0, train.co = train.co, se.train.co = se.train.co)
}



#' This functions performs bootstrap to compute the confidence intervals for the subgroup effect size estimates.
#' This function is only applicable to Random effects metaregression trees with 2 terminal nodes or more.
#' 
#' @param Metatree fitted tree of class \code{REmrt}.
#' @param nboot number of bootstrap samples.
#' 
#' @return tree containing the input tree, the Bootstrap estimates for the effect sizes and standard errors, Bootstrap estimate for tau2, and the Bootstrap bias correction.
#' 
#' @examples 
#' set.seed(12345) 
#' data(dat.BCT2009)
#' library(Rcpp)
#' REtree <- REmrt(g ~ T1 + T2+ T4 +T25, vi = vi, data = dat.BCT2009, c.pruning = 0)
#' BootTree<-BootCI(REtree, nboot = 3)
#' summary(BootTree)
#' 
#' @importFrom stats qnorm
#' @export
BootCI <- function(Metatree, nboot=50){
  suppressWarnings({
    
  Metadata<-Metatree$data
  #do we 
  fit <- Metatree
  #JVM when tree has no moderators
  if (length(fit$moderators) == 0) {
    res <- 0
  } else {
    Boot_out <- Bootstrap_bias_correction(fit = fit, data = Metadata, B = nboot) #works without test
    tau_out <- Recalc_tau_function(fit, Boot_out, Metadata) #works without test
    new_tau <- tau_out$tau.co #new tau
    SE_out <- Recalc_SE_function(fit, tau_out, Metadata) #work without test
    res_corrected <- data.frame(SE_out$train.co)
    ci.lb <- res_corrected$g.co - qnorm(0.975)*res_corrected$se.co
    ci.ub <- res_corrected$g.co + qnorm(0.975)*res_corrected$se.co
    CI <- cbind(ci.lb, ci.ub)
    Boot_res<-cbind(est = res_corrected$g.co, se = res_corrected$se.co, CI = CI)
    tree<-Metatree
    tree$boot_res<-Boot_res
    tree$boot_tau2 <- new_tau
    tree$boot_out <- Boot_out
    #res <- list(Boot_res=Boot_res, tau2 = new_tau, Boot_out = Boot_out)
    
  }
  return(tree) #res)
  
  })
}
