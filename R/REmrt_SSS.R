#' A function to fit the tree with look-ahead option
#'
#' @param mf the data.frame to grow the tree
#' @param maxL the maximum number of splits
#' @param minbucket the minimum number of the studies in a terminal node
#' @param minsplit the minimal number of studies in a parent node to be split
#' @param cp the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @param lookahead an argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @param alpha.endcut parameter used in the splitting algorithm to avoid the endcut preference problem.
#' @param a parameter used in the sss to determine the slope of the logistic function that replaces the indicator function.
#' @param multi.start boolean indicating whether multiple starts must be used
#' @param n.starts number of multiple starts
#' @return a list including a tree, the split points, the data, and the nodes after each split
#' @keywords internal
#' @importFrom stats terms model.response 
REmrt_SSS <- function(mf, maxL, minbucket, minsplit, cp, lookahead,
                      alpha.endcut = 0.02, a = 50,
                      multi.start=T, n.starts=3){
  #===================  Error message  ======================#
  if (alpha.endcut < 0 | alpha.endcut > 1) {
    stop("alpha.endcut should range between 0 and 1")
  }
  if (!is.logical(lookahead)) {
    stop("lookahead should be TRUE or FALSE")
  }
  if (lookahead & maxL < 2) {
    stop("the maximum number of splits should be larger than two when using lookahead")
  }
  if (nrow(mf) < minsplit) {stop("number of studies smaller than minsplit")}
  
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  original_data<-mods
  binary<-isBinary(mods) #Boolean vector. If TRUE, column is binary.
  num_vars<-which(sapply(mods, is.numeric)==TRUE)
  
  mods<-as.data.frame(sapply(mods,as.numeric))
  mods <- scale(mods[,sapply(mods, is.numeric)])
  
  # mus <- rep(NA, NCOL(mods))
  # sigmas <- rep(NA, NCOL(mods))
  # if(length(num_vars)>0){# STANDARIZE THE NUMERIC MODERATORS
  #   for (i in num_vars) {
  #     tempmods <- scale(mods[,i]) 
  #     mus[i] <- attr(tempmods, "scaled:center")
  #     sigmas[i] <- attr(tempmods, "scaled:scale")
  #     mods[,i] <-tempmods
  #   }
  # }
  mus <- attr(mods, "scaled:center")
  sigmas <- attr(mods, "scaled:scale")
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(rep(1, nrow(mf)))
  res.Qb <- 0
  res.tau2 <- (sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)   - length(y)+1)/(sum(1/vi)-sum(1/vi^2)/sum(1/vi)) #VERIFIED
  res.split <- NA
  res.mod <- NA
  res.pleaf <- NA
  delta.Q <- Inf
  Dev <- -Inf 
  
  
  if (lookahead) {
    i <- 1
    cnode <- nodemark[,i]    
	for (k1 in 1:nmod){
      xk1 <- mods[, k1]
####### START XINRU MODIFICATION 3 ############
      if (length(unique(xk1)) < 2) next #if number of categories less than 2 we go to the next mod
      for (k2 in 1:nmod){
        xk2 <- mods[, k2]
        if (length(unique(xk2)) < 2) next #if number of categories less than 2 we go to the next mod
        if ((length(unique(xk1)) == 2) & (length(unique(xk2)) == 2)) {
          ####### START XINRU MODIFICATION 1 ############
		  ## both binary
          if (k1 == k2) next
          first.splits <- lapply(list(mods[, k1]), function(x) make_first_split(x, 
                                                                    y, minbucket))
          second.splits <- lapply(list(k2), function(x) find_second_split(mods[, 
                                                                             x], first.splits, y, vi, minbucket, minsplit))
          tempQs <- sapply(1:length(second.splits), function(x) second.splits[[x]]$Q)
          x2.inx <- which.max(tempQs)
          x1.inx <- second.splits[[x2.inx]]$split1[1]
          tempQ2 <- second.splits[[x2.inx]]$Q 
          if (tempQ2 > Dev) { #old: tempQ2 < Dev
            Dev <- tempQ2
            #old:cstar1 <- second.splits[[x2.inx]]$split1[2]
			      cstar1 <- max(xk1)
			      cstar2 <- max(xk2)				
            node1 <- first.splits[[x1.inx]]$childNodes[, second.splits[[x2.inx]]$split1[3]]
            names(node1) <- NULL
            pleaf.inx <- node1 == second.splits[[x2.inx]]$split2$pleaf
            node2 <- node1
            #cstar2 <- names(second.splits[[x2.inx]]$split2$rank[second.splits[[x2.inx]]$split2$rank <= 
            #                                                      second.splits[[x2.inx]]$split2$cstar[1]])
            node2[pleaf.inx] <- ifelse(mods[pleaf.inx, x2.inx] %in% 
                                         cstar2, 4, 5)
			      sigma1 <- sigmas[k1]; mu1 <- mus[k1]
            sigma2 <- sigmas[k2]; mu2 <- mus[k2]								  
            #old: sigma2 <- sigmas[mod2]; mu2 <- mus[mod2]
            cpt1 <- cstar1*sigma1 + mu1
            cpt2 <- cstar2*sigma2 + mu2
            # msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
            # msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
            if(k1%in%num_vars){#old: mod1%in%num_vars){
              msplit1 <- paste(mods.names[k1], "<", cpt1, collapse = " ")
            }else{ #categorical
              #msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
              cstar1_cat <- levels(original_data[,k1])[1:(ceiling(cpt1)-1)]
              msplit1 <- paste(mods.names[k1], "=", paste(cstar1_cat, 
                                                            collapse = "/"), collapse = " ")
            }
            if(k2%in%num_vars){
              msplit2 <- paste(mods.names[k2], "<", cpt2, collapse = " ")
            }else{
              #msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
              cstar2_cat <- levels(original_data[,k2])[1:(ceiling(cpt2)-1)]
              msplit2 <- paste(mods.names[k2], "=", paste(cstar2_cat, 
                                                            collapse = "/"), collapse = " ")
            }
            
            pleaf <- second.splits[[x2.inx]]$split2$pleaf
            mod1.name <- mods.names[k1]
            mod2.name <- mods.names[k2]
          }#Added bracket
          # res.mod <- c(NA, mods.names[c(x1.inx, x2.inx)])
          # res.pleaf <- c(NA, 1, second.splits[[x2.inx]]$split2$pleaf)
          ####### END XINRU MODIFICATION 1 ############  
          
        } else if ((length(unique(xk1)) == 2) & (length(unique(xk2)) > 2)) {
          ## first binary second numerical
		
## start: first split 
        ####### START XINRU MODIFICATION 2 ############	
          
		#print("re.cutoff_cpp number1")      
        pleaf.inx1 <- cnode == 1
        temp1 <- re.cutoff_cpp(y, vi, xk1, pleaf.inx1, 
                                              cnode, minbucket)
        cstar1 <- temp1[1]
        cnode1 <- ifelse(xk1 <= cstar1, 
                           2 * i, 2 * i + 1)
        ## end: first split
        ## start: second split using SSS
        for (pl in 2:3){
          pleaf.inx2 <- cnode1 == pl
          cstar2 <- re_cutoff_SSS(y, vi, xk2[pleaf.inx2], cnode1, as.numeric(pl),
                                  a = 50,
                                  alpha.endcut = 0.02,
                                  multi.start=T, n.starts=3)
          
          temp.node2 <- cnode1
          temp.node2[pleaf.inx2] <- ifelse( xk2[pleaf.inx2] < cstar2, 2*(i+1), 2*(i+1)+1)
          temp2 <- comp_Qb_tau2(y, vi, temp.node2)
          if (is.na(temp2[1])) {
            Dev.new <- -Inf
          } else {
            Dev.new <- temp2[1]


          }
          if (Dev.new > Dev) {
            Dev <- Dev.new
            pleaf <- pl
            sigma1 <- sigmas[k1]; mu1 <- mus[k1] 
            sigma2 <- sigmas[k2]; mu2 <- mus[k2]
            cpt1 <- cstar1*sigma1 + mu1
            cpt2 <- cstar2*sigma2 + mu2
            #msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
            if(k1%in%num_vars){
              msplit1 <- paste(mods.names[k1], "<", cpt1, collapse = " ")
            }else{
              cstar_cat <- levels(original_data[,k1])[1:(ceiling(cpt1)-1)]
              msplit1 <- paste(mods.names[k1], "=", paste(cstar_cat, 
                                                          collapse = "/"), collapse = " ")
            }
            if(k2%in%num_vars){							 
              msplit2 <- paste(mods.names[k2], "<", cpt2, collapse = " ")
            }else{
              cstar_cat <- levels(original_data[,k2])[1:(ceiling(cpt2)-1)]
              msplit2 <- paste(mods.names[k2], "=", paste(cstar_cat, 
                                                          collapse = "/"), collapse = " ")										 
            }
            
           
            node1 <- ifelse(mods[, k1] <= cstar1, 2, 3) #Changed from < to <=
            node2 <- node1
            pleaf.inx <- node1 == pleaf
            node2[pleaf.inx] <- ifelse(mods[pleaf.inx, k2] < cstar2, 4, 5)
            mod1.name <- mods.names[k1]
            mod2.name <- mods.names[k2]

          }
          ## end: second split SSS
          ####### END XINRU MODIFICATION 2 ############
        }
      
		  
        } else if ((length(unique(xk1)) > 2) & (length(unique(xk2)) == 2)) {
          ## first numerical second binary
          #---first split sss
          cstar1 <- re_cutoff_SSS(y, vi, xk1, cnode, as.numeric(1),
                                               a = 50,
                                               alpha.endcut = 0.02,
                                               multi.start=T, n.starts=3)
          cnode1 <- cnode
          cnode1 <- ifelse( xk1 < cstar1, 2*i, 2*i+1)
          temp1 <- comp_Qb_tau2(y, vi, cnode1)
          #--- second split GS
          for (pl in 2:3) {
            pleaf.inx2 <- cnode1 == pl
            temp2 <- re.cutoff_cpp(y, vi, xk2[pleaf.inx2], pleaf.inx2, 
                                                cnode1, minbucket)
            if (is.null(temp2)) {
              Dev.new <- -Inf
            }
            else {
              Dev.new <- temp2[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp2[2]
              cstar2 <- temp2[1]
              pleaf <- pl
              sigma1 <- sigmas[k1]; mu1 <- mus[k1] 
              sigma2 <- sigmas[k2]; mu2 <- mus[k2]
              cpt1 <- cstar1*sigma1 + mu1
              cpt2 <- cstar2*sigma2 + mu2
              #msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
              if(k1%in%num_vars){
                msplit1 <- paste(mods.names[k1], "<", cpt1, collapse = " ")
              }else{
                cstar_cat <- levels(original_data[,k1])[1:(ceiling(cpt1)-1)]
                msplit1 <- paste(mods.names[k1], "=", paste(cstar_cat, 
                                                            collapse = "/"), collapse = " ")
              }
              if(k2%in%num_vars){
                msplit2 <- paste(mods.names[k2], "<", cpt2, collapse = " ")
              }else{
                cstar_cat <- levels(original_data[,k2])[1:(ceiling(cpt2)-1)]
                msplit2 <- paste(mods.names[k2], "=", paste(cstar_cat, 
                                                            collapse = "/"), collapse = " ")
              }
              
              
              node1 <- ifelse(mods[, k1] < cstar1, 2, 3)
              node2 <- node1
              pleaf.inx <- node1 == pleaf
              node2[pleaf.inx] <- ifelse(mods[pleaf.inx, k2] <= cstar2, 4, 5) #Changed from < to <=
              mod1.name <- mods.names[k1]
              mod2.name <- mods.names[k2]
            }
          } 
          ############# begin xinru modification 3############
        } else {
          ## both numerical -> SSS
          if (k1 == k2) {
            
            
            temp <- find_two_cutoffs_on_one_var(y, vi, xk1, a, assump = "re",
                                                multi.start = T, n.starts=3)
            
          } else {
            xk2 <- mods[, k2]
            temp <- find_two_cutoffs(y, vi, xk1, xk2, a, assump = "re",
                                     multi.start = T)				   
          }
		      cstar1_tmp <- temp[3]
          cstar2_tmp <- temp[4]
          tcnode1 <- ifelse(xk1 < cstar1_tmp, 2, 3)
          tpleaf.inx <- tcnode1 == as.integer(temp[1])
          tcnode2 <- tcnode1
          tcnode2[tpleaf.inx] <- ifelse(xk2[tpleaf.inx] < cstar2_tmp, 4, 5)
          tempQ <- comp_Qb_tau2(y, vi, tcnode2)
          if (is.null(tempQ)) {
            Dev.new <- -Inf
          }
          else {
            Dev.new <- tempQ[1] #Change Juan. Old: tempQ[2]
          }
          
          
          if(Dev.new > Dev) {
            Dev <- Dev.new														   
            pleaf <- as.integer(temp[1])
            cstar1 <- cstar1_tmp #either cstar1 or cstar2 should be 0.5 for dichotomous vars
            cstar2 <- cstar2_tmp
            mod1 <- k1 
            mod2 <- k2
            sigma1 <- sigmas[mod1]; mu1 <- mus[mod1]
            sigma2 <- sigmas[mod2]; mu2 <- mus[mod2]
            cpt1 <- cstar1*sigma1 + mu1
            cpt2 <- cstar2*sigma2 + mu2
            # msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
            # msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
            if(k1%in%num_vars){
              msplit1 <- paste(mods.names[k1], "<", cpt1, collapse = " ")
            }else{
              cstar_cat <- levels(original_data[,k1])[1:(ceiling(cpt1)-1)]
              msplit1 <- paste(mods.names[k1], "=", paste(cstar_cat, 
                                                        collapse = "/"), collapse = " ")
            }
            
            if(k2%in%num_vars){
              msplit2 <- paste(mods.names[k2], "<", cpt2, collapse = " ")
            }else{
              cstar_cat <- levels(original_data[,k2])[1:(ceiling(cpt2)-1)]
              msplit2 <- paste(mods.names[k2], "=", paste(cstar_cat, 
                                                          collapse = "/"), collapse = " ")
            }
            node1 <- tcnode1
            node2 <- tcnode2
            # pleaf.inx <- node1 == pleaf
            # node2[pleaf.inx] <- ifelse(mods[pleaf.inx, mod2] < cstar2, 4, 5)
            mod1.name <- mods.names[k1]
            mod2.name <- mods.names[k2]
            ## needed output: node1, node2, cpt1, cpt2, msplit1, msplit2, mod1, mod2, pleaf1=1, pleaf2=pleaf
          }
        }
        
        

      }
    }
    if (is.infinite(Dev)) {
      stop("no possible splits for look-ahead")
    } else {
      # sigma1 <- sigmas[mod1]; mu1 <- mus[mod1]
      # sigma2 <- sigmas[mod2]; mu2 <- mus[mod2]
      # cpt1 <- cstar1*sigma1 + mu1
      # cpt2 <- cstar2*sigma2 + mu2
      # msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
      # msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
      # # if(mod1%in%num_vars){
      # #   msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
      # # }else{ #categorical
      # #   #msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
      # #   cstar1_cat <- levels(original_data[,mod1])[1:(ceiling(cpt1)-1)]
      # #    msplit1 <- paste(mods.names[mod1], "=", paste(cstar1_cat, 
      # #                                                  collapse = "/"), collapse = " ")
      # # }
      # # if(mod2%in%num_vars){
      # #   msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
      # # }else{
      # #   #msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
      # #   cstar2_cat <- levels(original_data[,mod2])[1:(ceiling(cpt2)-1)]
      # #    msplit2 <- paste(mods.names[mod2], "=", paste(cstar2_cat, 
      # #                                                  collapse = "/"), collapse = " ")
      # # }

      # # msplit1 <- paste(mods.names[mod1], "<", round(cpt1,2), collapse = " ")
      # # msplit2 <- paste(mods.names[mod2], "<", round(cpt2,2), collapse = " ")
      # node1 <- ifelse(mods[, mod1] < cstar1, 2, 3)
      # node2 <- node1
      # pleaf.inx <- node1 == pleaf
      # node2[pleaf.inx] <- ifelse(mods[pleaf.inx, mod2] < cstar2, 4, 5)
      tempQ1 <- comp_Qb_tau2(y, vi, node1)
      tempQ2 <- comp_Qb_tau2(y, vi, node2)
      res.Qb <- c(res.Qb, tempQ1[1], tempQ2[1])
      res.tau2 <- c(res.tau2,  tempQ1[2], tempQ2[2])
      res.split <- c(res.split, msplit1, msplit2)
      res.mod <- c(res.mod, mod1.name, mod2.name ) #res.mod <- c(res.mod, mods.names[mod1], mods.names[mod2] )
      res.pleaf <- c(res.pleaf, 1, pleaf)
      nodemark <- cbind(nodemark, node1, node2)
      cpt[[1]] <- cpt1
      cpt[[2]] <- cpt2
      delta.Q <- res.Qb[3] - res.Qb[2]
      i <- 2
    }
####### END XINRU MODIFICATION 3 ############ 		
#End lookahead									  
  } else{
#here the non-lookahead starts									   
    for (i in 1) {
	  cnode <- nodemark[,i]
      
      Dev<- -Inf
      TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
      cnode <- nodemark[ ,i]
      len.node <- tapply(vi, cnode, length)
      nodes <- names(len.node) [len.node >= minsplit]
      for (pl in nodes) {
        pleaf.inx <- cnode == pl
        for (k in 1:nmod) {
          xk <- mods[pleaf.inx, k]
          c.splits <- unique(xk)
          if (length(c.splits) < 2) next
####### START XINRU MODIFICATION 1 ############
          if (length(c.splits) == 2) { # binary moderator
            # if (is.numeric(xk)) {
              #print("re.cutoff_cpp number2")
              temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, 
                                    cnode, minbucket)
              if (is.null(temp)) {
                Dev.new <- -Inf																				  
              }
              else {
                Dev.new <- temp[2]
              }
              if (Dev.new > Dev) {
                Dev <- Dev.new
              sigmak <- sigmas[k]; muk <- mus[k]
              c.star <- temp[1]
              cptk <- c.star*sigmak + muk
              #msplit <- paste(mods.names[k], "<=", c.star, collapse = " ")
              if(k%in%num_vars){
                msplit <- paste(mods.names[k], "<=", c.star, collapse = " ") #cptk or c.star? CHECK
              }else{
                cstar_cat <- levels(original_data[,k])[1:(ceiling(cptk)-1)]
                msplit <- paste(mods.names[k], "=", paste(cstar_cat, 
                                                          collapse = "/"), collapse = " ")
              }
              TQb = Dev.new
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse(xk <= c.star, 
                                            2 * i, 2 * i + 1)
              }
          } else if ((length(c.splits) > 2) & is.numeric(xk)) {
            # NUMERIC VARIABLE
            cstar <- re_cutoff_SSS(y, vi, xk, cnode, as.numeric(pl),
                                                a = 50,
                                                alpha.endcut = 0.02,
                                                multi.start=T, n.starts=3)
            
							   
									
									   
            temp.node <- cnode
              temp.node[pleaf.inx] <- ifelse( xk < cstar, 2*i, 2*i+1)
              temp <- comp_Qb_tau2(y, vi, temp.node)
              if (is.na(temp[1])) {
                Dev.new <- -Inf
              } else {
                Dev.new <- temp[1]
              }
              if (Dev.new > Dev) {
                Dev <- Dev.new
                sigmak <- sigmas[k]; muk <- mus[k]
                cptk <- cstar*sigmak + muk
                #msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
                if(k%in%num_vars){
                  msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
                }else{
                  cstar_cat <- levels(original_data[,k])[1:(ceiling(cptk)-1)]
                  msplit <- paste(mods.names[k], "=", paste(cstar_cat, 
                                                            collapse = "/"), collapse = " ")
                }
                #msplit <- paste(mods.names[k], "<", round(cptk,2), collapse = " ")
                TQb = Dev.new 
                Ttau2 = temp[2]
                Tsplit = msplit
                Tmod = mods.names[k]
                Tpleaf = as.numeric(pl)
                new.node <- temp.node
              }
            } else {
              stop("SSS for nominal variable is not finished yet")
            xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
            xk.ordinal <- xk.rank[as.character(xk)]
            temp <- re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
			   
            } else {
																									 
															 
													 
											 
																				   
								  
              Dev.new <- temp[2]
					  
								  
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- names(xk.rank[xk.rank < temp[1]])
              msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
			   
            }
          }
        }
      }
											 
      if (is.null(TQb)) {
        delta.Q <- -Inf
      } else {
        delta.Q <- abs(TQb - res.Qb[i])
      }
      
      nodemark <- cbind(nodemark, new.node)
      res.Qb <- c(res.Qb, TQb)
      res.tau2 <- c(res.tau2, Ttau2)
      res.split <- c(res.split, Tsplit)
      res.mod <- c(res.mod,Tmod)
      res.pleaf <- c(res.pleaf, Tpleaf)
      cpt[[i]] <- cptk}
  }#End no lookahead
  
  
  
  
  #-------          Continue with greedy procedure          -------#
  while(delta.Q >= cp & i < maxL) {
    i <- i+1
    TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    for (pl in nodes) {
      pleaf.inx <- cnode == pl
      for (k in 1:nmod) {
        xk <- mods[pleaf.inx, k]
        c.splits <- unique(xk)
        if (length(c.splits) < 2) next
####### START XINRU MODIFICATION 2############
        if (length(c.splits) == 2) { # binary moderator
          if (is.numeric(xk)) {
            #print("re.cutoff_cpp number4")

            temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx,
                                  cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            }
            else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- Dev.new
              sigmak <- sigmas[k]; muk <- mus[k]
              c.star <- temp[1]
              cptk <- c.star*sigmak + muk
              #msplit <- paste(mods.names[k], "<=", c.star, collapse = " ")
              if(k%in%num_vars){
                msplit <- paste(mods.names[k], "<=", c.star, collapse = " ") #cptk or c.star? CHECK
              }else{

                cstar_cat <- levels(original_data[,k])[1:(ceiling(cptk)-1)] #cptk or c.star? CHECK
                msplit <- paste(mods.names[k], "=", paste(cstar_cat,
                                                          collapse = "/"), collapse = " ")
              }
              TQb = Dev.new
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse(xk <= c.star,
                                            2 * i, 2 * i + 1)
            }
          } else {
            stop("GS binary: SSS for nominal variable is not finished yet")
          } # binary non-numerical moderator
        } else { # non-binary moderator
          if (is.numeric(xk)) {
            # NUMERIC VARIABLE

            cstar <- re_cutoff_SSS(y, vi, xk, cnode, as.numeric(pl),
                                   a = 50,
                                   alpha.endcut = 0.02,
                                   multi.start=T, n.starts=3)
            # print(cstar)
            temp.node <- cnode
            temp.node[pleaf.inx] <- ifelse( xk < cstar, 2*i, 2*i+1)
            temp <- comp_Qb_tau2(y, vi, temp.node)
            if (is.na(temp[1])) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[1]
            }
            if (Dev.new > Dev) {
              Dev <- Dev.new
              sigmak <- sigmas[k]; muk <- mus[k]
              cptk <- cstar*sigmak + muk

              #msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
              if(k%in%num_vars){
                msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
              }else{
                cstar_cat <- levels(original_data[,k])[1:(ceiling(cptk)-1)]
                msplit <- paste(mods.names[k], "=", paste(cstar_cat,
                                                          collapse = "/"), collapse = " ")
              }

              TQb = Dev.new
              Ttau2 = temp[2]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- temp.node
            }
          } else {
            stop("GS nominal >3 cat: SSS for nominal variable is not finished yet")
            xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
            xk.ordinal <- xk.rank[as.character(xk)]
            #print("re.cutoff_cpp number5")
            temp <- re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- names(xk.rank[xk.rank < temp[1]])
              msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
            }
          }
        }
      }
      # if (is.null(TQb)) {
        # delta.Q <- -Inf
      # } else {
        # delta.Q <- abs(TQb - res.Qb[i])
      # }
    }
####### END XINRU MODIFICATION 2############	
    # nodemark <- cbind(nodemark, new.node)
    # res.Qb <- c(res.Qb, TQb)
    # res.tau2 <- c(res.tau2, Ttau2)
    # res.split <- c(res.split, Tsplit)
    # res.mod <- c(res.mod,Tmod)
    # res.pleaf <- c(res.pleaf, Tpleaf)
    # cpt[[i]] <- cptk
  # }
  # #-------          Continue with greedy procedure          -------#
  # while(delta.Q >= cp & i < maxL) {
    # i <- i+1
    # TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
    # cnode <- nodemark[ ,i]
    # len.node <- tapply(vi, cnode, length)
    # nodes <- names(len.node) [len.node >= minsplit]
    # for (pl in nodes) {
      # pleaf.inx <- cnode == pl
      # for (k in 1:nmod) {
        # xk <- mods[pleaf.inx, k]
        # c.splits <- unique(xk)
        # if (length(c.splits) < 2) next
        # if (is.numeric(xk)) {
          # # NUMERIC VARIABLE
          # cstar <- re_cutoff_SSS(y, vi, xk, cnode, as.numeric(pl),
                                 # a = 50,
                                 # alpha.endcut = 0.02,
                                 # multi.start=T, n.starts=3)

          # temp.node <- cnode
          # temp.node[pleaf.inx] <- ifelse( xk < cstar, 2*i, 2*i+1)
          # temp <- comp_Qb_tau2(y, vi, temp.node)
          # if (is.na(temp[1])) {
            # Dev.new <- -Inf
          # } else {
            # Dev.new <- temp[1]
          # }
          # if (Dev.new > Dev) {
            # Dev <- Dev.new
            # sigmak <- sigmas[k]; muk <- mus[k]
            # cptk <- cstar*sigmak + muk
            # msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
            # # if(k%in%num_vars){
            # #   msplit <- paste(mods.names[k], "<", cptk, collapse = " ")
            # # }else{
            # #   cstar_cat <- levels(original_data[,k])[1:(ceiling(cptk)-1)]
            # #   msplit <- paste(mods.names[k], "=", paste(cstar_cat, 
            # #                                             collapse = "/"), collapse = " ")
            # # }
            
            # #msplit <- paste(mods.names[k], "<", round(cptk,2), collapse = " ")
            # TQb = Dev.new
            # Ttau2 = temp[2]
            # Tsplit = msplit
            # Tmod = mods.names[k]
            # Tpleaf = as.numeric(pl)
            # new.node <- temp.node
          # }
        # } else {
          # stop("SSS for nominal variable is not finished yet")
          # xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
          # xk.ordinal <- xk.rank[as.character(xk)]
          # temp <- re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
          # if (is.null(temp)) {
            # Dev.new <- -Inf
          # } else {
            # Dev.new <- temp[2]
          # }
          # if (Dev.new > Dev) {
            # Dev <- temp[2]
            # c.star <- names(xk.rank[xk.rank < temp[1]])
            # msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
            # TQb = temp[2]
            # Ttau2 = temp[3]
            # Tsplit = msplit
            # Tmod = mods.names[k]
            # Tpleaf = as.numeric(pl)
            # new.node <- cnode
            # new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
          # }
        # }
      # }
    # }		 
    if (is.null(TQb)) {
      delta.Q <- -Inf
    } else {
      nodemark <- cbind(nodemark, new.node)
      res.Qb <- c(res.Qb, TQb)
      res.tau2 <- c(res.tau2, Ttau2)
      res.split <- c(res.split, Tsplit)
      res.mod <- c(res.mod,Tmod)
      res.pleaf <- c(res.pleaf, Tpleaf)
      cpt[[i]] <- cptk
      delta.Q <- abs(TQb - res.Qb[i])
    }

    #}
  }
  

  
  return(list(tree = data.frame(Qb = res.Qb, tau2 = res.tau2, split = res.split,
                         mod = res.mod, pleaf = res.pleaf, stringsAsFactors = FALSE),
       node.split = nodemark, cpt = cpt, data = mf))

}




isBinary<-function(df){
  #df is a dataframe
  f_df<-as.data.frame(apply(df,2,as.factor), stringsAsFactors = TRUE)
  Binary<-rep(FALSE, NCOL(df))
  
  for (i in 1:NCOL(df)) {
    
    if(is.factor(f_df[,i])){#if column is factor
      
      n_level<-length(levels(df[,i]))
      if(n_level==2){ #if factor and two categories
        Binary[i]<-TRUE
      }else{ #if factor and more than 2 categories
        Binary[i]<-FALSE
      }
      
    }else{#if is not factor
      
      Binary[i]<-FALSE
      
    }
  }
  
  return(Binary)
}

