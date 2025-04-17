expit <- function(x) (tanh(x/2)+1)/2 
Sright <- function(cpt, xk, a) expit(a*(xk-cpt))
Iright <- function(cpt, xk) {1*(xk > cpt)}
comp_Q_within <- function(y, vi, i.node){
  sum(y^2*i.node/vi) - (sum(y*i.node/vi))^2/sum(i.node/vi)
}
comp_C <- function(vi, i.node){
  sum(i.node/vi^2)/sum(i.node/vi)
}
comp_tau2 <- function(vi, sumQ, sumC, df) {
  max((sumQ - df)/(sum(1/vi) - sumC),0)
}

comp_Qb_tau2 <- function(y, vi, node) {
  uni.node <- unique(node)
  sumQC <- compute_left_(y, vi, node, uni.node)
  tau2 <- comp_tau2(vi, sumQC[1], sumQC[2], df = length(y)- length(uni.node))
  vi.star <- vi + tau2
  sumQ.star <- compute_left_(y, vi.star, node, uni.node)[1]
  Qb.star <- sum(y^2/vi.star) - (sum(y/vi.star))^2/sum(1/vi.star) - sumQ.star
  c(Qb.star, tau2)
}

# test
# rebetQ<- function(yi, vi, mods){
#   wts = 1/vi
#   wy = wts*yi
#   wy2 = wts * yi^2
#   Q <- tapply(wy2, mods, sum) - tapply(wy, mods, function(x) (sum(x))^2)/tapply(wts, mods, sum)
#   df <- tapply(wy, mods, length)-1
#   C <- tapply(wts, mods, sum) - tapply(wts, mods, function(x) sum(x^2))/ tapply(wts, mods, sum)
#   tau2 <- (sum(Q) - sum(df))/sum(C)
#   tau2 <- max(0, tau2)
#   wstar = 1/(vi+tau2)
#   wystar = wstar*yi
#   wy2star = wstar*yi^2
#   Qstar <- tapply(wy2star, mods, sum) - tapply(wystar, mods, function(x) (sum(x))^2)/tapply(wstar, mods, sum)
#   Qstar.total <- sum(wy2star) - (sum(wystar))^2/sum(wstar)
#   Qbet <- Qstar.total - sum(Qstar)
#   if (is.na(Qbet)) {
#     Qbet <- Inf
#   }
#   return(c(Qbet, tau2))
#   
# }
# y = dat$efk
# vi = dat$vark
# xk1 = dat$x1
# xk2 = dat$x2
# cpt1 = 0.5
# cpt2 = 0.5
# mods <- sample(1:17, length(y), replace = T)
# comp_Qb_tau2(y, vi, mods)
# rebetQ(y, vi, mods)
# system.time(for(i in 1:100) comp_Qb_tau2(y, vi, mods))
# system.time(for(i in 1:100) rebetQ(y, vi, mods))
