# library(TailRank)
# library(tidyverse)
# 
# plot(function(x) dnorm(x), -10, 10,
#      main = "Normal density")
# 
# ggplot(data.frame(x=c(-5, 5)), aes(x=x)) + 
#   stat_function(fun = function(x) dnorm(x) ) +  
#   geom_vline(xintercept = qnorm(c(0.025,0.975)), colour="red", linetype=2)
# 
# # #________ beta binomial function ______
# # 
# # N = 100
# # a = 30 
# # b = 70
# # 
# # ggplot(data.frame(x=c(0, N)), aes(x=x)) + 
# #   stat_function(fun = function(x) dbb(x, N, a, b) )
# # 
# # N = 1000
# # a = 300
# # b = 700
# # 
# # ggplot(data.frame(x=c(0, N)), aes(x=x)) + 
# #   stat_function(fun = function(x) dbb(x, N, a, b) )
# 
# #================= beta distribution function model for mutations ====================
# 
# N= 30 # number of trials (depth / total reads)
# k = 4 # number of successes (mutated reads)
# 
# alpha = k + 1
# beta = N - k + 1
# 
# ggplot(data.frame(x=c(0, 1)), aes(x=x)) + 
#   stat_function(fun = function(x) pbeta(x, alpha, beta) ) + geom_vline(xintercept = k/N, colour="blue", linetype=2)
# 
# CI <- qbeta(c(0.025,0.975), alpha, beta)
# 
# ggplot(data.frame(x=c(0, 1)), aes(x=x)) + 
#   stat_function(fun = function(x) dbeta(x, alpha, beta) ) + 
#   geom_vline(xintercept = k/N, colour="blue", linetype=2) +
#   geom_vline(xintercept = CI, colour="red", linetype=2)
# 
# 
# qnorm(0.954)
# 
# 2 * (34.1 + 13.6)
# 
# pnorm(2)
# dnorm(2)
# qnorm(0.9772)


######################################################

# require(TailRank)
require(tidyverse)

# depth <- 30
# alt <- 25

mut_VAF_CI <- Vectorize(function(depth, alt, plot=FALSE){
  
  alpha = alt + 1
  beta = depth - alt + 1
  
  
  CI <- qbeta(c(0.025,0.975), alpha, beta)
  
  VAF_hat <- alt/depth
  
  CI_low <- CI[1]
  CI_up <- CI[2]
  
  if(plot==TRUE){
    
    gg <- ggplot(data.frame(x=c(0, 1)), aes(x=x)) + 
      stat_function(fun = function(x) dbeta(x, alpha, beta) ) + 
      geom_vline(xintercept = VAF_hat, colour="blue", linetype=2) +
      geom_vline(xintercept = CI, colour="red", linetype=2)
    print(gg)
  }
  
  out <- c(CI_low= CI_low, VAF_hat= VAF_hat, CI_up= CI_up)
  out
  
})


# mut_VAF_CI(59, 10)
# 
# df <- data.frame(DEP=c(50,67,43,100,10), ALT= c(3,56,34,10,4))
# 
# mut_VAF_CI(df$DEP, df$ALT)
