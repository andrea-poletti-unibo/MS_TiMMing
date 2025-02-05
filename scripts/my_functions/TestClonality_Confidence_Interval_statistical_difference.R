
# =========== CHECK FOR STATISTICAL DIFFERENCE BETWEEN pl1 (score and CI) and pl2 (score and CI) ===============

# https://statisticsbyjim.com/hypothesis-testing/confidence-intervals-compare-means/
# https://psycnet.apa.org/fulltext/2005-01817-003.pdf?auth_token=855585f5601dde59579c2903b21824e130c0ff13

# 1) check the Cumming & Finch (2005) assumption: CI-1 and CI-2 not differing by more than a factor of 2 
# 2) compute "proportion overlap" 


# define data
pl1_score <- 0.75
pl2_score <- 0.46

CI_1 <- 0.10
CI_2 <- 0.20


TestClonality <- function(pl1_score, pl2_score, CI_1, CI_2, outplot=FALSE){
  
  require(tidyverse)
  
  # 1) check assumption
  if(CI_1/CI_2 > 2 | CI_1/CI_2 < 0.5 ) {

        # if no assumption met, use instead the Goldstein and Healy (1995) rule: for barely non-overlapping intervals 
    # to represent a 95% significant difference between two means, use an 83% confidence interval 
    # of the mean for each group and check if they overlap.
    
    # compute the 83% Confidence Intervals
    CI_1s <- CI_1 * (83/95)
    CI_2s <- CI_2 * (83/95)
    
    # check the overlap between the 83% CIs
    overlap83 <- if(pl1_score >= pl2_score) {
      p1_low <- pl1_score - CI_1s 
      p2_upp <- pl2_score + CI_2s
      p1_low < p2_upp
    } else if (pl2_score > pl1_score) {
      p2_low <- pl2_score - CI_2s 
      p1_upp <- pl1_score + CI_1s
      p2_low < p1_upp
    } else {"error"}
    
    alternative_significative <- !overlap83
  } else {
    alternative_significative <- FALSE
  }
  
  
  
  
  t_df <- data.frame(Players = c("player1" ,"player2"),
                     clonality_score = c(pl1_score ,pl2_score),
                     CI= c(CI_1, CI_2), 
                     CI_83 = c(CI_1 * (83/95) ,CI_2 * (83/95)) %>% round(3)
  )
  
  # 2) compute the "proportion overlap"

  average_margin_error <- t_df$CI %>% mean
  
  overlap <- if(pl1_score >= pl2_score) { (pl2_score + CI_2) - (pl1_score - CI_1) 
  } else if (pl2_score > pl1_score) { (pl1_score + CI_1) - (pl2_score - CI_2)}
  
  proportion_overlap <- (overlap/average_margin_error) %>% round(3)
  
  
  # 3) - optional - plot the comparison
  
  if (outplot==TRUE) {
    gg <- t_df %>% ggplot(aes(Players, clonality_score)) + 
      geom_point() + 
      geom_errorbar(aes(ymin=clonality_score-CI, ymax=clonality_score+CI), width=.1) +
      geom_errorbar(aes(ymin=clonality_score-CI_83, ymax=clonality_score+CI_83), width=.05, colour="red", alpha=0.5) +
      
      # P2 > P1 score
      {if (pl2_score > pl1_score){ 
        ggplot2::geom_segment(aes(x = 1, xend = 2, 
                                  y = c(clonality_score[1] + CI[1], clonality_score[2] - CI[2]), 
                                  yend = c(clonality_score[1] + CI[1], clonality_score[2] - CI[2])),
                              linetype=2) }  } +
      {if (pl2_score > pl1_score & proportion_overlap>0){
        ggplot2::geom_rect( aes(xmin = 1, xmax = 2, 
                                ymin = clonality_score[2] - CI[2], 
                                ymax = clonality_score[1] + CI[1]), 
                            alpha=0.1, fill= ifelse(proportion_overlap<0.5,"blue","red")) } } +
      
      
      # P1 > P2 score
      {if (pl1_score > pl2_score){ 
        ggplot2::geom_segment(aes(x = 1, xend = 2, 
                                  y = c(clonality_score[2] + CI[2], clonality_score[1] - CI[1]), 
                                  yend = c(clonality_score[2] + CI[2], clonality_score[1] - CI[1])),
                              linetype=2)}  } + 
      geom_hline(yintercept = c(0,1)) +
      
      {if (pl1_score > pl2_score & proportion_overlap>0){
        ggplot2::geom_rect( aes(xmin = 1, xmax = 2, 
                                ymin = clonality_score[1] - CI[1], 
                                ymax = clonality_score[2] + CI[2]), 
                            alpha=0.1, fill=ifelse(proportion_overlap<0.5,"blue","red")) } } +
      
      annotate("text", x=1.5, y= max(mean(pl1_score, pl2_score)-0.40, 0.10), 
               label= paste0("proportion overlap = ", proportion_overlap,"\n", ifelse(proportion_overlap<0.5,"p<0.05","n.s")))  
    
    print(gg)
    
    
  }
  
  
  res <- ifelse(alternative_significative==T|proportion_overlap<=0.5, "significative p<0.05", "not significative" )
  
  out <- list(test_result=res,
              table=t_df,
              alterntive=alternative_significative,
              proportion_overlap=proportion_overlap,
              overlap=overlap,
              average_margin_error=average_margin_error)
  
  return(out)
  
}


TestClonality(pl1_score = 0.75, 
              pl2_score = 0.27, 
              CI_1 = 0.16, 
              CI_2 = 0.10, 
              outplot=T)
# 
# TestClonality(pl1_score = 0.8, 
#               pl2_score = 0.4, 
#               CI_1 = 0.28, 
#               CI_2 = 0.2, 
#               outplot=T)
# 
# TestClonality(pl1_score = 0.2, 
#               pl2_score = 0.4, 
#               CI_1 = 0.1, 
#               CI_2 = 0.1,
#               outplot=T)
# 
# # alt
# TestClonality(pl1_score = 0.2, 
#               pl2_score = 0.6, 
#               CI_1 = 0.1, 
#               CI_2 = 0.25, 
#               outplot=T)
# 
# TestClonality(pl1_score = 0.4, 
#               pl2_score = 0.6, 
#               CI_1 = 0.2, 
#               CI_2 = 0.1,
#               outplot=T)





#=========================== summarized T-test approach =============================

computeSDfromCI <- function(CI, N){
  
  require(tidyverse)
  
  # factor = 3.920 if N > 60 (assume normality and use normal-distribution: 1.96 * 2)
  # factor = 4.128 if N < 60 (use t-distibution: 2.064 * 2 )
  
  factor <- ifelse(N>60, 3.920, 4.128 )

  SD <- sqrt(N)*(CI*2)/factor
  SD
}

# computeSDfromCI(0.1, 100)


# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

# res <- t.test2(0.36,0.5,0.5,0.5,100,100)
# res
# 
# res2 <- BSDA::tsum.test(mean.x = 0.36, s.x = 0.5, n.x = 100,  
#                 mean.y = 0.5, s.y = 0.5, n.y = 100)
# 
# res2$p.value
# 
# SD_1 <- computeSDfromCI(CI_1, 100)
# SD_2 <- computeSDfromCI(CI_2, 100)




######################### aggregate FUNCTION ############################

# N_1 <- 100
# N_2 <- 80
# 
# testClonality <- function(pl1_score, pl2_score, 
#                           CI_1, CI_2, 
#                           N_1, N_2) {
#   
#   
#   res_viz <- CI_viz_test(pl1_score, pl2_score, CI_1, CI_2)
#   
#   
#   SD_1 <- computeSDfromCI(CI_1, N_1)
#   SD_2 <- computeSDfromCI(CI_2, N_2)
#   
#   res_t <- res2 <- BSDA::tsum.test(mean.x = pl1_score, s.x = SD_1, n.x = N_1,  
#                                    mean.y = pl2_score, s.y = SD_2, n.y = N_2)
#   
#   res_t
#   
#   out <- list(viz=res_viz, t.test= res_t)
#   out
# }
# 
# 
# RES <- testClonality(0.3, 0.6, 0.1, 0.1, 100, 100)
# 
# RES$t.test$p.value
# RES$viz$test_result
# 
# 
# 
# RES <- testClonality(pl1_score = 0.3, pl2_score = 0.4, CI_1 = 0.1, CI_2 = 0.1, N_1 = 40, N_2 = 60)
# RES
# RES$t.test$p.value
# 
# RES <- testClonality(0.8, 0.6, 0.15, 0.10, 100, 100)
# RES
# RES$t.test$p.value



##################### for presentation ######################

# t_df <- data.frame(Players = c("Alteration A" ,"Alteration B"),
#                    clonality_score = c(0.75, 0.46),
#                    CI= c(0.10, 0.20)
# )

# gg <- t_df %>% ggplot(aes(Players, clonality_score)) +
#   geom_errorbar(
#     aes(
#       colour = Players,
#       ymin = clonality_score - CI, 
#       ymax = clonality_score + CI), size =1, width = .2) +
#   geom_point(size=5, aes(fill=Players), shape=21, stroke=2) +
#   geom_hline(yintercept = c(0, 1)) +
#   geom_hline(yintercept = c(0.2,0.4,0.6,0.8), linetype = 3) +
#   annotate("text",
#     x = 1.5, y = 1.1,
#     label = paste0("t-test p-value", "\n", "p<0.05")
#   ) +
#   theme_minimal() +
#   # disable y grid lines
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.major.y = element_blank(), 
#     panel.grid.minor.y = element_blank(),
#     # label size 
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     # axis label size
#     axis.title.y = element_text(size = 14),
#     # disable legend
#     legend.position = "none"
#     ) +
#   # y axix labels ticks frequency
#   scale_y_continuous( breaks = seq(0, 1, by = 0.2)) +
#   # labels 
#   labs(
#     x = "PLAYERS",
#     y = "Clonality level"
#   )
  


# print(gg)
