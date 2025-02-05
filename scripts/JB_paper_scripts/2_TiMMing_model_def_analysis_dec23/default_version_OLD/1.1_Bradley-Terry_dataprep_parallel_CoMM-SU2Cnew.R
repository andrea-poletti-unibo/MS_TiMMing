library(data.table)
library(tidyverse)

#============== load quality prepared data ================


calls <- fread("results/JB_paper/calls/clonalityCalls_broad_focal_CoMM-SU2Cnew.txt")

purity_file <- fread("data/CoMMpass_acs/CoMMpass_ABSOLUTE_tables_969samp.txt")
purity <- purity_file %>% select(sample, purity)

CI95 <- fread("results/JB_paper/calls/allCalls_CN_CI_broad_focal_CoMM-SU2Cnew.txt")


dfCI <- apply(CI95, c(1,2), function(x) str_remove(x,".*\\|"))

dfCI_m <- dfCI %>% as.data.frame() %>% 
  reshape2::melt(id.vars="ID") %>% 
  setNames(c("sample","chrarm","CI95"))



# melt the data => 1 row per event observation 
melt <- melt(calls, id.vars= "sample")

# add purity calls
melt <- left_join(melt, purity %>% dplyr::select(sample, purity), by=c("sample"))

# compute call thresholds based on purity levels
melt$thresh <- 0.05/melt$purity

melt$call <- ifelse(melt$value > melt$thresh,1,0) 


# add 95% Confidence Intervals

melt$chrarm <- melt$variable %>% str_remove("amp_|del_")

melt <- left_join(melt, dfCI_m, by=c("sample","chrarm"))

melt$CI95 <- as.numeric(melt$CI95)


#====================== generate HD call for each sample ====================

melt$HD_chr <- ifelse(melt$variable %>% str_detect("amp_chr_[3579]|amp_chr_1[159]|amp_chr_21"), "yes","no")

# melt %>% group_by(variable) %>% summarise(A=unique(HD_chr)) %>% View

cohort <- unique(melt$sample)

melt_HD <- data.frame()

pt=cohort[1]

# CALL HyperDiploidy CLONALITY (median value if n HD arms >= 3)
for ( pt in cohort){
  
  # message(which(pt==cohort))
  
  p <- melt %>% filter(sample == pt)
  p_HD <- p %>% filter(HD_chr=="yes" & call == 1)
  
  if(nrow(p_HD)>3){
    
    HD_value <- median( p_HD$value)
    HD_CI <- median( p_HD$CI95 )
    
    # modify value on arms to avoid redundancy
    p$value[p$HD_chr=="yes" & p$call == 1] <- -999
    p$call[p$HD_chr=="yes" & p$call == 1] <- -999
    
    
  } else {  HD_value <- 0 }
  
  HD_player <- data.frame(sample=pt,
                          variable="HyperDiploidy",
                          value=HD_value,
                          purity=p$purity %>% unique, 
                          thresh=p$thresh %>% unique,
                          call= ifelse(HD_value > 0 , 1, 0),
                          chrarm=NA,
                          CI95= HD_CI,
                          HD_chr=NA)
  
  p <- rbind(p, HD_player)
  
  melt_HD <- rbind(melt_HD,p)
}


# save melted callset
write_tsv(melt_HD, "workfiles/callset_CoMM-SU2Cnew_HDcall.txt")


#================= create the players ====================

players <- c(melt_HD$variable %>% unique() %>% as.character()) 
players

#================= create the matches =========================

matches <- t(combn(players, 2))

matches.df <- data.frame(player1 = matches[,1], 
                         player2 = matches[,2])

matches.df$player1 <- as.character(matches.df$player1)
matches.df$player2 <- as.character(matches.df$player2)
matches.df %>% str

matches.df$win1 <- 0
matches.df$win2 <- 0

# delete nosense matches (same chrarm - amp vs del and mut)
same.check <- function(var1, var2){ 
  v1 <- var1 %>% str_extract("[0-9]+[pq]|Hyper|^mut.*") 
  v2 <- var2 %>% str_extract("[0-9]+[pq]|Hyper|^mut.*") 
  same <- v1 != v2
  same %>% table()
  same
}


senseIDX <- same.check(matches.df$player1, matches.df$player2)
senseIDX %>% table

matches.dfREV <- matches.df[!senseIDX,]

matches.df <- matches.df[senseIDX,]

matches.df$matchID <- 1:nrow(matches.df)



##################### statistical Clonal Difference (ClonDiff) method ##################

# load TestClonality function
source("scripts/my_functions/TestClonality_Confidence_Interval_statistical_difference.R")

cohort <- unique(melt_HD$sample)

#================ PARALLEL LOOP over samples > tournaments ====================

#______parallel version___

library(parallel)
library(doParallel)
library(foreach)
library(doSNOW)

N_cores <- parallel::detectCores() - 4

cl <- makeCluster(N_cores)
registerDoSNOW(cl)


iterations <- length(cohort)

# iterations <- 3
# progress bar ------------------------------------------------------------
library(progress)
pb <- progress_bar$new(
  format = "letter = :letter [:bar] :elapsed | eta: :eta",
  total = iterations,    # 100 
  width = 60)
progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar

# allowing progress bar to be used in foreach -----------------------------
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 

opts <- list(progress = progress)

i <- 1

parRes <- foreach(i=1:iterations, 
                  .packages = "tidyverse", 
                  .options.snow = opts) %dopar% {
                    
  pt <- cohort[i]
  tournament.pt <- melt_HD %>% filter(sample==pt) 
  matches_res <- matches.df
  
  matches_res$performed <- 0
  matches_res$sample <- pt
  
  m=1
  for(m in 1:nrow(matches_res)){
    
    if(m %% 500 ==0) message(m)
    
    mtc <- matches_res[m,1:2]
    pl1 <- mtc$player1
    pl2 <- mtc$player2
    
    pl1_stats <- tournament.pt %>% filter(variable == pl1)
    pl2_stats <- tournament.pt %>% filter(variable == pl2)
    
    #_____________ if both event are present: perform a match ________________
    
    if(pl1_stats$call == 1 & pl2_stats$call == 1) {
      
      # check if distributions are significantly different based on value and CI95% (method 1 - TestClonality)
      testclon_res <- TestClonality(pl1_score = pl1_stats$value,
                                    pl2_score = pl2_stats$value, 
                                    CI_1 = pl1_stats$CI95, 
                                    CI_2 = pl2_stats$CI95,
                                    outplot = F)
      
      # do actual match
      if(testclon_res$test_result == "significative p<0.05") { 
        
        # declare winner
        winner <- (which.max( c(pl1_stats$value, pl2_stats$value)  ))
        
        # declare points = proportional to the clonality difference (ClonDiff method)
        points <- round((abs(pl1_stats$value - pl2_stats$value) * 100) / 20)
        
        # assign points to winner and record 
        colnames(matches_res)
        if(winner==1 & points>0){
          matches_res[m,"win1"] <- matches_res[m,"win1"] + points
          matches_res[m,"performed"] <- 1
          
        } else if(winner==2 & points>0){
          matches_res[m,"win2"] <- matches_res[m,"win2"] + points 
          matches_res[m,"performed"] <- 1
        }
      }
    }
  }
  
  matches_res
}

stopCluster(cl)




#___aggregate all tables scores____

all <- Reduce(rbind, parRes)

expscorestab <- all %>% group_by(matchID) %>% 
  summarise(player1= unique(player1), 
            player2=unique(player2), 
            win1= sum(win1), 
            win2=sum(win2),
            total_matches_performed=sum(performed))


# save single tables with points ## WARNING too big! 150+ Mb
# write_tsv(all,"workfiles/Bus_matches_tables_ClonDiffPoints.txt")

# save aggregated table with points
write_tsv(expscorestab,"workfiles/CoMM-SU2Cnew_matches_results_ClonDiffPoints.txt")