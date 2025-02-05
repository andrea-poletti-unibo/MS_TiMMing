library(BradleyTerry2)
library(qvcalc)
library(data.table)
library(tidyverse)


import <- fread("workfiles/Phylogic_GoodQualSamp_matches_results_ClonDiffPoints.txt") 


import2 <- import %>% select(CNA1 = player1, CNA2= player2, win1, win2)


# round scores
import2$win1 <- import2$win1 %>% round()
import2$win2 <- import2$win2 %>% round()


matchesB <- import2 %>% as.data.frame()

# exclude <- c("amp_chr_5q", "amp_chr_7p")
exclude <- c()

matchesB <- matchesB %>% filter((!CNA1 %in% exclude) & (!CNA2 %in% exclude))


all_levels <- c(matchesB$CNA1,matchesB$CNA2) %>% unique %>% sort


matchesB$CNA1 <- factor(matchesB$CNA1, levels = all_levels)
matchesB$CNA2 <- factor(matchesB$CNA2, levels = all_levels)

B.CNA.Model <- BTm(cbind(win1, win2), 
                   CNA1, CNA2, 
                   ~ CNA, 
                   id = "CNA", 
                   refcat = "HyperDiploidy", # reference category
                   br = T, # bias reduced
                   data = matchesB)

B.CNA.Model
BTabilities(B.CNA.Model) %>% as.data.frame() 

DF <- BTabilities(B.CNA.Model) %>% as.data.frame()
nas <- DF[is.na(DF$ability), ] %>% rownames()

exclude2 <- c(exclude, nas) %>% unique()


matchesB <- matchesB %>% filter((!CNA1 %in% exclude2) & (!CNA2 %in% exclude2))


all_levels <- c(matchesB$CNA1,matchesB$CNA2) %>% unique %>% sort


matchesB$CNA1 <- factor(matchesB$CNA1, levels = all_levels)
matchesB$CNA2 <- factor(matchesB$CNA2, levels = all_levels)

B.CNA.Model <- BTm(cbind(win1, win2), 
                   CNA1, CNA2, 
                   ~ CNA, 
                   id = "CNA", 
                   refcat = "HyperDiploidy", # reference category
                   br = T, # bias reduced
                   data = matchesB)



ab.qv <- BradleyTerry2::qvcalc(BTabilities(B.CNA.Model))


plot(ab.qv)

ab_df <- as.data.frame(ab.qv$qvframe) %>% rownames_to_column("cna") 
# ab_df <- ab_df %>% filter(cna!= "chr_13q_amp" )

ORD <- ab_df$estimate %>% order


ab_df <- as.data.frame(ab.qv$qvframe) %>% rownames_to_column("cna") 
ab_df$CNA_type <- ifelse(grepl("amp|Hyp",ab_df$cna),"Amp","Del")

ORD <- ab_df$estimate %>% order

ab_df$cna <- factor(ab_df$cna, levels = ab_df$cna[ORD] )

ab_df$timing <- 0- ab_df$estimate



#================ plot timing MAP ===================

gg <- ab_df %>% 
  ggplot(aes(x = cna, y = timing )) + 
  geom_errorbar(aes(ymin=timing-quasiSE, ymax=timing+quasiSE)) +
  geom_point(aes(color=CNA_type)) + 
  xlab("CNA") + ylab("Timing Estimate") +
  ggtitle("Bradley-Terry CNAs timing estimates on SU2C good qual samps", subtitle = "HD call and Clonal difference points" )+
  coord_flip() + theme(legend.position = "bottom") + ylim(-2,12)

gg


ggsave("plots/timing_maps/Phylogic_comparison/GoodQualSamp_timing_map.png", gg, width = 7, height = 10, units = "in", dpi = 300)

write_tsv(ab_df, "workfiles/Phylogic_comparison/BT_model_timing_GoodQualSamp.txt")

Sys.info()
