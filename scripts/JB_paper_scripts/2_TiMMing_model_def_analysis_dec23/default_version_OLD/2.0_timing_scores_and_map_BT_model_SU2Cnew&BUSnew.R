
library(BradleyTerry2)
library(qvcalc)
library(data.table)
library(tidyverse)


import <- fread("workfiles/SU2Cnew&BUSnew_matches_results_ClonDiffPoints.txt") 


import2 <- import %>% select(CNA1 = player1, CNA2= player2, win1, win2)


# round scores
import2$win1 <- import2$win1 %>% round()
import2$win2 <- import2$win2 %>% round()


matchesB <- import2 %>% as.data.frame()

exclude <- c("amp_chr_5q", "amp_chr_7p")

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
  ggtitle("Bradley-Terry CNAs timing estimates on SMM SU2C&Bustoros cohort", subtitle = "HD call and Clonal difference points" )+
  coord_flip() + theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm"))

gg

#_____ load callset _____

cs <- fread("workfiles/callset_SU2Cnew&BUSnew_HDcall.txt")

cs2 <- cs %>% filter(call != -999)
cs_w <- reshape2::dcast(cs2, variable ~ call ) %>% setNames(c("variable","no_call","calls"))
cs_w$variable <- factor(cs_w$variable, levels = ab_df$cna[ORD], ordered = T)
cs_w$type <- ifelse(cs_w$variable %>% str_detect("del"), "del","amp")

ggF <- cs_w %>% filter(!is.na(variable)) %>% 
  ggplot(aes(variable, calls, fill=type, label= paste0(calls, "(", round((calls/(calls+no_call))*100,1) ,"%)"))) + 
  geom_bar(stat="identity") + 
  geom_text(stat="identity", size=2, vjust=0.5, hjust= -0.1, position="stack")+
  coord_flip() + 
  theme(legend.position = "bottom", 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(), 
        plot.margin = unit(c(0,1,0,0), "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(cs_w$calls) + 60)) +
  xlab("") + ggtitle(label = "", subtitle = "Frequency of alterations")

ggF

#_____ composite plot _____

library(cowplot)

grid <- plot_grid(gg, ggF, ncol = 2, nrow=1, align = "h", rel_widths = c(2,1))

grid

# save plot
ggsave(plot = grid, filename =  "plots/timing_maps/SU2Cnew&BUSnew_timing_map_HDcall_clonDiff.png", width = 14, height = 12, units = "in")


#======================== only frequent alterations ==============================

hifreq <- cs_w$variable[cs_w$calls>10] %>% as.character()
hifreq

gg <- ab_df %>% filter(cna %in% hifreq) %>% 
  ggplot(aes(x = cna, y = timing )) + 
  geom_errorbar(aes(ymin=timing-quasiSE, ymax=timing+quasiSE)) +
  geom_point(aes(color=CNA_type)) + 
  xlab("CNA") + ylab("Timing Estimate") +
  ggtitle("Bradley-Terry CNAs timing estimates on SMM SU2C&Bustoros cohort", subtitle = "HD call and Clonal difference points - HighFreq > 10" )+
  coord_flip() + theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm"))

gg

ggF <- cs_w %>%  filter(variable %in% hifreq) %>% 
  ggplot(aes(variable, calls, fill=type, label= paste0(calls, "(", round((calls/(calls+no_call))*100,1) ,"%)"))) + 
  geom_bar(stat="identity") + 
  geom_text(stat="identity", size=2, vjust=0.5, hjust= -0.1, position="stack")+
  coord_flip() + 
  theme(legend.position = "bottom", 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(), 
        plot.margin = unit(c(0,1,0,0), "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(cs_w$calls) + 60)) +
  xlab("") + ggtitle(label = "", subtitle = "Frequency of alterations")


ggF


library(cowplot)

grid <- plot_grid(gg, ggF, ncol = 2, nrow=1, align = "h", rel_widths = c(2,1))
grid

# save plot
ggsave(plot = grid, filename = "plots/timing_maps/SU2Cnew&BUSnew_timing_map_HiFreq_HDcall_clonDiff.png", width = 14, height = 12, units = "in")


# save dataset of timing and freq
ab_df$cna <- ab_df$cna %>% as.character()
timing_freq <- left_join(ab_df, cs_w, by=c("cna"="variable")) %>% select(-type) %>% mutate(observ=calls + no_call)

write_tsv(timing_freq,"workfiles/BTmodel_SU2Cnew&BUSnew_HDcall_clonDiff.txt") 





##################### additional analysis #########################


# correlation between frequency and timing
merge <- left_join(ab_df %>% mutate(cna=as.character(cna)), 
                   cs_w %>% mutate(variable=as.character(variable)) , 
                   by=c("cna"="variable"))

merge$Frequency <- merge$calls / (merge$calls + merge$no_call)

merge %>% ggplot(aes(timing, Frequency , colour=type)) + geom_point() + 
  geom_smooth(method = lm, alpha=0.1, size =0.3) +
  ggpubr::stat_cor() + xlab("Timing Estimate (TE)") + ggtitle("Correlation between Timing Estimates and Frequencies of alterations")

cor.test(merge$timing, merge$calls)

ggsave("plots/TE_Freq_correlation_SU2Cnew.png", width = 7, height = 7)


# clustering of timing windows

library(NbClust)

ab_df$cna <- factor(ab_df$cna, levels = ab_df$cna[ORD] )

ab_df %>% filter(cna %in% hifreq) %>% .$timing

NB <- ab_df %>% filter(cna %in% hifreq) %>% .$timing %>% NbClust(method="kmeans" ,index = "silhouette")

time_window <- NB$Best.partition %>% as.factor()

gg <- ab_df %>% filter(cna %in% hifreq) %>%
  ggplot(aes(x = cna, y = timing )) + 
  geom_errorbar(aes(ymin=timing-quasiSE, ymax=timing+quasiSE)) +
  geom_point(aes(color=time_window)) + 
  xlab("CNA") + ylab("Timing Estimate") +
  ggtitle("Bradley-Terry CNAs timing estimates", subtitle = "Time-windows intervals" )+
  coord_flip() + theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm"))

gg


