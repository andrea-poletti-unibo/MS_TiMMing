library(data.table)
library(tidyverse)
library(BradleyTerry2)
library(qvcalc)

files <- list.files("workfiles/saturation_analysis_10boot/matches", full.names = T)

f=files[1]

res_all <- data.frame()

for (f in files) {
  sub <- f %>%
    str_extract("n_[0-9]+") %>%
    str_remove("n_") %>%
    as.numeric()

  boot <- f %>%
    str_extract("boot[0-9]+") %>%
    str_remove("boot") %>%
    as.numeric()
  message("sub: ", sub," / boot: ", boot)


  import <- fread(f)

  import2 <- import %>% select(CNA1 = player1, CNA2 = player2, win1, win2)

  # round scores
  import2$win1 <- import2$win1 %>% round()
  import2$win2 <- import2$win2 %>% round()

  matchesB <- import2 %>% as.data.frame()

  exclude <- c("amp_chr_5q", "amp_chr_7p")

  matchesB <- matchesB %>% filter((!CNA1 %in% exclude) & (!CNA2 %in% exclude))

  all_levels <- c(matchesB$CNA1, matchesB$CNA2) %>%
    unique() %>%
    sort()

  matchesB$CNA1 <- factor(matchesB$CNA1, levels = all_levels)
  matchesB$CNA2 <- factor(matchesB$CNA2, levels = all_levels)

  B.CNA.Model <- BTm(cbind(win1, win2),
    CNA1, CNA2,
    ~CNA,
    id = "CNA",
    refcat = "HyperDiploidy", # reference category
    br = T, # bias reduced
    data = matchesB
  )

  B.CNA.Model
  DF <- BTabilities(B.CNA.Model) %>% as.data.frame()
  nas <- DF[is.na(DF$ability), ] %>% rownames()

  exclude2 <- c(exclude, nas) %>% unique()

  matchesB <- matchesB %>% filter((!CNA1 %in% exclude2) & (!CNA2 %in% exclude2))

  all_levels <- c(matchesB$CNA1, matchesB$CNA2) %>%
    unique() %>%
    sort()

  matchesB$CNA1 <- factor(matchesB$CNA1, levels = all_levels)
  matchesB$CNA2 <- factor(matchesB$CNA2, levels = all_levels)

  B.CNA.Model <- BTm(cbind(win1, win2),
    CNA1, CNA2,
    ~CNA,
    id = "CNA",
    refcat = "HyperDiploidy", # reference category
    br = T, # bias reduced
    data = matchesB
  )

  ab.qv <- BradleyTerry2::qvcalc(BTabilities(B.CNA.Model))

  ab_df <- as.data.frame(ab.qv$qvframe) %>% rownames_to_column("cna")

  ORD <- ab_df$estimate %>% order()

  ab_df <- as.data.frame(ab.qv$qvframe) %>% rownames_to_column("cna")
  ab_df$CNA_type <- ifelse(grepl("amp|Hyp", ab_df$cna), "Amp", "Del")

  ORD <- ab_df$estimate %>% order()

  ab_df$cna <- factor(ab_df$cna, levels = ab_df$cna[ORD])

  ab_df$timing <- 0 - ab_df$estimate

  ab_df$sub <- sub

  ab_df$boot <- boot

  res_all <- rbind(res_all, ab_df)

#   gg <- ab_df %>%
#     ggplot(aes(x = cna, y = timing)) +
#     geom_errorbar(aes(ymin = timing - quasiSE, ymax = timing + quasiSE)) +
#     geom_point(aes(color = CNA_type)) +
#     xlab("CNA") +
#     ylab("Timing Estimate") +
#     ggtitle("Bradley-Terry CNAs timing estimates on CoMMpass new cohort", subtitle = paste("subsample n = ", sub)) +
#     coord_flip() +
#     theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))

#   gg

#   dir.create("plots/timing_maps/saturations/", showWarnings = F)
#   ggsave(plot = gg, filename = paste0("plots/timing_maps/saturations/sub_", sub, "CoMM-SU2Cnew_timing_map_HDcall_clonDiff.png"), width = 14, height = 12, units = "in")
}


dir.create("results/saturation_analysis_10boot/", showWarnings = F)
write_tsv(res_all, "results/saturation_analysis_10boot/CoMM-SU2Cnew_timings_saturation10boot.txt")



res_all <- fread("results/saturation_analysis_10boot/CoMM-SU2Cnew_timings_saturation10boot.txt")

res_all <- res_all %>% filter(timing > -10) # remove outliers (model errors in amp_chr_3p at very low sub size)

res_all$cna_f <- factor(res_all$cna, levels = str_sort(res_all$cna %>% unique, numeric = T ))

res_all <- res_all %>% group_by(cna, sub ) %>% mutate(sd_group = sd(timing), median_group = median(timing))

res_all$sd_group[res_all$sd_group>3] <- 3

res_all %>% filter(cna_f == "del_chr_13q") %>% 
    ggplot(aes(sub, timing, group=sub)) + 
    geom_boxplot(aes(fill=sd_group,colour=sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.8, width=5) +
    facet_wrap(~cna_f, scales = "fixed") +
    theme(legend.position = "bottom") + ylim(-2,12) +
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1)



map <- res_all %>% 
    ggplot(aes(sub, timing, group=sub)) + 
    geom_boxplot(aes(fill=sd_group, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    facet_wrap(~cna_f, scales = "fixed") +
    theme(legend.position = "bottom") + ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1)

map

dir.create("plots/saturations_10boot/", showWarnings = F)
# ggsave(plot = map, filename = "plots/saturations_10boot/CoMM_saturation_10boot.png" , width = 36, height = 36, units = "in")




#_____ load callset _____

cs <- fread("workfiles/callset_CoMM-SU2Cnew_HDcall.txt")

cs2 <- cs %>% filter(call != -999)
cs_w <- reshape2::dcast(cs2, variable ~ call ) %>% setNames(c("variable","no_call","calls"))
cs_w$type <- ifelse(cs_w$variable %>% str_detect("del"), "del","amp")

cs_w$tot <- cs_w$no_call + cs_w$calls
cs_w$perc <- cs_w$calls / cs_w$tot




HI_freq <- cs_w %>% filter(perc > 0.23) %>% pull(variable)

satgg2 <- res_all %>% filter(cna %in% HI_freq) %>% 
ggplot(aes(sub, timing, color = cna)) + 
geom_errorbar(aes(ymin=timing-quasiSE, ymax=timing+quasiSE), width=3, size = 0.2) +
geom_point() + 
geom_line(size = 0.3) + theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm"))

satgg2

# ggsave(plot = satgg2, filename = "plots/saturations/CoMM-SU2Cnew_saturation_HIfreq.png", width = 16, height = 10, units = "in")



LOW_freq <- cs_w %>% filter(perc < 0.10 & perc> 0.08) %>% pull(variable)

satgg3 <- res_all %>% filter(cna %in% LOW_freq) %>% 
ggplot(aes(sub, timing, color = cna)) + 
geom_errorbar(aes(ymin=timing-quasiSE, ymax=timing+quasiSE), width=3, size = 0.2) +
geom_point() + 
geom_line(size = 0.3) + theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm"))

satgg3

# ggsave(plot = satgg3, filename = "plots/saturations/CoMM-SU2Cnew_saturation_LOfreq.png", width = 16, height = 10, units = "in")



#_____ merge with callset _____

merge <- res_all %>% left_join(cs_w, by = c("cna" = "variable"))

dat_anno <- merge %>% group_by(cna_f) %>%mutate(perc = mean(perc)) %>% ungroup() %>% select(cna_f, perc) %>% distinct()

map <- merge %>% 
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    facet_wrap(~cna_f, scales = "fixed") +
    geom_text(data=dat_anno, mapping = aes(x = 500, y = -1, label = paste0("freq = ", perc %>% round(3))) , size=3) +
    theme(legend.position = "bottom") + ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1)

map 

# ggsave(plot = map, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq.png", width = 36, height = 36, units = "in")

# ggsave(plot = map, filename = "plots/saturations/CoMM-SU2Cnew_saturation_freq.png", width = 18, height = 18, units = "in")




#_____ find saturation cutoffs _____

cuts <- data.frame()
a=dat_anno$cna_f[1]

for (a in dat_anno$cna_f){

  message(a)
  df <- merge %>% filter(cna_f == a) %>% arrange(sub) %>% group_by(sub) %>% summarise(sd=mean(sd_group))
  
  i=1
  for (i in 1:nrow(df)){
    s <- df$sub[i]
    message(s)
    df2 <- df %>% filter(sub >= s)

    if (all(df2$sd < 0.5)){
      res <- data.frame(cna_f = a, sub = s)
      cuts <- rbind(cuts, res)
      break
    }
  }

}

export <- left_join(cuts, dat_anno, by = "cna_f") %>% select(alteration= cna_f, stability_cutoff=sub, freq= perc)

write_tsv(export, "results/saturation_analysis_10boot/Freq_vs_saturation_cutoffs.txt")




#_____ plot saturation cutoffs - DEF FIGURE _____



map <- merge %>% 
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    geom_vline(data= cuts, mapping = aes(xintercept = sub), color = "purple") +
    facet_wrap(~cna_f, scales = "fixed") +
    geom_text(data=dat_anno, mapping = aes(x = 400, y = -1, label = paste0("freq = ", perc %>% round(3))) , size=5) +
    geom_text(data=cuts, mapping = aes(x = sub+125, y = 9, label = paste0("sat = ", sub)) , size=5) +
    theme(legend.position = "bottom") + ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    labs(x = "Subsample size", y = "BT-score") 




ggsave(plot = map, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs.png", width = 36, height = 36, units = "in")
ggsave(plot = map, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs.pdf", width = 36, height = 36, units = "in")

ggsave(plot = map, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs.svg", width = 36, height = 36, units = "in")



write_tsv(merge, "workfiles/saturation_cutoffs_plot_boot.txt")




################### after revision version: better readability ####################

map2 <- merge %>% 
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    geom_vline(data= cuts, mapping = aes(xintercept = sub), color = "purple") +
    facet_wrap(~cna_f, scales = "fixed", ncol = 7) +
    geom_text(data=dat_anno, mapping = aes(x = 400, y = -1, label = paste0("freq = ", perc %>% round(3))) , size=5) +
    geom_text(data=cuts, mapping = aes(x = 400, y = 10, label = paste0("sat = ", sub)) , size=5) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 10), 
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 10, face = "bold")) + 
    ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    labs(x = "Subsample size", y = "BT-score")


ggsave(plot = map2, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v2.pdf", width = 12, height = 20, units = "in")
ggsave(plot = map2, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v2.svg", width = 12, height = 20, units = "in")



################ JB suggestions ####################

mergeF <- merge %>% filter(cna_f != "HyperDiploidy")

#____ rename labels of CNAs _____
mergeF$cna_f  %>% table
mergeF$cna_f2 <- mergeF$cna_f %>% str_replace("del_chr_","-") %>% str_replace("amp_chr_","+") %>% str_replace("_del","") %>% str_replace("_amp","")
mergeF$cna_f2  %>% table


dat_anno$cna_f2 <- dat_anno$cna_f %>% str_replace("del_chr_","-") %>% str_replace("amp_chr_","+") %>% str_replace("_del","") %>% str_replace("_amp","")
cuts$cna_f2 <- cuts$cna_f %>% str_replace("del_chr_","-") %>% str_replace("amp_chr_","+") %>% str_replace("_del","") %>% str_replace("_amp","")

#___ split amps and dels _____
m_amp <- mergeF %>% filter(CNA_type %>% str_detect("Amp")) 
m_del <- mergeF %>% filter(CNA_type %>% str_detect("Del"))

m_amp$cna_f2 <- factor(m_amp$cna_f2 , levels = str_sort(m_amp$cna_f2  %>% unique, numeric = T ))
m_del$cna_f2 <- factor(m_del$cna_f2 , levels = str_sort(m_del$cna_f2  %>% unique, numeric = T ))

dat_anno_amp <- dat_anno %>% filter(cna_f %>% str_detect("amp"))
dat_anno_del <- dat_anno %>% filter(cna_f %>% str_detect("del"))

dat_anno_amp$cna_f2 <- factor(dat_anno_amp$cna_f2 , levels = str_sort(dat_anno_amp$cna_f2  %>% unique, numeric = T ))
dat_anno_del$cna_f2 <- factor(dat_anno_del$cna_f2 , levels = str_sort(dat_anno_del$cna_f2  %>% unique, numeric = T ))

cuts_amp <- cuts %>% filter(cna_f %>% str_detect("amp"))
cuts_del <- cuts %>% filter(cna_f %>% str_detect("del"))

cuts_amp$cna_f2 <- factor(cuts_amp$cna_f2 , levels = str_sort(cuts_amp$cna_f2  %>% unique, numeric = T ))
cuts_del$cna_f2 <- factor(cuts_del$cna_f2 , levels = str_sort(cuts_del$cna_f2  %>% unique, numeric = T ))


#_______ single alteration plot ________
m_amp$cna_f2 %>% unique
ALT <- "+3q_TERC"

m_amp %>% filter(cna_f2 == ALT) %>%
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    geom_jitter(size=0.8, width=5, colour="grey30", alpha=0.7) +
    geom_vline(data= cuts_amp %>% filter(cna_f2==ALT), mapping = aes(xintercept = sub), color = "purple") +
    geom_text(data=cuts_amp %>% filter(cna_f2==ALT), mapping = aes(x = sub+100, y = 6, label = paste0("cutoff = ", sub)) , size=5) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(face = "bold", size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")) + 
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    labs(x = "Sample size", y = "BT-score", colour = "Standard Deviation", fill = "Standard Deviation") +
    ggtitle(paste(ALT), 
      subtitle = paste0("Alteration frequency = ", dat_anno_amp %>% filter(cna_f2==ALT) %>% pull(perc) %>% round(3) %>% `*`(100),"%"))


#_______________________________________

map3a <- m_amp %>% 
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    # geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    geom_vline(data= cuts_amp, mapping = aes(xintercept = sub), color = "purple") +
    facet_wrap(~cna_f2, scales = "fixed", ncol = 7) +
    geom_text(data=dat_anno_amp, mapping = aes(x = 400, y = -1, label = paste0("frequency = ", perc %>% round(3))) , size=3) +
    geom_text(data=cuts_amp, mapping = aes(x = 400, y = 10, label = paste0("cutoff = ", sub)) , size=3) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10), 
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 10, face = "bold")) + 
    ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    labs(x = "Sample size", y = "BT-score", colour = "Standard Deviation", fill = "Standard Deviation")

ggsave(plot = map3a, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v3_amp.pdf", width = 10, height = 10, units = "in")
ggsave(plot = map3a, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v3_amp.svg", width = 10, height = 10, units = "in")


map3d <- m_del %>% 
    ggplot(aes(sub, timing)) + 
    geom_boxplot(aes(fill=sd_group, group=sub, colour= sd_group), alpha=0.2, outlier.shape = NA) +
    # geom_jitter(size=0.5, width=5, colour="grey30", alpha=0.5) +
    geom_vline(data= cuts_del, mapping = aes(xintercept = sub), color = "purple") +
    facet_wrap(~cna_f2, scales = "fixed", ncol = 7) +
    geom_text(data=dat_anno_del, mapping = aes(x = 400, y = -1, label = paste0("frequency = ", perc %>% round(3))) , size=3) +
    geom_text(data=cuts_del, mapping = aes(x = 400, y = 10, label = paste0("cutoff = ", sub)) , size=3) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10), 
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 10, face = "bold")) + 
    ylim(-2,12) +  
    scale_color_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    scale_fill_gradient2(low="#079e07", mid="yellow", high="red", midpoint = 1) +
    labs(x = "Sample size", y = "BT-score", colour = "Standard Deviation", fill = "Standard Deviation")

ggsave(plot = map3d, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v3_del.pdf", width = 10, height = 10, units = "in")
ggsave(plot = map3d, filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_v3_del.svg", width = 10, height = 10, units = "in")






#=========================== further analysis ===========================

m2 <- merge %>% left_join(cuts, by = c("cna_f" = "cna_f"))

s2 <- m2 %>% group_by(cna_f) %>% mutate(cut=unique(sub.y), freq=unique(perc)) %>% ungroup() %>% select(cna_f, cut, freq) %>% distinct()

s2 %>% filter(cna_f != "HyperDiploidy")  %>% ggplot(aes(freq, cut)) + geom_point() + geom_smooth(method = "lm") + ggpubr::stat_cor(method = "pearson", label.x = 0.25, label.y = 50, size=5) + 
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = cna_f), size = 3) +
  xlab("CNA frequency") + ylab("Saturation cutoff") + 
  ggtitle("CoMM-SU2Cnew saturation cutoffs vs CNA frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm")) + ylim(-50, NA) + xlim(-0.03,NA)

# ggsave(plot = last_plot(), filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_correlation.png", width = 14, height = 12, units = "in")






s2 %>% filter(cna_f != "HyperDiploidy")  %>% ggplot(aes(freq, cut)) + geom_point() + 
geom_smooth(method="lm", formula= (y ~ log(x))) +
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = cna_f), size = 3) +
  xlab("CNA frequency") + ylab("Saturation cutoff") + 
  ggtitle("CoMM-SU2Cnew saturation cutoffs vs CNA frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm")) + ylim(-50, NA) + xlim(-0.03,NA)



logM <- lm(log(cut) ~ freq, data = s2) %>% summary %>% .$r.squared
linM <- lm(cut ~ freq, data = s2) %>% summary %>% .$r.squared

s2 %>% filter(cna_f != "HyperDiploidy")  %>% ggplot(aes(freq, cut)) + geom_point() + 
geom_smooth(method="lm", formula= (y ~ log(x))) +
geom_smooth(method="lm", formula= (y ~ x), color = "orange") +
annotate("text", x = 0.3, y = 600, label = paste(" R^2 linear =", round(linM, 3)), colour="orange", hjust=0, size=7) +
annotate("text", x = 0.3, y = 550, label = paste(" R^2 log =", round(logM, 3)), colour="blue", hjust=0, size=7) +
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = cna_f), size = 3) +
  xlab("CNA frequency") + ylab("Saturation cutoff") + 
  ggtitle("CoMM-SU2Cnew saturation cutoffs vs CNA frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm")) + ylim(0, NA) + xlim(0,NA)

# ggsave(plot = last_plot(), filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_log_vs_lin_Model.png", width = 10, height = 10, units = "in")





s2$label <- s2$cna_f %>% str_replace("del_chr_", "-") %>% str_replace("amp_chr_", "+") %>% str_replace("_del", "") %>% str_replace("_amp", "")

s2 %>% filter(cna_f != "HyperDiploidy")  %>% ggplot(aes(freq, cut)) + geom_point() + 
geom_smooth(method="lm", formula= (y ~ log(x))) +
annotate("text", x = 0.3, y = 550, label = paste("R^2 ==", round(logM, 3)), colour="blue", hjust=0, size=7, parse=T) +
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = label), size = 2.5) +
  xlab("CNA frequency") + ylab("Saturation cutoff") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  # theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm")) + 
  ylim(0, NA) + xlim(0,NA)




ggsave(plot = last_plot(), filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_logModel.pdf", width = 6, height = 6, units = "in")

ggsave(plot = last_plot(), filename = "plots/saturations_10boot/CoMM_saturation_10boot_freq_cutoffs_logModel.svg", width = 6, height = 6, units = "in")






#========================= Ranking analysis =========================

# compute a ranking of CNAs based on their timing in different subsamples groups
gr <- res_all %>% group_by(sub) %>% mutate(Rank = rank(timing)) %>% ungroup()

#___________ animated ranking ___________

library(gganimate)
library(tweenr)
library(magick)
library(gifski)
library(grid)




p <- ggplot(gr) + #, aes(x = Rank, group = sub, country = as.factor(Code)
        geom_col(aes(x = Rank, y = timing, fill = CNA_type), width = 0.5, color = "black") + # Columns
        coord_flip(clip = "off", expand = FALSE) + # Flip
        labs(title='{closest_state}', x = "", y = "Timing Estimate (TE)") + # Labels
        theme_minimal() + # Theme
        geom_text(aes(x = Rank, y = -3 , label = cna, color = CNA_type), hjust = 1) + # Names
        geom_text(aes(x = Rank, y = timing + 0.3, label = as.character(timing %>% round(3))), hjust = 0, color = "black") + # Values  
        # scale_y_continuous(labels = scales::comma) + # Format y-axis values
        scale_x_reverse() + # Highest values on top
        transition_states(sub, transition_length = 4, state_length = 6) + # Animate
        theme(
            plot.title = element_text(hjust = 0, size = 20),
            plot.margin = margin(0,2,0,3,"cm"),
            axis.text.y  = element_blank()
        )
p

animate(p, fps = 25, duration = 50, width = 800, height = 1000)

anim_save("anim.gif", animation = last_animation(), path = "C:/Users/andre/OneDrive/Desktop/")


#___________ animated ranking of alterations > 5% ___________
gr2 <- gr %>% left_join(cs_w, by = c("cna" = "variable"))
gr2f <- gr2 %>% filter(perc > 0.05)

p <- ggplot(gr2f) + #, aes(x = Rank, group = sub, country = as.factor(Code)
        geom_col(aes(x = Rank, y = timing, fill = CNA_type), width = 0.5, color = "black") + # Columns
        coord_flip(clip = "off", expand = FALSE) + # Flip
        labs(title='{closest_state}', x = "", y = "Timing Estimate (TE)") + # Labels
        theme_minimal() + # Theme
        geom_text(aes(x = Rank, y = -3 , label = cna, color = CNA_type), hjust = 1) + # Names
        geom_text(aes(x = Rank, y = timing + 0.3, label = as.character(timing %>% round(3))), hjust = 0, color = "black") + # Values  
        # scale_y_continuous(labels = scales::comma) + # Format y-axis values
        scale_x_reverse() + # Highest values on top
        transition_states(sub, transition_length = 4, state_length = 6) + # Animate
        theme(
            plot.title = element_text(hjust = 0, size = 20),
            plot.margin = margin(0,2,0,3,"cm"),
            axis.text.y  = element_blank()
        )
p

animate(p, fps = 25, duration = 50, width = 800, height = 1000)

anim_save("anim_m05.gif", animation = last_animation(), path = "C:/Users/andre/OneDrive/Desktop/")


gr3f <- gr2 %>% filter(perc > 0.20)


gr3f  %>% ggplot(aes(sub, Rank, color = cna)) + geom_point() + geom_line()

gr3f  %>% ggplot(aes(sub, timing, color = cna)) + geom_point() + geom_line() + theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))
