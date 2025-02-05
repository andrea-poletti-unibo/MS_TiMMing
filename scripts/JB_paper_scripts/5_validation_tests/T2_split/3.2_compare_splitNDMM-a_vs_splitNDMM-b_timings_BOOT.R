library(tidyverse)
library(data.table)
library(ggrepel)
library(gridExtra)


com <- fread("results/JB_paper/def_aug23/validation_tests/T2_split_NDMM/splitNDMM-a_timings_20boot.txt")

com_boot <- com %>% group_by(cna) %>% summarise(med_timing=median(timing), mean_timing=mean(timing), sd = sd(timing), min= min(timing), max= max(timing), first_q= quantile(timing)[2], third_q= quantile(timing)[4],  n=n()) %>% arrange(med_timing)


smm <- fread("results/JB_paper/def_aug23/validation_tests/T2_split_NDMM/splitNDMM-b_timings_20boot.txt")

smm_boot <- smm %>% group_by(cna) %>% summarise(med_timing=median(timing), mean_timing=mean(timing), sd = sd(timing), min= min(timing), max= max(timing), first_q= quantile(timing)[2], third_q= quantile(timing)[4],  n=n()) %>% arrange(med_timing)

t1 <- paste0(com$cna,"_",com$boot)
t2 <- paste0(smm$cna,"_",smm$boot)

setdiff(t1, t2)


# compare the smm and comm norm timing estimates

m <- merge(com_boot, smm_boot, by = "cna", suffixes = c("_com", "_smm"))

# m <- m %>% filter(cna != "del_19p")

m %>% ggplot(aes(med_timing_com, med_timing_smm)) + 
geom_point() + 
geom_smooth(method = "lm") + 
geom_abline(intercept = 0, slope = 1, color = "red") + ggpubr::stat_cor()



#________________ HI FREQ alterations only ________________

# SMM

cs <- fread("workfiles/callset_SU2Cnew&BUSnew_HDcall.txt")

cs2 <- cs %>% filter(call != -999)

# modify cutoff at 10% instead of 5%
cs2$thresh <- cs2$thresh * 2
cs2$call <- ifelse(cs2$value > cs2$thresh, 1, 0)
cs2$call_FISH <- ifelse(cs2$value > cs2$thresh * 6, 1, 0) %>% as.factor() 
cs2$cohort <- "SMM"

cs_w <- reshape2::dcast(cs2, variable ~ call ) %>% setNames(c("variable","no_call","calls"))
cs_w$type <- ifelse(cs_w$variable %>% str_detect("del"), "del","amp")
cs_w$frequency <- cs_w$calls / (cs_w$calls + cs_w$no_call)
cs_w$cohort <- "SMM"


hifreq <- cs_w$variable[cs_w$calls>20] %>% as.character()
hifreq


# hifreq <- cs_w$variable[cs_w$frequency>0.10] %>% as.character()
# hifreq


# MM

mmcs <- fread("workfiles/callset_CoMM-SU2Cnew_HDcall.txt")

mmcs2 <- mmcs %>% filter(call != -999)

# modify cutoff at 10% instead of 5%
mmcs2$thresh <- mmcs2$thresh * 2
mmcs2$call <- ifelse(mmcs2$value > mmcs2$thresh, 1, 0)
mmcs2$call_FISH <- ifelse(mmcs2$value > mmcs2$thresh * 6, 1, 0) %>% as.factor() 
mmcs2$cohort <- "NDMM"


mmcs_w <- reshape2::dcast(mmcs2, variable ~ call ) %>% setNames(c("variable","no_call","calls"))
mmcs_w$type <- ifelse(mmcs_w$variable %>% str_detect("del"), "del","amp")
mmcs_w$frequency <- mmcs_w$calls / (mmcs_w$calls + mmcs_w$no_call)
mmcs_w$cohort <- "NDMM"


allfreq <- rbind(cs_w, mmcs_w)

allfreq2 <- allfreq %>% filter(variable %in% hifreq)

# plot a cohort grouped histogram of frequencies per variable

allfreq2 %>% ggplot(aes(variable, frequency, fill=cohort)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("orangered", "turquoise3")) + theme_bw() + theme(legend.position = "none")

allfreq3 <- rbind(cs2, mmcs2)

allfreq3 %>% ggplot(aes(variable, call, fill=call_FISH, group=cohort)) + geom_bar(stat="count", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



#_____________ compute p-values and q-values ________________

i=1

pvals <- data.frame()
for ( i in seq_along(m$cna)){

    message(i, " - ", m$cna[i])

    smm_alt <- smm %>% filter(cna == m$cna[i]) %>% .$timing

    com_alt <- com %>% filter(cna == m$cna[i]) %>% .$timing

    p <- wilcox.test(smm_alt, com_alt) %>% .$p.value

    pvals <- rbind(pvals, data.frame(cna = m$cna[i], pval = p))

}

#apply fdr coorection

pvals$qval <- p.adjust(pvals$pval, method = "bonferroni")

signif_p <- pvals %>% filter(pval < 0.05) %>% pull(cna)
signif_q <- pvals %>% filter(qval < 0.05) %>% pull(cna)



m %>% ggplot(aes(sd_com ,sd_smm)) + geom_point() + geom_abline(intercept = 0, slope = 1, color = "red") + ggpubr::stat_cor()

d <- m$sd_smm - m$sd_com
d %>% summary()

m <- m %>% left_join(pvals, by = "cna")

m$timing_diff <- m$med_timing_com - m$med_timing_smm

m2 <- m %>% filter(cna %in% hifreq)

m2$lab <- m2$cna %>% str_remove("_chr") %>% str_remove("_del|_amp")

ml <- m2 %>% select(cna, med_timing_com, med_timing_smm, sd_com, sd_smm) %>%
pivot_longer(cols = c(med_timing_com, med_timing_smm), names_to = "model", values_to = "timing") %>% 
mutate(model = factor(model, levels = c("med_timing_com", "med_timing_smm")))

ml$sd <- ifelse(ml$model == "med_timing_com", ml$sd_com, ml$sd_smm)

ml$cna_order <- factor(ml$cna, levels = m$cna[order(m$med_timing_com, decreasing = T)])

m2$qval_signif <- ifelse(m2$qval < 0.05, "yes", "no")


#____________________ VIZ 1 _________________________________

m2 %>% ggplot(aes(med_timing_com, med_timing_smm)) + 
geom_abline(intercept = c(1,0,-1), slope = 1, color = "grey50", linetype=2) +
geom_errorbar(aes(ymin= med_timing_smm - sd_com , ymax= med_timing_smm  + sd_com ), colour="grey50", size =0.3) +
geom_errorbar(aes(xmin = med_timing_com - sd_smm, xmax= med_timing_com + sd_smm), colour="grey50", size =0.3) +
geom_point(data= m2 %>% filter(qval<0.05), aes(med_timing_com, med_timing_smm, fill=timing_diff), stroke=1, size =2, shape=21, colour="orange") +
geom_point(data= m2 %>% filter(qval>=0.05), aes(med_timing_com, med_timing_smm, fill=timing_diff), stroke=0, size =2, shape=21, colour="black") +
# label points more than 2 from the diagonal
geom_text_repel(data= m2 %>% filter(qval<0.05 & abs(timing_diff)> 1), aes(label = lab), size = 3, min.segment.length = 3) +
geom_text_repel(data= m2 %>% filter(med_timing_com < 1.5 & med_timing_smm < 1.5), aes(label = lab), size = 3) +
# scale colour gradient mid point at 0, high point at 1 and low point at -1
scale_fill_gradient2(low = "blue", mid = "black", high = "red", 
midpoint = 0, limits=c(-2,2), oob= scales::squish ) +
ggpubr::stat_cor(data=m2 %>% filter(timing_diff<5)) +
xlab("MM timing estimate") + ylab("SMM timing estimate") +
theme_bw()

# save plot
ggsave("results/JB_paper/def_aug23/validation_tests/T2_split_NDMM/compare_Viz1_timings_BOOTS.png", width = 7, height = 6, dpi = 300)



#____________________ VIZ 2 _________________________________

ml %>% ggplot(aes(cna_order, timing, colour = model)) +
    geom_errorbar(aes(ymin = timing - 2 * sd, ymax = timing + 2 * sd, group = model),
        colour = "grey50", size = 0.3, position = position_dodge(width = 0.8)) +
    geom_point(position = position_dodge(width = 0.8), size = 1) +
    geom_smooth(method = "lm") +
    theme_bw() +
    ylab("timing estimate") +
    scale_color_discrete(labels = c("MM", "SMM")) +
    coord_flip()

# save plot
ggsave("results/JB_paper/def_aug23/validation_tests/T2_split_NDMM/compare_Viz2_timings_BOOTS.png", width = 7, height = 6, dpi = 300)
