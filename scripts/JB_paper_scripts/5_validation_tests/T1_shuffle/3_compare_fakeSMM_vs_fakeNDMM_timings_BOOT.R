library(tidyverse)
library(data.table)
library(ggrepel)
library(gridExtra)


com <- fread("results/JB_paper/def_aug23/validation_tests/T1_shuffle/Fake_NDMM_boots/fakeNDMM_timings_20boot.txt")

com_boot <- com %>% group_by(cna) %>% summarise(med_timing=median(timing), mean_timing=mean(timing), sd = sd(timing), min= min(timing), max= max(timing), first_q= quantile(timing)[2], third_q= quantile(timing)[4],  n=n()) %>% arrange(med_timing)


smm <- fread("results/JB_paper/def_aug23/validation_tests/T1_shuffle/Fake_SMM_boots/fakeSMM_timings_20boot.txt")

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


hifreq <- cs_w$variable[cs_w$calls>15] %>% as.character()
hifreq

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
ggsave("results/JB_paper/def_aug23/validation_tests/T1_shuffle/compare_Viz1_timings_BOOTS.png", width = 7, height = 6, dpi = 300)



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
ggsave("results/JB_paper/def_aug23/validation_tests/T1_shuffle/compare_Viz2_timings_BOOTS.png", width = 7, height = 6, dpi = 300)




# #____________________ VIZ 3 _________________________________

# m2$cna_order <- m2$cna %>% factor(levels = m2$cna[order(m2$timing_diff)])

# m2 <- left_join(m2, pvals, by = "cna")
# m2$qval_signif <- ifelse(m2$qval < 0.05, "yes", "no")

# # BAR PLOT
# p1 <- m2 %>% ggplot(aes( cna_order ,abs(timing_diff))) + 
# geom_bar(stat = "identity", aes(fill=qval_signif), colour="black") +
# theme_bw()+
# theme(axis.text.x = element_blank(), 
# axis.ticks.x = element_blank(),
# axis.title.x = element_blank(),
# axis.title.y = element_text(size = 10),
# axis.text.y = element_text(size = 9),
# legend.position = "none") + ylab("abs Clonality difference") +
# # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
# # set plot limit in y to 5 without changing the data
# coord_cartesian(ylim = c(0, 2.7), expand = FALSE) 

# p1


# # Frequency plot
# allfreq2$cna_order <- allfreq2$variable %>% factor(levels = m2$cna[order(m2$timing_diff)])
# allfreq2$cohort <- factor(allfreq2$cohort, levels = c("SMM", "NDMM"))


# p2 <- allfreq2 %>% ggplot(aes(cna_order, frequency, fill=cohort)) + 
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     theme_bw() + 
#     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
#     # remove x lab
#     axis.title.x = element_blank(),
#     ) + 
#     scale_fill_manual(values = c( "grey40", "grey80")) + 
#     theme(legend.position = "top") +
#     # y tick every 0.1
#     scale_y_continuous(breaks = seq(0, 1, 0.1)) +
#     # remove x labels
#     theme() +
#     # start plot area at exactly 0
#     coord_cartesian(expand = FALSE, ylim = c(0, 0.63))
    
# p2


# # BAR PLOT
# ml2 <- ml
# ml2$cna_order <- factor(ml2$cna, levels = m2$cna[order(m2$timing_diff)])

# ml2 <- ml2 %>% left_join(pvals, by = "cna") 


# ml2 <- ml2 %>% group_by(cna) %>% mutate(timing_diff=diff(timing))

# ml2$qval_signif <- ifelse(ml2$qval < 0.05, "yes", "no")
# ml2$qval_signif[is.na(ml2$qval_signif)] <- "no"

# p3 <- ml2 %>% ggplot(aes(cna_order, timing)) +
#     geom_errorbar(aes(ymin = timing - 2 * sd, ymax = timing + 2 * sd, fill = model, colour = qval_signif),
#         colour = "grey50", size = 0.3, position = position_dodge(width = 0.8)
#     ) +
#     geom_point(position = position_dodge(width = 0.8), aes(shape=model, fill = model, colour = qval_signif), size=1.7) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
#         legend.position = "bottom"
#     ) +
#     # scale_colour_gradient2(low = "red", mid = "gray50", high = "blue", 
#     #     midpoint = 0, limits=c(-2,2), oob= scales::squish ) +
#     scale_shape_manual(values = c(1, 19), labels = c("MM", "SMM")) +
#     scale_fill_discrete(labels = c("MM", "SMM")) +
#     xlab("Alteration") + ylab("BT-score") +
#     # remove "chr from x labels"
#     scale_x_discrete(labels = function(x) str_replace_all(x, "_chr|_amp|_del", "") %>% str_replace("del_","-") %>% str_replace("amp_","+")) 
    
# p3

# # create a grid plot

# # grid.arrange( p2, p1, p3, ncol = 1, heights = c(1, 1, 2))

# GR <- grid.arrange(p1, p3, ncol = 1, heights = c(1, 1, 2))
# GR

# # ggsave("results/JB_paper/def_aug23/boot_approach/compare_Viz3_SMM_vs_CoMMpass_timings_BOOT.png", plot = GR, width = 10, height = 7, dpi = 300)


# library(cowplot)
# GR3 <- plot_grid(p2, p1, p3, ncol = 1, align = "v", rel_heights = c(1,0.8,1.5))
# GR3

# # save GR3
# # ggsave("results/JB_paper/def_aug23/boot_approach/combined3plots_SMM_vs_CoMMpass_timings_BOOT.png", plot = GR3, width = 8, height = 7, dpi = 300)


# #____________________ VIZ 4 _________________________________

# m2 %>% ggplot(aes(timing_diff, -log(qval.x))) + geom_point()

# m2 %>% ggplot(aes(timing_diff, -log(qval.x))) + geom_point(aes(colour=timing_diff)) + geom_hline(yintercept = -log(0.05), linetype=2) + geom_vline(xintercept = c(-1,1), linetype=2) +
# scale_colour_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0, limits=c(-2,2), oob= scales::squish ) + geom_label_repel(data= m2 %>% filter(qval<0.05 & abs(timing_diff)>1), aes(label = lab), size = 3, min.segment.length = 3) +
# theme_bw() + xlab("Timing difference") + ylab("-log(q-value)") 

# # ggsave("results/JB_paper/def_aug23/boot_approach/compare_Viz4_SMM_vs_CoMMpass_timings_BOOTS.png", width = 6, height = 5, dpi = 300)



# #____________________ VIZ 5: diffs ______________________________

# v5f <- allfreq2 %>% dcast(variable ~ cohort, value.var = "frequency") %>% mutate(diff = NDMM - SMM) %>% arrange(desc(diff))

# m2

# v5 <- left_join(v5f, m2, by = c("variable" = "cna"))

# v5 %>% ggplot(aes(diff, timing_diff)) + 
#     geom_smooth(method = "lm", alpha=0.3, colour="grey50") + 
#     ggpubr::stat_cor( method = "spearman", label.x = 0.08, label.y = 2) +
#     geom_point(aes(colour=qval_signif)) + 
#     geom_hline(yintercept = 0, linetype=2) + 
#     geom_vline(xintercept = 0, linetype=2) + 
#     # geom_abline(intercept = 0, slope = 10, color = "red") +
#     geom_text_repel(data= v5 %>% filter(qval_signif=="yes" | abs(diff)>0.05 | abs(timing_diff)>1), aes(label = lab), size = 3, max.overlaps=100)  + 
#     theme_bw() + 
#     xlab("Frequency difference") + ylab("Timing difference") +
#     # legend bottom
#     theme(legend.position = "bottom") + 
#     ggtitle("Frequency difference vs Timing difference (diff = NDMM - SMM)", subtitle = "events with a frequency > 5% in at least one cohort are shown (n=44)")

# # ggsave("results/JB_paper/def_aug23/boot_approach/compare_Viz5_SMM_vs_CoMMpass_timings_BOOTS.png" , width = 6, height = 6, dpi = 300)

# # save diff data
# write_tsv(v5, "workfiles/data_diff-timing_vs_diff-frequency_SMM_CoMMpass.txt")


# # strict frequency version

# v5_strict_hifreq <- v5 %>% filter( NDMM > 0.05 & SMM > 0.05)

# v5_strict_hifreq %>% ggplot(aes(diff, timing_diff)) + 
#     geom_smooth(method = "lm", alpha=0.3, colour="grey50") + 
#     ggpubr::stat_cor( method = "spearman", label.x = 0.08, label.y = 2) +
#     geom_point(aes(colour=qval_signif)) + 
#     geom_hline(yintercept = 0, linetype=2) + 
#     geom_vline(xintercept = 0, linetype=2) + 
#     # geom_abline(intercept = 0, slope = 10, color = "red") +
#     geom_text_repel(data= v5_strict_hifreq %>% filter(qval_signif=="yes" | abs(diff)>0.05 | abs(timing_diff)>1), aes(label = lab), size = 3, max.overlaps=100)  + 
#     theme_bw() + 
#     xlab("Frequency difference") + ylab("Timing difference") +
#     # legend bottom
#     theme(legend.position = "bottom") +
#     ggtitle("Frequency difference vs Timing difference (diff = NDMM - SMM)", subtitle = "STRINGENT frequency: events > 5% in BOTH cohorts are shown (n=32)")
    
# # ggsave("results/JB_paper/def_aug23/boot_approach/compare_Viz5.2_SMM_vs_CoMMpass_timings_BOOTS.png" , width = 6, height = 6, dpi = 300)




# v5$timing_diff %>% density() %>% plot
# v5$diff %>% density() %>% plot

# # test normality
# shapiro.test(v5$diff)
# shapiro.test(v5$timing_diff)


# #___________ export definitive table with freq, timing and p/q vals for paper ________

# # Genomic alteration	NDMM freq	SMM freq	freq diff	BT-score NDMM	BT-score SMM	Clonality difference	p-val	q-val (BH)	q-val signif

# v5 %>% names

# exptab <- v5 %>% select(`Genomic alteration` = variable, 
#     `NDMM freq` = NDMM, 
#     `SMM freq` = SMM, 
#     `freq diff` = diff, 
#     `BT-score NDMM` = med_timing_com, 
#     `BT-score SMM` = med_timing_smm, 
#     `Clonality difference` = timing_diff, 
#     `p-val`= pval.x, 
#     `q-val (BH)`=qval.x, 
#     `q-val signif`= qval_signif)

# write_tsv(exptab, "results/JB_paper/all_alterations_freq_timing_diff_table.txt")
