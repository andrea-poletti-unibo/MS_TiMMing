library(tidyverse)
library(data.table)
library(ggrepel)


com <- fread("workfiles/BTmodel_CoMM-SU2Cnew_HDcall_clonDiff.txt")


com$timing_norm <- com$timing %>% scale %>% as.vector
com$timing_rescaled <- scales::rescale(com$timing)

com %>% ggplot(aes(timing, timing_norm)) + geom_point() + 
geom_smooth(method = "lm") + 
geom_abline(intercept = 0, slope = 1, color = "red") + 
theme_bw() + 
theme(axis.text.x = element_text(angle = 90, hjust = 1))


smm <- fread("workfiles/BTmodel_SU2Cnew&BUSnew_HDcall_clonDiff.txt")

smm$timing_norm <- smm$timing %>% scale %>% as.vector
smm$timing_rescaled <- scales::rescale(smm$timing)

smm %>% ggplot(aes(timing_rescaled , timing)) + geom_point() 



# compare the smm and comm norm timing estimates

m <- merge(com, smm, by = "cna", suffixes = c("_com", "_smm"))

m <- m %>% filter(cna != "del_19p")

m %>% ggplot(aes(timing_com, timing_smm)) + 
geom_point() + 
geom_smooth(method = "lm") + 
geom_abline(intercept = 0, slope = 1, color = "red") + ggpubr::stat_cor()

m$cna <- m$cna %>% str_remove("chr_")

#================== rescaling ==================

m$timing_diff_rescaled <- m$timing_rescaled_com - m$timing_rescaled_smm

# identify outliers in the timing difference based on IQR rule
IQR <- m$timing_diff_rescaled %>% IQR %>% print

# compute the upper and lower limits
upper <- IQR * 1.5
lower <- IQR * -1.5





# sort cna by timing difference, order factor
m$cna_order <- factor(m$cna, levels = m$cna[order(m$timing_diff_rescaled)])


# compute a z score of the timing difference
m$timing_diff_rescaled_z <- m$timing_diff_rescaled %>% scale %>% as.vector
m$timing_diff_rescaled_z %>% density %>% plot



#compute a p value for the z score
m$timing_diff_rescaled_p <- 2 * pnorm(-abs(m$timing_diff_rescaled_z)) %>% round(6)

2 * pnorm(-abs(1.96)) %>% round(6)

# label the cna with a p value < 0.05
m$label <- ifelse(m$timing_diff_rescaled_p < 0.05, "signif", "ns")


m %>% ggplot(aes( cna_order ,timing_diff_rescaled_z)) + 
geom_bar(stat = "identity", aes(fill=timing_diff_rescaled), colour="black") +
# rotate labels x axis
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# change colour
scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
# plot significance lines
geom_hline(yintercept = c(1.96, -1.96), linetype = "dashed")




#____________________ VIZ 1 _________________________________

m$timing_diff <- m$timing_com - m$timing_smm

m2 <- m %>% filter(calls_com>30 & calls_smm > 10)

m2 %>% ggplot(aes(timing_com, timing_smm)) + 
geom_abline(intercept = 0, slope = 1, color = "grey50", linetype=2) +
geom_errorbar(aes(ymin= timing_smm - quasiSE_smm , ymax= timing_smm  + quasiSE_smm ), colour="grey50", size =0.3) +
geom_errorbar(aes(xmin = timing_com -  quasiSE_com, xmax= timing_com + quasiSE_com), colour="grey50", size =0.3) +
geom_point(aes(colour=timing_diff)) +
# label points more than 2 from the diagonal
geom_label_repel(data= m2 %>% filter(abs(timing_com - timing_smm) > 1.5), aes(label = cna), size = 3) +
geom_label_repel(data= m2 %>% filter(timing_com < 1.5 & timing_smm < 1.5), aes(label = cna), size = 3) +
# scale colour gradient mid point at 0, high point at 1 and low point at -1
scale_colour_gradient2(low = "blue", mid = "black", high = "red", 
midpoint = 0, limits=c(-2,2), oob= scales::squish ) +
ggpubr::stat_cor() +
xlab("MM timing estimate") + ylab("SMM timing estimate") +
theme_bw()


dir.create("results/JB_paper/def_aug23", showWarnings = FALSE)

# save plot
ggsave("results/JB_paper/def_aug23/compare_Viz1_SMM_vs_CoMMpass_timings.png", width = 8, height = 7, dpi = 300)


#____________________ VIZ 2 _________________________________

# plot the normal timing estimates for the two models with a scatterplot, transform the data in long format 

ml <- m2 %>% select(cna, timing_com, timing_smm, quasiSE_com, quasiSE_smm) %>%
pivot_longer(cols = c(timing_com, timing_smm), names_to = "model", values_to = "timing") %>% 
mutate(model = factor(model, levels = c("timing_com", "timing_smm")))


ml$quasiSE <- ifelse(ml$model == "timing_com", ml$quasiSE_com, ml$quasiSE_smm)

# create a test ggplot with categorical variable in x and continuos variable in y with associated error bars
# shift the points in the plot bwtween the two models



ml$cna_order <- factor(ml$cna, levels = m$cna[order(m$timing_com, decreasing = T)])

ml %>% ggplot(aes(cna_order, timing, colour = model)) +
    geom_errorbar(aes(ymin = timing - 2 * quasiSE, ymax = timing + 2 * quasiSE, group = model),
        colour = "grey50", size = 0.3, position = position_dodge(width = 0.8)) +
    geom_point(position = position_dodge(width = 0.8), size = 1) +
    geom_smooth(method = "lm") +
    theme_bw() +
    ylab("timing estimate") +
    scale_color_discrete(labels = c("MM", "SMM")) +
    coord_flip()

# save plot
ggsave("results/JB_paper/def_aug23/compare_Viz2_SMM_vs_CoMMpass_timings.png", width = 8, height = 8, dpi = 300)



#____________________ VIZ 3 _________________________________

m2$cna_order <- m2$cna %>% factor(levels = m2$cna[order(m2$timing_diff)])


# BAR PLOT
p1 <- m2 %>% ggplot(aes( cna_order ,abs(timing_diff))) + 
geom_bar(stat = "identity", aes(fill=timing_diff), colour="black") +
cowplot::theme_minimal_hgrid()+
theme(axis.text.x = element_blank(), 
axis.ticks.x = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_text(size = 10),
axis.text.y = element_text(size = 9),
legend.position = "none") + ylab("Timing difference") +
scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 

p1


# BAR PLOT
ml2 <- ml

ml2$cna_order <- factor(ml2$cna, levels = m2$cna[order(m2$timing_diff)])

ml2 <- ml2 %>% group_by(cna) %>% mutate(timing_diff=diff(timing))

p3 <- ml2 %>% ggplot(aes(cna_order, timing, fill = model, colour = timing_diff)) +
    geom_errorbar(aes(ymin = timing - 2 * quasiSE, ymax = timing + 2 * quasiSE),
        colour = "grey50", size = 0.3, position = position_dodge(width = 0.8)
    ) +
    geom_point(position = position_dodge(width = 0.8), aes(shape=model), size=1.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom"
    ) +
    scale_colour_gradient2(low = "red", mid = "gray50", high = "blue", 
        midpoint = 0, limits=c(-2,2), oob= scales::squish ) +
    scale_shape_manual(values = c(1, 19), labels = c("MM", "SMM")) +
    scale_fill_discrete(labels = c("MM", "SMM")) +
    xlab("Alteration") + ylab("Timing estimate")
    
p3


# create a grid with plots p1 and p3
library(gridExtra)

grid.arrange(p1, p3, ncol = 1, heights = c(1, 2))

GR <- grid.arrange(p1, p3, ncol = 1, heights = c(1, 2))
GR
ggsave("results/JB_paper/def_aug23/compare_Viz3_SMM_vs_CoMMpass_timings.png", plot = GR, width = 10, height = 7, dpi = 300)

