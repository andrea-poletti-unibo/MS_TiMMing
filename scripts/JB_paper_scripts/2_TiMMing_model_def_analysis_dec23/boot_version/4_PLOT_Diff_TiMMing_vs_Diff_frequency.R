library(tidyverse)
library(data.table)
library(ggrepel)
library(gridExtra)


#____________________ VIZ 5: diffs ______________________________

v5 <- fread("workfiles/data_diff-timing_vs_diff-frequency_SMM_CoMMpass.txt")
v5$type <- ifelse(v5$variable  %>% str_detect("Hyper|amp"), "Amplification", "Deletion")

v5$label <- v5$variable %>% 
    str_replace("del","-") %>% 
    str_replace("amp","+") %>% 
    str_replace("_chr_"," ") %>% 
    str_replace_all("_"," ") %>% 
    str_replace(" amp","")


# invert for better interpretation
v5$timing_diff_inv <- -v5$timing_diff


v5 %>% ggplot(aes(diff*100, timing_diff_inv)) + 
    geom_hline(yintercept = 0, linetype=1, alpha=0.4) + 
    geom_vline(xintercept = 0, linetype=1, alpha=0.4) + 
    geom_smooth(method = "lm", alpha=0.2, colour="grey50") + 
    ggpubr::stat_cor(method = "spearman", label.x = -9, label.y = 1.8) +
    geom_point(aes(colour=type)) + 
    # geom_abline(intercept = 0, slope = 10, color = "red") +
    geom_text_repel(
        aes(label = label), 
        size = 2.3, 
        max.overlaps=8)  + 
    theme_classic() + 
    xlab("Frequency increase from MGUS/SMM to NDMM (%)") + 
    ylab("Difference of BT-scores \n (NDMM minus MGUS/SMM)") +
    coord_cartesian(xlim=c(-9, 16)) +
    ylim(-2.6,2.6) +
    scale_x_continuous(n.breaks = 10) +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = c("#fb8072", "#80b1d3")) +
    labs(colour = "Copy number abnormality")

ggsave("results/JB_paper/def_aug23/boot_approach/TiMMing-diff_vs_freq-diff_paper_plot.svg" , width = 6, height = 6, dpi = 300)
ggsave("results/JB_paper/def_aug23/boot_approach/TiMMing-diff_vs_freq-diff_paper_plot.pdf" , width = 6, height = 6, units="in")




v5_strict_hifreq <- v5 %>% filter( NDMM > 0.05 & SMM > 0.05)


v5_strict_hifreq %>% ggplot(aes(diff, timing_diff_inv)) + 
    geom_hline(yintercept = 0, linetype=1, alpha=0.4) + 
    geom_vline(xintercept = 0, linetype=1, alpha=0.4) + 
    geom_smooth(method = "lm", alpha=0.2, colour="grey50") + 
    ggpubr::stat_cor(method = "spearman", label.x = -0.10, label.y = 1.8) +
    geom_point(aes(colour=type)) + 
    # geom_abline(intercept = 0, slope = 10, color = "red") +
    geom_text_repel(
        aes(label = label), 
        size = 3, 
        max.overlaps=100)  + 
    theme_classic() + 
    xlab("Frequency increase from SMM to NDMM") + 
    ylab("Timing difference (NDMM - SMM)") +
    xlim(-0.12, 0.16) + ylim(-1.8,2.2) +
    # legend bottom
    theme(legend.position = "bottom") 

ggsave("results/JB_paper/def_aug23/boot_approach/TiMMing-diff_vs_freq-diff_HIFREQ_paper_plot.png" , width = 6, height = 6, dpi = 300)
ggsave("results/JB_paper/def_aug23/boot_approach/TiMMing-diff_vs_freq-diff_HIFREQ_paper_plot.pdf" , width = 6, height = 6, units="in")
