library(tidyverse)
library(data.table)

ndt <- fread("data/Phylogic_NDT_League_Model/CoMMpass_phylogic_league_model_JB_250124/AllCoMMpass.log_odds.tsv")
setnames(ndt, c("cohort", "event", "event_split", "type", "perm_run", "log_odds_early","NA"))

s_ndt <- ndt %>% group_by(event) %>% summarise(
    mean_log_odds= mean(type),
    median_log_odds = median(type),
    sd_log_odds = sd(type),
    n = n()
)

s_ndt$mean_log_odds_rev <- -s_ndt$mean_log_odds


s_ndt %>% ggplot(aes(x=mean_log_odds, y=event)) + geom_point()
# sort by mean_log_odds
s_ndt %>% arrange(mean_log_odds) %>% ggplot(aes(x=mean_log_odds, y=event)) + geom_point()

# sort in ggplot axis
s_ndt %>% ggplot(aes(x=mean_log_odds, y=reorder(event, mean_log_odds, decreasing=T))) + geom_point()

s_ndt %>% ggplot(aes(x=mean_log_odds_rev, y=reorder(event, mean_log_odds_rev, decreasing=T))) + geom_point()

ggsave("plots/Phylogic_NDT_comparison/CoMMpass_PhylogicNDT_log_odds.png", width=9, height=8, units="in")



######### load TiMMing data #########

tim <- fread("results/JB_paper/def_aug23/CoMMpass_20boots/CoMMpass_timings_20boot.txt")

s_tim <- tim %>% group_by(cna) %>% summarise(
    mean_timing= mean(timing),
    median_timing = median(timing),
    sd_timing = sd(timing),
    n = n()
)

s_tim$mean_timing_rev <- -s_tim$mean_timing

# sort by mean_timing
s_tim %>% ggplot(aes(x=mean_timing, y=reorder(cna, mean_timing, decreasing=T))) + geom_point()

ggsave("plots/Phylogic_NDT_comparison/CoMMpass_TiMMing_scores.png", width=9, height=12, units="in")

#============= rename the CNA events in common format =============

s_tim$cna 
s_ndt$event

s_tim$event <- s_tim$cna %>% str_replace_all("_chr", "") %>% 
    str_replace_all("amp", "gain") %>% 
    str_replace_all("del", "loss") %>%
    str_extract("gain_[0-9]+[pq]|loss_[0-9]+[pq]")

s_tim$event 


merge <- left_join(s_tim, s_ndt, by = c("event" = "event"))


dir.create("plots/Phylogic_NDT_comparison/", showWarnings = F, recursive = T)



merge %>% filter(mean_timing<4.5) %>%
    ggplot(aes(x=mean_timing, y=mean_log_odds_rev)) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    ggpubr::stat_cor(method = "spearman") +
    ggrepel::geom_text_repel(aes(label = event)) +
    ggtitle("TiMMing scores vs PhylogicNDT log odds ratios")

ggsave("plots/Phylogic_NDT_comparison/CoMMpass_TiMMing_vs_PhylogicNDT_all_events.png", width=9, height=8, units="in")



merge %>% filter(mean_timing<4.5) %>%  
    ggplot(aes(x=mean_timing, y=mean_log_odds_rev)) + geom_point() + 
    geom_smooth(method = "lm") +
    ggpubr::stat_cor(method = "spearman") +
    ggrepel::geom_text_repel(data = merge %>% filter(
        (mean_log_odds_rev > mean_timing - 1.5) |
        (mean_log_odds_rev < mean_timing - 3.5)),
        # (mean_timing < 2 & mean_log_odds_rev > 0)
        aes(label = event)) +
    ggrepel::geom_text_repel(data = merge %>% filter(
        (mean_log_odds_rev < mean_timing - 1.5) &
        (mean_log_odds_rev > mean_timing - 3.5)),
        # (mean_timing < 2 & mean_log_odds_rev > 0)
        aes(label = event), color="blue") +
    geom_abline(intercept = -2.5, slope = 1, color="red") +
    geom_abline(intercept = -3.5, slope = 1, color="red", linetype="dashed") +
    geom_abline(intercept = -1.5, slope = 1, color="red", linetype="dashed") +
    ggtitle("TiMMing scores vs PhylogicNDT log odds ratios - interpretation")


ggsave("plots/Phylogic_NDT_comparison/CoMMpass_TiMMing_vs_PhylogicNDT_all_events_interpretation.png", width=9, height=8, units="in")


merge %>% filter(
        (mean_log_odds_rev < mean_timing - 1.5) &
        (mean_log_odds_rev > mean_timing - 3.5)) %>%  
    ggplot(aes(x=mean_timing, y=mean_log_odds_rev)) + geom_point() + 
    geom_smooth(method = "lm") +
    ggrepel::geom_text_repel(aes(label = event)) +
    ggpubr::stat_cor(method = "spearman") +
    ggtitle("TiMMing scores vs PhylogicNDT log odds ratios - outliers excluded") 


ggsave("plots/Phylogic_NDT_comparison/CoMMpass_TiMMing_vs_PhylogicNDT_outliers_excluded.png", width=9, height=8, units="in")

