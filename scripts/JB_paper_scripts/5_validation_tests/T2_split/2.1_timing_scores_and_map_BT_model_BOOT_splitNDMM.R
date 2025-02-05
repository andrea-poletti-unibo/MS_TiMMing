library(data.table)
library(tidyverse)
library(BradleyTerry2)
library(qvcalc)

# set the output directory
outdir <- "results/JB_paper/def_aug23/validation_tests/T2_split_NDMM/"
dir.create(outdir, showWarnings = F, recursive = T)


################################## split A #####################################

workpath <- "workfiles/validation_tests/T2_split_NDMM/split_NDMM-a/"


# read in the data
files <- list.files(workpath, full.names = T, pattern = "matches")

f=files[1]

res_all <- data.frame()

for (f in files) {

  boot <- f %>%
    str_extract("BOOT_[0-9]+") %>%
    str_remove("BOOT_") %>%
    as.numeric()
  message("boot: ", boot)

  import <- fread(f)

  import2 <- import %>% select(CNA1 = player1, CNA2 = player2, win1, win2)

  # round scores
  import2$win1 <- import2$win1 %>% round()
  import2$win2 <- import2$win2 %>% round()

  matchesB <- import2 %>% as.data.frame()

# exclude <- c("amp_chr_5q", "amp_chr_7p")

  exclude <- c("")

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

  ab_df$boot <- boot

  res_all <- rbind(res_all, ab_df)

  # gg <- ab_df %>%
  #   ggplot(aes(x = cna, y = timing)) +
  #   geom_errorbar(aes(ymin = timing - quasiSE, ymax = timing + quasiSE)) +
  #   geom_point(aes(color = CNA_type)) +
  #   xlab("CNA") +
  #   ylab("Timing Estimate") +
  #   ggtitle("Bradley-Terry CNAs timing estimates on CoMMpass new cohort", subtitle = paste("boot n = ", boot)) +
  #   coord_flip() +
  #   theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))

}


write_tsv(res_all, paste0(outdir, "/splitNDMM-a_timings_20boot.txt"))
res_all <- fread(paste0(outdir, "splitNDMM-a_timings_20boot.txt"))


res_all_boots <- res_all %>% 
  group_by(cna) %>% 
  summarise(
    med_timing=median(timing), 
    mean_timing=mean(timing), 
    sd = sd(timing), 
    min= min(timing), 
    max= max(timing), 
    first_q= quantile(timing)[2], 
    third_q= quantile(timing)[4],  
    n=n()) %>% 
  arrange(med_timing)

exclude <- res_all_boots %>% filter(n < 15) %>% pull(cna)

res_all_boots <- res_all_boots %>% filter(!cna %in% exclude)
res_all <- res_all %>% filter(!cna %in% exclude)


write_tsv(res_all, paste0(workpath, "/BTmodels_20BOOTS_splitNDMM-a.txt"))

write_tsv(res_all_boots, paste0(workpath, "/Aggregate_BTmodels_20BOOTS_splitNDMM-a.txt"))


ORD <- res_all_boots$med_timing %>% order(decreasing = T)

res_all_boots$CNA_type <- ifelse(grepl("amp|Hyp", res_all_boots$cna), "Amp", "Del")
res_all_boots$cna_f <- factor(res_all_boots$cna, levels = res_all_boots$cna[ORD], ordered = T)

res_all$CNA_type <- ifelse(grepl("amp|Hyp", res_all$cna), "Amp", "Del")
res_all$cna_f <- factor(res_all$cna, levels = res_all_boots$cna[ORD],  ordered = T)


p1 <- res_all %>% ggplot(aes(x=timing, y = cna_f)) + 
  geom_jitter(size=0.5, alpha=0.5, height = 0.2, aes(colour=CNA_type)) + 
  # geom_boxplot(outlier.shape = NA, fill=NA) +
  geom_violin(fill=NA, scale = "width", trim = T, adjust = .5, colour="grey70") +
  geom_point(data = res_all_boots, aes(x=med_timing, y = cna_f, colour=CNA_type), size=2)+
  geom_errorbar(data = res_all_boots, aes(x=med_timing, y = cna_f, xmin= (mean_timing-(sd)), xmax=(mean_timing+(sd))), width=0.2) +
  theme_bw() +
  theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm")) +
  xlab("Timing Estimate") +
  ylab("CNA") +
  ggtitle("Bradley-Terry CNAs timing estimates on split-NDMM cohort", subtitle = "20 bootstraps runs")

p1

dir.create("plots/timing_maps/validation_tests/T2_split/", showWarnings = F, recursive = T)
# save plot
ggsave(plot = p1, filename =  "plots/timing_maps/validation_tests/T2_split/splitNDMM-a_timing_map_20BOOTS.png", width = 14, height = 12, units = "in")





################################## split B #####################################

workpath <- "workfiles/validation_tests/T2_split_NDMM/split_NDMM-b/"


# read in the data
files <- list.files(workpath, full.names = T, pattern = "matches")

f=files[1]

res_all <- data.frame()

for (f in files) {

  boot <- f %>%
    str_extract("BOOT_[0-9]+") %>%
    str_remove("BOOT_") %>%
    as.numeric()
  message("boot: ", boot)

  import <- fread(f)

  import2 <- import %>% select(CNA1 = player1, CNA2 = player2, win1, win2)

  # round scores
  import2$win1 <- import2$win1 %>% round()
  import2$win2 <- import2$win2 %>% round()

  matchesB <- import2 %>% as.data.frame()

# exclude <- c("amp_chr_5q", "amp_chr_7p")

  exclude <- c("")

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

  ab_df$boot <- boot

  res_all <- rbind(res_all, ab_df)

  # gg <- ab_df %>%
  #   ggplot(aes(x = cna, y = timing)) +
  #   geom_errorbar(aes(ymin = timing - quasiSE, ymax = timing + quasiSE)) +
  #   geom_point(aes(color = CNA_type)) +
  #   xlab("CNA") +
  #   ylab("Timing Estimate") +
  #   ggtitle("Bradley-Terry CNAs timing estimates on CoMMpass new cohort", subtitle = paste("boot n = ", boot)) +
  #   coord_flip() +
  #   theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))

}


write_tsv(res_all, paste0(outdir, "/splitNDMM-b_timings_20boot.txt"))
res_all <- fread(paste0(outdir, "splitNDMM-b_timings_20boot.txt"))


res_all_boots <- res_all %>% 
  group_by(cna) %>% 
  summarise(
    med_timing=median(timing), 
    mean_timing=mean(timing), 
    sd = sd(timing), 
    min= min(timing), 
    max= max(timing), 
    first_q= quantile(timing)[2], 
    third_q= quantile(timing)[4],  
    n=n()) %>% 
  arrange(med_timing)


exclude <- res_all_boots %>% filter(n < 15) %>% pull(cna)

res_all_boots <- res_all_boots %>% filter(!cna %in% exclude)
res_all <- res_all %>% filter(!cna %in% exclude)

write_tsv(res_all, paste0(workpath, "/BTmodels_20BOOTS_splitNDMM-b.txt"))

write_tsv(res_all_boots, paste0(workpath, "/Aggregate_BTmodels_20BOOTS_splitNDMM-b.txt"))


ORD <- res_all_boots$med_timing %>% order(decreasing = T)

res_all_boots$CNA_type <- ifelse(grepl("amp|Hyp", res_all_boots$cna), "Amp", "Del")
res_all_boots$cna_f <- factor(res_all_boots$cna, levels = res_all_boots$cna[ORD], ordered = T)

res_all$CNA_type <- ifelse(grepl("amp|Hyp", res_all$cna), "Amp", "Del")
res_all$cna_f <- factor(res_all$cna, levels = res_all_boots$cna[ORD],  ordered = T)



p1 <- res_all %>% ggplot(aes(x=timing, y = cna_f)) + 
  geom_jitter(size=0.5, alpha=0.5, height = 0.2, aes(colour=CNA_type)) + 
  # geom_boxplot(outlier.shape = NA, fill=NA) +
  geom_violin(fill=NA, scale = "width", trim = T, adjust = .5, colour="grey70") +
  geom_point(data = res_all_boots, aes(x=med_timing, y = cna_f, colour=CNA_type), size=2)+
  geom_errorbar(data = res_all_boots, aes(x=med_timing, y = cna_f, xmin= (mean_timing-(sd)), xmax=(mean_timing+(sd))), width=0.2) +
  theme_bw() +
  theme(legend.position = "bottom", plot.margin = unit(c(0.5,0,0,0.5), "cm")) +
  xlab("Timing Estimate") +
  ylab("CNA") +
  ggtitle("Bradley-Terry CNAs timing estimates on split-NDMM-b cohort", subtitle = "20 bootstraps runs")

p1

dir.create("plots/timing_maps/validation_tests/T2_split/", showWarnings = F, recursive = T)
# save plot
ggsave(plot = p1, filename =  "plots/timing_maps/validation_tests/T2_split/splitNDMM-b_timing_map_20BOOTS.png", width = 14, height = 12, units = "in")
