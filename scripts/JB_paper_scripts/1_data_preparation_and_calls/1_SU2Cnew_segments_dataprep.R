library(tidyverse)
library(data.table)
library(readxl)

seg <- fread("data/JB_Paper/Gistic_newSU2C_all_cohorts_180samp_v2.seg")

seg <- BioNerds::format_segments(seg) %>% dplyr::rename(ID = Sample)

# data is in log2R form
seg$Seg.CN %>% density(na.rm=T) %>% plot()

# compute width
seg$width <- seg$end - seg$start

# filter out NA values
seg$Seg.CN %>% is.na %>% table

seg <- seg %>% filter(!is.na(Seg.CN))

# compute CN from log2r
seg$CN <- 2^(seg$Seg.CN + 1)

# visual check
viz <- seg
viz$CN_viz <- viz$CN
viz$CN_viz[viz$CN_viz>5] <- 5
viz$CN_viz[viz$CN_viz< 0] <- 0

viz %>% ggplot(aes(CN_viz, weight = width)) + 
  geom_density() +
  geom_vline(xintercept = 1, colour="blue") + 
  geom_vline(xintercept = 3, colour="red")


samples <- seg$ID %>% unique

anno <- read_xlsx("data/JB_Paper/DuttaAlberge_CleanClinicalRef_noMRN_29-06-23.xlsx")

anno$Disease_Status %>% table

MM_samples <- anno %>% filter(Disease_Status %in% c("MM","NDMM")) %>% dplyr::select(ID= Reference_Pair) %>% pull(ID)

# filter out MM samples from seg

seg2 <- seg %>% filter(!ID %in% MM_samples)

# from 180 to 164 after filtering out MM samples !
seg2$ID %>% unique %>% length

pur <- anno %>% dplyr::select(ID= Reference_Pair, purity)

pur$ID %in% samples %>% table
samples %in% pur$ID %>% table

# data.frame(samples, check = samples %in% pur$ID) %>% View


seg3 <- left_join(seg2, pur, by="ID")


#========================== ComphyNumber =======================================

# estimate SD from a MAD of 0.1

MAD <- 0.1

SD <- MAD / (sqrt(2/pi))

seg3$SD <- SD

# 95% Confidence Interval FORMULA
# CI = 1.959*(SD/sqrt(n_probes))

seg3$CI95 <-1.959*(seg3$SD / sqrt(seg3$n_probes) ) 

# correct CI per purity

seg3$CI95_purCorr <- seg3$CI95 / seg3$purity

seg3$CI95_purCorr %>% summary
seg3$CI95 %>% summary

write_tsv(seg3, "data/JB_Paper/dataprep_segments_SU2CNew_164.txt")