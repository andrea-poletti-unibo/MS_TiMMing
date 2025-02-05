library(tidyverse)
library(data.table)

#============ Bustoros ============

bus_CN <- fread("results/JB_Paper/calls/clonalityCalls_broad_focal_BUSnew.txt")
bus_CI <- fread("results/JB_Paper/calls/allCalls_CN_CI_broad_focal_BUSnew.txt")

bus_samples <- bus_CN %>%
  select(sample) %>%
  distinct() %>%
  pull()


# because of doublets with WGS
remove.samples <- paste0("SMM_", c("017", "040"), "_Tumor")
# beacuse of no myeloma / IGH mutation or cnv in actual genome data (only relies on FISH from clinic?)
remove.samples <- c(remove.samples, paste0("SMM_", str_pad( c(9, 29, 37, 38, 42), 3,side = "left", pad= "0"), "_Tumor"))
# because of non human contamination
remove.samples <- c(remove.samples, "SMM_074_Tumor")

# additional doublets discovered:
#
# 'SMM-076-Tumor' / SMM102 from British Cohort == MBp01
# 'SMM-081-Tumor' / SMM108 from British Cohort == MBp09
# 'SMM-093-Tumor' / SMM127 from British Cohort == MBp30

remove.samples <- c(remove.samples, paste0("SMM_", c("076", "081", "093"), "_Tumor"))

remove.samples <- str_replace(remove.samples, "SMM_", "SMM-") %>% str_replace("_Tumor", "_pair")

remove.samples %in% bus_samples %>% table


bus_samples_keep <- bus_samples[!bus_samples %in% remove.samples]


bus_CN_keep <- bus_CN %>%
  filter(sample %in% bus_samples_keep)

bus_CI_keep <- bus_CI %>%
  filter(ID %in% bus_samples_keep)




#============ SU2C ============

suc_CN <- fread("results/JB_Paper/calls/clonalityCalls_broad_focal_SU2Cnew.txt")

suc_CI <- fread("results/JB_Paper/calls/allCalls_CN_CI_broad_focal_SU2Cnew.txt")

suc_samples <- suc_CN %>%
  select(sample) %>%
  distinct() %>%
  pull()

franG <- c(
"PD47561a",
"PD47563a",
"PD47567a",
"PD47570a",
"PD47572a",
"PD47573a",
"PD47574a",
"PD47575a",
"PD47576a",
"PD47577a",
"PD47578a",
"PD47579a",
"PD47580a",
"PD47581a",
"PD47582a")

franG %in% suc_samples %>% table

w <- c("IID_H196061_T01",
"IID_H196062_T01",
"IID_H196063_T01",
"IID_H196064_T01")

w %in% suc_samples %>% table


annot_aug23 <- readxl::read_xlsx("data/JB_Paper/final_Aug23/Annotations_SW.xlsx", sheet = 1)

suc_samples %in% annot_aug23$Reference_Pair

suc_samples %in% annot_aug23$Reference_Pair %>% table

# which samples are false

suc_samples[!suc_samples %in% annot_aug23$Reference_Pair]

annot_aug23$Reference_Pair %in%  suc_samples


# which samples are missing

annot_aug23$Reference_Pair[!annot_aug23$Reference_Pair %in% suc_samples]


#============ CoMMpass ============

com_portal <- fread("data/JB_Paper/final_Aug23/mmrf_final_annot.tsv")


com_samples <- com_portal %>%
  select(SPECTRUM_SEQ) %>% 
    distinct() %>% 
      pull()

com_CN <- fread("results/JB_Paper/calls/clonalityCalls_broad_focal_CoMM-SU2Cnew.txt")

com_CI <- fread("results/JB_Paper/calls/allCalls_CN_CI_broad_focal_CoMM-SU2Cnew.txt")

s <- com_CN$sample %>% unique() %>% str_remove("_BM_.*")


s %in% com_samples %>% table

com_samples %in% s %>% table

# which are false

com_samples[!com_samples %in% s]


com_CN$SPECTRUM_SEQ <- com_CN$sample %>% str_remove("_BM_.*")
com_CI$SPECTRUM_SEQ <- com_CI$ID %>% str_remove("_BM_.*")


com_CN_keep <- com_CN %>%
  filter(SPECTRUM_SEQ %in% com_samples)

com_CI_keep <- com_CI %>%
  filter(SPECTRUM_SEQ %in% com_samples)



######### export ###########

dir.create("results/JB_paper/calls/CN_only", showWarnings = FALSE, recursive = TRUE)

write_tsv(com_CN_keep, "results/JB_paper/calls/CN_only/CoMMpass_CN_762_samples.txt")

write_tsv(suc_CN, "results/JB_paper/calls/CN_only/SU2C_CN_164_samples.txt")

write_tsv(bus_CN_keep, "results/JB_paper/calls/CN_only/Bustoros_CN_164_samples.txt")
