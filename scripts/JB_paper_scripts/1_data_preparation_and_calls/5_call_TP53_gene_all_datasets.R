################### CN CALL TP53 - all datasets ###################

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(BioNerds)


# load SU2C data prep segments
segm_SU2C <- fread("data/JB_Paper/dataprep_segments_SU2CNew_164.txt")
segm_SU2C_f <- segm_SU2C %>% select(chr, start, end, width, ID, CN, CI95_purCorr)

# load BUS data prep segments
segm_BUS_TN <- fread("data/Bustoros_WES/dataprep_segments_Bus-TN_83.seg")
segm_BUS_TO <- fread("data/Bustoros_WES/dataprep_segments_Bus-TO_92.seg")
segm_BUS_TN$added <- F
segm_BUS <- rbind(segm_BUS_TN, segm_BUS_TO) # merge TN and TO

segm_BUS_f <- segm_BUS %>% select(chr, start, end, width, ID, CN, CI95_purCorr = seg_sigma)

# load CoMM data prep segments
segm_CoMM <- fread("data/CoMMpass_acs/dataprep_segments_CoMM_847.seg")

segm_CoMM_f <- segm_CoMM %>% select(chr, start, end, width, ID, CN = BOB_CN, CI95_purCorr = seg_sigma)

# merge all segments


names(segm_SU2C)
names(segm_BUS)
names(segm_CoMM)


# define TP53 gene locus manually # coordinates from https://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53 

all_focal_df <- data.frame(
  hgnc_symbol = "TP53",
  chromosome_name = 17,
  start_position =7661779,
  end_position = 7687538,
  type="gene_del",
  refgenome="hg38")




#======================= perform Focal Calls SU2C ========================
segm <- segm_SU2C_f

# create samples vector and SORT IT (CRITICAL!)
samples <- segm$ID %>% unique %>% sort


# create a Granges object from segments
Gsegm <- makeGRangesFromDataFrame(segm, keep.extra.columns = T)

# create a Granges object from the genes
Gfoc_genes <- all_focal_df %>% 
  dplyr::rename(chr=chromosome_name, 
                start=start_position, 
                end=end_position) %>%  
  makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T) 



#__________ Gene-level calls algorithm _____________

all_focal_df$hgnc_symbol[1]
i=1

genesCallsRes <- data.frame()

for ( i in 1:nrow(all_focal_df)) {
  G_gene <- Gfoc_genes[i] # select a single gene
  message(G_gene$hgnc_symbol)
  wid_gene <- G_gene %>% width() # save the length of the selected gene
  overlap <- pintersect( Gsegm, G_gene) # comupte overlaps between the gene and the segment-gene positions
  res.gene <- overlap %>% as.data.frame() %>% filter(hit==TRUE) # generate the GENE-level results
  res.gene$percLen <- res.gene$width/wid_gene %>% round(3) # create the percentage of segment overlap with the gene
  res.gene$sample %>% unique() %>% length()
  res.gene <- dplyr::filter(res.gene, percLen >= 0.10) # exclude those segments that overlap the gene for less then 10%
  res.gene <- res.gene %>% arrange(ID, desc(abs(2- CN ))) # ! SUPER IMPORTANT ORDERING: per sample and then per descending CN change -> the first entry per sample will be the most changed in CN!
  res.gene.extreme <- dplyr::distinct(res.gene, ID, .keep_all = TRUE ) # KEEP only the first entry for each sample (this will end up always with a number of rows equal to number of samples) 
  CNvalues<- res.gene.extreme$CN %>% round(4)# extract the CN values 
  CI95values<- res.gene.extreme$CI95_purCorr %>% round(4) # extract the CI95 values corrected for purity
  dataExport <- paste0(CNvalues, "|", CI95values)
  dfExport <- res.gene.extreme %>% mutate(gene = G_gene$hgnc_symbol, 
                                          CN_CI= paste0(CN %>% round(4),"|", CI95_purCorr %>% round(4) )) %>% 
    dplyr::select(ID, gene, CN_CI)
  genesCallsRes <- rbind(genesCallsRes, dfExport)
}

genesCallsFinal <- reshape2::dcast(genesCallsRes, ID ~ gene )

# check number of NA
genesCallsFinal %>% reshape2::melt(id.vars = "ID") %>% 
  group_by(ID) %>% summarise(nNA=sum(is.na(value))) %>% 
  arrange(desc(nNA)) %>% View

# eventually fill NAs
genesCallsFinal[is.na(genesCallsFinal)] %>% length()
genesCallsFinal[is.na(genesCallsFinal)] <- "2.009|0.009"
genesCallsFinal[is.na(genesCallsFinal)] %>% length()


TP53_SU2C <- genesCallsFinal %>% dplyr::select(ID, TP53_CN = TP53) %>% 
  dplyr::mutate(dataset="SU2C", TP53_CN=TP53_CN %>% str_remove("\\|.*") %>% as.numeric() )





#======================= perform Focal Calls BUS ========================

segm <- segm_BUS_f

# create samples vector and SORT IT (CRITICAL!)
samples <- segm$ID %>% unique %>% sort

# create a Granges object from segments
Gsegm <- makeGRangesFromDataFrame(segm, keep.extra.columns = T)

# create a Granges object from the genes
Gfoc_genes <- all_focal_df %>% 
  dplyr::rename(chr=chromosome_name, 
                start=start_position, 
                end=end_position) %>%  
  makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T) 



#__________ Gene-level calls algorithm _____________

all_focal_df$hgnc_symbol[1]
i=1

genesCallsRes <- data.frame()

for ( i in 1:nrow(all_focal_df)) {
  
  G_gene <- Gfoc_genes[i] # select a single gene
  message(G_gene$hgnc_symbol)
  wid_gene <- G_gene %>% width() # save the length of the selected gene
  overlap <- pintersect( Gsegm, G_gene) # comupte overlaps between the gene and the segment-gene positions
  res.gene <- overlap %>% as.data.frame() %>% filter(hit==TRUE) # generate the GENE-level results
  res.gene$percLen <- res.gene$width/wid_gene %>% round(3) # create the percentage of segment overlap with the gene
  res.gene$sample %>% unique() %>% length()
  res.gene <- dplyr::filter(res.gene, percLen >= 0.10) # exclude those segments that overlap the gene for less then 10%
  res.gene <- res.gene %>% arrange(ID, desc(abs(2- CN ))) # ! SUPER IMPORTANT ORDERING: per sample and then per descending CN change -> the first entry per sample will be the most changed in CN!
  res.gene.extreme <- dplyr::distinct(res.gene, ID, .keep_all = TRUE ) # KEEP only the first entry for each sample (this will end up always with a number of rows equal to number of samples) 
  CNvalues<- res.gene.extreme$CN %>% round(4)# extract the CN values 
  CI95values<- res.gene.extreme$CI95_purCorr %>% round(4) # extract the CI95 values corrected for purity
  dataExport <- paste0(CNvalues, "|", CI95values)
  dfExport <- res.gene.extreme %>% mutate(gene = G_gene$hgnc_symbol, 
                                          CN_CI= paste0(CN %>% round(4),"|", CI95_purCorr %>% round(4) )) %>% 
    dplyr::select(ID, gene, CN_CI)
  genesCallsRes <- rbind(genesCallsRes, dfExport)
}

genesCallsFinal <- reshape2::dcast(genesCallsRes, ID ~ gene )

# check number of NA
genesCallsFinal %>% reshape2::melt(id.vars = "ID") %>% 
  group_by(ID) %>% summarise(nNA=sum(is.na(value))) %>% 
  arrange(desc(nNA)) %>% View

# eventually fill NAs
genesCallsFinal[is.na(genesCallsFinal)] %>% length()
genesCallsFinal[is.na(genesCallsFinal)] <- "2.009|0.009"
genesCallsFinal[is.na(genesCallsFinal)] %>% length()


TP53_BUS <- genesCallsFinal %>% dplyr::select(ID, TP53_CN = TP53) %>% 
  dplyr::mutate(dataset="BUS", TP53_CN=TP53_CN %>% str_remove("\\|.*") %>% as.numeric() )


#======================= perform Focal Calls CoMM ========================

segm <- segm_CoMM_f


# create samples vector and SORT IT (CRITICAL!)
samples <- segm$ID %>% unique %>% sort

# create a Granges object from segments
Gsegm <- makeGRangesFromDataFrame(segm, keep.extra.columns = T)

# create a Granges object from the genes
Gfoc_genes <- all_focal_df %>% 
  dplyr::rename(chr=chromosome_name, 
                start=start_position, 
                end=end_position) %>%  
  makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T) 



#__________ Gene-level calls algorithm _____________

all_focal_df$hgnc_symbol[1]
i=1

genesCallsRes <- data.frame()

for ( i in 1:nrow(all_focal_df)) {
  G_gene <- Gfoc_genes[i] # select a single gene
  message(G_gene$hgnc_symbol)
  wid_gene <- G_gene %>% width() # save the length of the selected gene
  overlap <- pintersect( Gsegm, G_gene) # comupte overlaps between the gene and the segment-gene positions
  res.gene <- overlap %>% as.data.frame() %>% filter(hit==TRUE) # generate the GENE-level results
  res.gene$percLen <- res.gene$width/wid_gene %>% round(3) # create the percentage of segment overlap with the gene
  res.gene$sample %>% unique() %>% length()
  res.gene <- dplyr::filter(res.gene, percLen >= 0.10) # exclude those segments that overlap the gene for less then 10%
  res.gene <- res.gene %>% arrange(ID, desc(abs(2- CN ))) # ! SUPER IMPORTANT ORDERING: per sample and then per descending CN change -> the first entry per sample will be the most changed in CN!
  res.gene.extreme <- dplyr::distinct(res.gene, ID, .keep_all = TRUE ) # KEEP only the first entry for each sample (this will end up always with a number of rows equal to number of samples) 
  CNvalues<- res.gene.extreme$CN %>% round(4)# extract the CN values 
  CI95values<- res.gene.extreme$CI95_purCorr %>% round(4) # extract the CI95 values corrected for purity
  dataExport <- paste0(CNvalues, "|", CI95values)
  dfExport <- res.gene.extreme %>% mutate(gene = G_gene$hgnc_symbol, 
                                          CN_CI= paste0(CN %>% round(4),"|", CI95_purCorr %>% round(4) )) %>% 
    dplyr::select(ID, gene, CN_CI)
  genesCallsRes <- rbind(genesCallsRes, dfExport)
}

genesCallsFinal <- reshape2::dcast(genesCallsRes, ID ~ gene )

# check number of NA
genesCallsFinal %>% reshape2::melt(id.vars = "ID") %>% 
  group_by(ID) %>% summarise(nNA=sum(is.na(value))) %>% 
  arrange(desc(nNA)) %>% View

# eventually fill NAs
genesCallsFinal[is.na(genesCallsFinal)] %>% length()
genesCallsFinal[is.na(genesCallsFinal)] <- "2.009|0.009"
genesCallsFinal[is.na(genesCallsFinal)] %>% length()


TP53_CoMM<- genesCallsFinal %>% dplyr::select(ID, TP53_CN = TP53) %>% 
  dplyr::mutate(dataset="CoMM", TP53_CN=TP53_CN %>% str_remove("\\|.*") %>% as.numeric() )


#======================= merge all datasets ========================

all <- rbind(TP53_SU2C, TP53_BUS, TP53_CoMM)

all %>% write_tsv("results/JB_paper/calls/CN_only/CN_TP53_all_datasets.txt")



########################### compare with broad del 17p calls #################################


library(data.table)
library(tidyverse)
library(cowplot)
library(grid)
library(flextable)


foc <- fread("results/JB_paper/calls/CN_only/CN_TP53_all_datasets.txt")

foc$dataset %>% table

results/JB_paper/calls/allCalls_CN_CI_broad_focal_BUSnew.txt

br_suc <- fread("results/JB_paper/calls/allCalls_CN_CI_broad_focal_SU2Cnew.txt")
br_bus <- fread("results/JB_paper/calls/allCalls_CN_CI_broad_focal_BUSnew.txt")
br_com <- fread("results/JB_paper/calls/allCalls_CN_CI_broad_focal_CoMM-SU2Cnew.txt")

br <- rbind(br_com, br_suc, br_bus, fill=T)

merge <- left_join(foc, br, by=c("ID"))


merge$del_TP53 <- ifelse(merge$TP53_CN < 2, 2-merge$TP53_CN  ,0)

merge$del_TP53 %>% sort %>% plot


merge$del_chr_17p

merge$chr_17p_CN <- merge$chr_17p %>% str_remove("\\|.*") %>% as.numeric()

merge$chr_17p_CN %>% sort %>% plot


merge$del_chr_17p <- ifelse(merge$chr_17p_CN < 2, 2-merge$chr_17p_CN  ,0)

merge$del_chr_17p %>% sort %>% plot


merge %>% ggplot(aes(del_chr_17p, del_TP53, color=dataset)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_vline(xintercept = 0.1, linetype=1) +
  geom_hline(yintercept = 0.1, linetype=1) +
  theme_minimal() + 
  facet_wrap(~dataset, scales="free") +
  theme(legend.position = "bottom") + 
  labs(x="del chr17p", y="del TP53")

m2 <- merge %>% 
  mutate( concordance= case_when(
    del_TP53 >= 0.1 & del_chr_17p >= 0.1 ~ "concordant_pos",
    del_TP53 < 0.1 & del_chr_17p < 0.1 ~ "concordant_neg",
    del_TP53 >= 0.1 & del_chr_17p < 0.1 ~ "only_TP53",
    del_TP53 < 0.1 & del_chr_17p >= 0.1 ~ "only_chr17p"

  ))

m2$concordance %>% table

install.packages("crosstable")

library(crosstable)
ft <- crosstable(m2, c(concordance), by = "dataset", percent_pattern = "{n} ({p_col})") %>% 
  as_flextable() 

ft

ft_raster <- as_raster(ft)

mgg <- merge %>% ggplot(aes(del_chr_17p, del_TP53, color=dataset)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_vline(xintercept = 0.1, linetype=1) +
  geom_hline(yintercept = 0.1, linetype=1) +
  theme_minimal() + 
  facet_wrap(~dataset, scales="fixed") +
  theme(legend.position = "bottom") + 
  labs(x="del chr 17p (CN units)", y="del focal TP53 (CN units)") 

gflextable <- ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ft_raster), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)


plot_grid(mgg, gflextable, nrow = 2, ncol = 1, rel_heights = c(4, 1) )

ggsave("plots/del_TP53_vs_del_chr17p_concordance.pdf", width=8, height=5.5, dpi=300)
