################### BROAD AND FOCAL CALLS - Bustoros dataset ###################

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(BioNerds)
library(biomaRt)

#=================== load prepared data ====================

# load data prep segments
segm_TN <- fread("data/Bustoros_WES/dataprep_segments_Bus-TN_83.seg")
segm_TO <- fread("data/Bustoros_WES/dataprep_segments_Bus-TO_92.seg")


# merge the two sub-cohorts (tumor-normal pairs (TN) and Tumor-only (TO))
segm_TN$added <- F
segm <- rbind(segm_TN, segm_TO)



# visual check
viz <- segm
viz$CN[viz$CN>6] <- 6

viz %>% ggplot(aes(CN, weight = width)) + 
  geom_density() + xlim(0,4.5) + 
  geom_vline(xintercept = 1, colour="blue") + 
  geom_vline(xintercept = 3, colour="red")


# load focal genes (hg38) with coordinates from BioMart
all_focal_df <- fread("data/focal_loci_hg38_JBpaper.txt")



#======================= perform Focal Calls ========================

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

all_focal_df$hgnc_symbol[3]
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
  
  CI95values<- res.gene.extreme$seg_sigma %>% round(4) # extract the CI95 values corrected for purity
  
  dataExport <- paste0(CNvalues, "|", CI95values)
  
  
  dfExport <- res.gene.extreme %>% mutate(gene = G_gene$hgnc_symbol, 
                                          CN_CI= paste0(CN %>% round(4),"|", seg_sigma %>% round(4) )) %>% 
    dplyr::select(ID, gene, CN_CI)
  
  genesCallsRes <- rbind(genesCallsRes, dfExport)
  
}

table(genesCallsRes$ID, genesCallsRes$gene)



genesCallsFinal <- reshape2::dcast(genesCallsRes, ID ~ gene )

# check number of NA
genesCallsFinal %>% reshape2::melt(id.vars = "ID") %>% 
  group_by(ID) %>% summarise(nNA=sum(is.na(value))) %>% 
  arrange(desc(nNA)) %>% View


# fill NAs

genesCallsFinal[is.na(genesCallsFinal)] %>% length()

genesCallsFinal[is.na(genesCallsFinal)] <- "2.009|0.009"

genesCallsFinal[is.na(genesCallsFinal)] %>% length()


genesCallsFinal %>% names

genesCallsFinal <- genesCallsFinal %>% 
 rename(BCL6 = "chr_3q_BCL6",
                                              BIRC2 = "chr_11q_BIRC2",
                                              CDKN2C = "chr_1p_CDKN2C",
                                              CYLD = "chr_16q_CYLD",
                                              EVI5 = "chr_1p_EVI5",
                                              MAX = "chr_14q_MAX",
                                              MCL1 = "chr_1q_MCL1",
                                              MYC_amp = "chr_8q_MYC_amp",
                                              MYC_del = "chr_8q_MYC_del",
                                              NSD2 = "chr_4p_NSD2",
                                              SP140 = "chr_2q_SP140",
                                              TENT5C = "chr_1p_TENT5C", 
                                              TERC = "chr_3q_TERC",
                                              TNFRSF17 = "chr_16p_BCMA",
                                              XBP1 = "chr_22q_XBP1"
)




#============== call arm level calls ================

segm_arm <- segm
segm_arm$ID <- segm_arm$sample

# import Purity file, needed to compute normal threshold for segments

# import Purity file, needed to compute threshold for normal/altered segments
pur_tab_TN <- fread("data/Bustoros_WES/TN_RAPH_purity_table.txt")
pur_tab_TO <- fread("data/Bustoros_WES/TO_RAPH_purity_table.txt")

pur_tab <- rbind(pur_tab_TN, pur_tab_TO)
pur_tab$purity[is.na(pur_tab$purity)] <- 1

pur_tab <- pur_tab %>% dplyr::select(ID= sample, purity)


# CALL arm level events
#_______ loop per sample ________
res_armcalls <- data.frame()

i=samples[1]
for(i in samples){
  
  message(which(i==samples), "/",length(samples)," - ",i)
  
  segments <- segm_arm %>%  filter(ID==i)
  
  Purity <- pur_tab %>% filter(ID == i) %>% .$purity
  
  thresh <- 0.05/Purity # ADJUST LOD threshold based on Purity
  
  segments$CN_abs_dev <- abs(segments$CN - 2)
  
  ARMS <- segments$chrarm %>% unique()
  
  #____ loop per arm ____
  
  a=ARMS[2]
  CN_CHR <- data.frame()
  
  for(a in ARMS){
    
    seg_arm <- segments %>% filter(chrarm==a)
    alt_seg <- seg_arm %>% filter(CN_abs_dev > thresh & width > 5*10^6 )
    
    if(nrow(alt_seg)>0) {
      wm_CN <- weighted.mean(alt_seg$CN, w = alt_seg$width)
      wm_CI <- weighted.mean(alt_seg$seg_sigma, w = alt_seg$width)
      
    } else {
      wm_CN <- weighted.mean(seg_arm$CN, w = seg_arm$width)
      wm_CI <- weighted.mean(seg_arm$seg_sigma, w = seg_arm$width)
      
    }
    
    CN_CHR <- rbind(CN_CHR, data.frame(chrarm=a, weighted_mean_CN = wm_CN, CI95 = wm_CI ))
  }
  
  CN_CHR$ID <- i
  
  res_armcalls <- rbind(res_armcalls, CN_CHR)
}



res_armcalls$chrarm <- paste0("chr_", res_armcalls$chrarm)
res_armcalls$value <- paste0(res_armcalls$weighted_mean_CN %>% round(4),"|", res_armcalls$CI95 %>% round(4))

res_armcalls_table <- res_armcalls %>% dplyr::select(chrarm, value, ID) %>% 
  reshape2::dcast(ID~chrarm, value.var = "value")






# remove satellite arms (chr_13p, chr_22p, chr_21p)
apply(res_armcalls_table, 2, function(x) sum(is.na(x)))

res_armcalls_table2 <- res_armcalls_table %>% dplyr::select(-chr_21p )

apply(res_armcalls_table2, 2, function(x) sum(is.na(x)))

# remove few NAs
res_armcalls_table2[is.na(res_armcalls_table2)] <- "2.009|0.009"
res_armcalls_table2[is.na(res_armcalls_table2)]


ALL_CALLS <- left_join(res_armcalls_table2, genesCallsFinal)







ALL_CALLS <- left_join(res_armcalls_table, genesCallsFinal)

# save
write_tsv(ALL_CALLS, "results/JB_Paper/calls/allCalls_CN_CI_broad_focal_BUSnew.txt" )


#================ compute AMP and DEL clonality levels =======================

rownames(ALL_CALLS) <- ALL_CALLS$sample

ACN <- apply(ALL_CALLS, c(1,2), function(x) {str_remove(x, "\\|.*")} %>% as.numeric) %>% as.data.frame()

# str(ACN)
# ACN[!complete.cases(ACN[,-1]), ] %>% View
# is.na(ACN[,-1]) %>% table


# AMP clonality level calls
ALL_CALLS_AMP <- apply(ACN[,-1],
                       MARGIN = c(1,2),
                       function(x){
                         y <- x-2
                         if(y<0){0} else if(y>1){1} else {y %>% round(3)}
                       }) %>% as.data.frame()

names(ALL_CALLS_AMP) <- paste0("amp_", names(ALL_CALLS_AMP))


# DEL clonality level calls
ALL_CALLS_DEL <- apply(ACN[,-1],
                       MARGIN = c(1,2),
                       function(x){
                         y <- 2-x
                         if(y<0){0} else if(y>1){1} else {y %>% round(3)}
                       }) %>% as.data.frame()

names(ALL_CALLS_DEL) <- paste0("del_", names(ALL_CALLS_DEL))

ALL_CLON_CALLS <- cbind(ALL_CALLS_AMP, ALL_CALLS_DEL)

rownames(ALL_CLON_CALLS) <- ALL_CALLS$ID


#_______ take off redundant / collinear calls _______

genesCallsFinal %>% names
all_focal_df 

ALL_CLON_CALLS %>% names

ALL_CLON_CALLS_def <- ALL_CLON_CALLS %>%
  dplyr::select(-c(  
                   # remove broad AMP arms with focal event
                   amp_chr_1q,  # MCL1               
                   amp_chr_8q,  # MYC_amp               
                   amp_chr_3q,  # TERC               
                   amp_chr_16p, # TNFRSF17              
                   
                   # remove broad DEL arms with focal event
                   del_chr_3q,   # BCL6
                   del_chr_11q, # BIRC2              
                   del_chr_1p,  # CDKN2C, TENT5C, EVI5         
                   del_chr_16q, # CYLD
                   del_chr_14q, # MAX
                   del_chr_8q,  # MYC_del
                   del_chr_4p,  # NSD2
                   del_chr_2q,  # SP140
                   del_chr_22q, # XBP1
                   
                   # remove not appropriate focal AMPS
                   amp_chr_3q_BCL6,
                   amp_chr_11q_BIRC2, 
                   amp_chr_1p_CDKN2C, 
                   amp_chr_16q_CYLD, 
                   amp_chr_1p_EVI5, 
                   amp_chr_14q_MAX, 
                   amp_chr_8q_MYC_del, 
                   amp_chr_4p_NSD2, 
                   amp_chr_2q_SP140, 
                   amp_chr_1p_TENT5C, 
                   amp_chr_22q_XBP1,
                   
                   # remove not appropriate focal DELS
                   del_chr_1q_MCL1, 
                   del_chr_8q_MYC_amp, 
                   del_chr_3q_TERC, 
                   del_chr_16p_BCMA)) 

ALL_CLON_CALLS %>% names

ALL_CLON_CALLS_def <- ALL_CLON_CALLS_def %>% dplyr::select(str_sort(names(ALL_CLON_CALLS_def), numeric = T)) %>% rownames_to_column(var = "sample")
ALL_CLON_CALLS_def %>% names

write_tsv(ALL_CLON_CALLS_def, "results/JB_Paper/calls/clonalityCalls_broad_focal_BUSnew.txt")
