################# generate focal genes loci from Biomart ####################

library(data.table)
library(tidyverse)
library(biomaRt)

#======== define focal genes (hg38) and retrieve coordinates from BioMart =============

focal_genes_del_38 <- c( "CDKN2C", 
                         "TENT5C", # FAM46C
                         "EVI5",
                         "SP140",
                         "BCL6",
                         "NSD2",
                         "MYC",
                         "BIRC2",
                         "MAX",
                         "CYLD",
                         "XBP1"
                         )

focal_genes_amp_38 <-c("MYC",
                       "MCL1",
                       "TERC", 
                       "TNFRSF17" # BCMA
                       ) 


mart38 = useMart("ensembl", 
                 dataset="hsapiens_gene_ensembl")

# listFilters(mart38)
# listAttributes(mart38) 

hg38_focal_amps <- getBM(attributes =  c('hgnc_symbol',
                                         "chromosome_name", 
                                         "start_position", 
                                         "end_position") , 
                         filters = c('hgnc_symbol'),
                         values = focal_genes_amp_38, 
                         mart = mart38) %>% mutate(type="focal_amp")

hg38_focal_dels <- getBM(attributes =  c('hgnc_symbol',
                                         "chromosome_name", 
                                         "start_position", 
                                         "end_position") , 
                         filters = c('hgnc_symbol'),
                         values = focal_genes_del_38, 
                         mart = mart38) %>% mutate(type="focal_del")


# merge
all_focal_hg38_df <- rbind(hg38_focal_amps, hg38_focal_dels)
all_focal_hg38_df$refgenome <- "hg38"


# differentiate MYC amp e del names
all_focal_hg38_df$hgnc_symbol[all_focal_hg38_df$hgnc_symbol == "MYC" & all_focal_hg38_df$type == "focal_amp"] <- "MYC_amp"

all_focal_hg38_df$hgnc_symbol[all_focal_hg38_df$hgnc_symbol == "MYC" & all_focal_hg38_df$type == "focal_del"] <- "MYC_del"

# save
write_tsv(all_focal_hg38_df, "data/focal_loci_hg38_JBpaper.txt")
