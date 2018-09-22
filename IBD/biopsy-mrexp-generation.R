# MRexperiment creation for biopsy data

setwd("/Users/domebraccia/Desktop/UMD/hcbravolab/hmp2-metaviz/IBD")

library(metagenomeSeq)

hmp2_metadata <- read.csv("data/hmp2_metadata.csv")
biopsy_biom <- load_biom("data/biopsy_taxonomic_profiles.biom")

# removing two problem samples
ibd_16s_samples <- colnames(assayData(biopsy_biom)$counts[,-which(colnames(assayData(biopsy_biom)$counts) %in% c("206672", "206718"))])
ibd_metadata_biopsy_16s <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"),]
ibd_16s_metadata <- ibd_metadata_biopsy_16s[which(ibd_metadata_biopsy_16s$External.ID %in% ibd_16s_samples),]

# removing duplicates from the metadata
remove_duplicates_metadata <- ibd_16s_metadata

# updating rownames of metadata and biom objs with dups removed
rownames(remove_duplicates_metadata) <- remove_duplicates_metadata$External.ID
remove_duplicates_biopsy_biom <- biopsy_biom[,rownames(remove_duplicates_metadata)]

# ordering sample names in metadata obj and mrexp obj
ordered_removed_dup_metadata <- remove_duplicates_metadata[colnames(assayData(remove_duplicates_biopsy_biom)$counts),]
phenoData(remove_duplicates_biopsy_biom) <- AnnotatedDataFrame(ordered_removed_dup_metadata)

mr_feature_data <- fData(remove_duplicates_biopsy_biom)

uu <- rep("a", nrow(mr_feature_data['taxonomy']))
for(k in seq(1, nrow(mr_feature_data['taxonomy']))){
  mm <- strsplit(as.character(mr_feature_data['taxonomy'][k,]), ";")
  uu[k] <- mm
}
feature_order <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

ii <- as.data.frame(do.call(rbind, uu))
ii$species <- rownames(mr_feature_data)
colnames(ii) <- feature_order

ii_df <- as.data.frame(ii)
fData(remove_duplicates_biopsy_biom) = ii[,feature_order]
rownames(featureData(remove_duplicates_biopsy_biom)) <- rownames(mr_feature_data)

# filtering and normalization
filtered_remove_duplicates_ibd_16s_biom <- filterData(remove_duplicates_biopsy_biom, present = 5)
biopsy_mrexp <- cumNorm(filtered_remove_duplicates_ibd_16s_biom, p = .75)