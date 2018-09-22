# MRexperiment creation for pilot data

setwd("/Users/domebraccia/Desktop/UMD/hcbravolab/hmp2-metaviz/IBD")

library(metagenomeSeq)

hmp2_metadata <- read.csv("data/hmp2_metadata.csv")
pilot_biom <- load_biom("data/pilot_taxonomic_profiles.biom")

# removal of 1 sample (probably control sample)
assayData(pilot_biom)$counts <- assayData(pilot_biom)$counts[, -which(colnames(assayData(pilot_biom)$counts) == "ESM5MEET")]

# initializing some objects
features <- rownames(pilot_biom)
df <- data.frame()
df <- rbind(df, strsplit(features[1], split = ";")[[1]])
uu <- rep("a", nrow(MRcounts(pilot_biom)))

# making list uu contain all of the taxonomic data from the `features` object
for (i in 1:nrow(MRcounts(pilot_biom))) {
  mm <- strsplit(as.character(features[i]), split = "; ")
  uu[i] <- mm
}

# changing the fData to make them more descriptive/organized?
ii <- as.data.frame(do.call(rbind, uu))
feature_order <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
colnames(ii) <- feature_order
r_names <- ii[,"species"]
fData(pilot_biom) <- ii
rownames(fData(pilot_biom)) <- r_names
rownames(assayData(pilot_biom)$counts) <- r_names

newFData <- ii
rownames(newFData) <- r_names

##################################################################

# getting relevant metadata
stool_subset <- hmp2_metadata[which(hmp2_metadata$data_type == "stool_16S"),]
subset_metadata <- stool_subset[which(stool_subset$External.ID %in% colnames(MRcounts(pilot_biom))),]
subset_samples_with_metadata <- assayData(pilot_biom)$counts[,which(colnames(assayData(pilot_biom)$counts) %in% as.character(subset_metadata$External.ID))]

# adding phenoData to the pilot stool mrexp
sample_metadata_fields <- c("External.ID","diagnosis", "Age.at.diagnosis")
phenoData <- subset_metadata[,sample_metadata_fields]
rownames(phenoData) <- phenoData$External.ID
ordering <- sampleNames(assayData(pilot_biom))
phenoData <- phenoData[ordering,]
pData(pilot_biom) <- phenoData

# creadtion of the new MRexperiment object to be filtered, normed, and analyzed below
nn <- metagenomeSeq::newMRexperiment(counts = assayData(pilot_biom)$counts, phenoData = AnnotatedDataFrame(phenoData), featureData = AnnotatedDataFrame(newFData))

# filtering and normalization
filtered_remove_duplicates_ibd_16s_biom <- filterData(nn, present = 5)
pilot_mrexp <- cumNorm(filtered_remove_duplicates_ibd_16s_biom, p = .75)