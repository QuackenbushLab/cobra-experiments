setwd("/home/soel/cobra/")
library(DESeq2)
library(fastDummies)
library(r2r)
library(recount3)

# Extracting data from TCGA
data <- recount3::create_rse_manual(
  project = "THCA",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
# Gene expression pre-processing, for more information see 
# https://github.com/netZoo/netbooks/blob/main/netbooks/netZooR/gene_expression_for_coexpression_nets.ipynb
G <- transform_counts(data, by = "mapped_reads")
G <- G[data@rowRanges@elementMetadata@listData$gene_type == "protein_coding",]
G <- G[-which(rowSums(G) <= 1),] # Filtering: remove genes with no counts
countMat=SummarizedExperiment::assay(DESeqDataSetFromMatrix(G, data.frame(row.names=seq_len(ncol(G))), ~1), 1)
vsd <- vst(countMat, blind=FALSE)

# Extract metadata that we want to include in COBRA
metadata_url <- locate_url(
  "THCA",
  "data_sources/tcga")
metadata <- read_metadata(file_retrieve(url = metadata_url))

sex <- metadata$tcga.gdc_cases.demographic.gender
race <- metadata$tcga.gdc_cases.demographic.race
stage <- metadata$tcga.gdc_cases.diagnoses.tumor_stage
batch <- metadata$tcga.cgc_case_batch_number
age <- metadata$tcga.cgc_case_age_at_diagnosis
cancer <- metadata$tcga.gdc_cases.samples.sample_type
cancer <- ifelse(cancer == "Solid Tissue Normal", 0, 1)

X <- cbind(rep(1, dim(G)[2]), cancer, sex == 'female', age, as.matrix(dummy_cols(race)[-1]), as.matrix(dummy_cols(stage)[-1]), as.matrix(dummy_cols(batch)[-1]))
write.csv(vsd, "thca_expression.csv")
write.csv(X, "X.csv")