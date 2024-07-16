setwd("/home/soel/cobra/")
library(AnnotationDbi)
library(cluster)
library(clusterProfiler)
library(DESeq2)
library(fastDummies)
library(flashClust)
library(ggplot2)
library(GOstats)
library(igraph)
library(KEGGandMetacoreDzPathwaysGEO)
library(netZooR)
library(org.Hs.eg.db)
library(psych)
library(r2r)
library(NetworkToolbox)
library(recount3)
library(WGCNA)
library(glue)
library(MASS)

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
start.time <- Sys.time()
c <- cobra(X, vsd)
Sigma_D <- abs(c$Q %*% diag(c$psi[2,]) %*% t(c$Q))
print(paste("COBRA in ",round(as.numeric(difftime(Sys.time(), start.time,units = "secs")),1), "seconds"))

start.time <- Sys.time()
C <- cor(t(vsd))
end.time <- Sys.time()
print(paste("Co-expression in ",round(as.numeric(difftime(Sys.time(), start.time,units = "secs")),1), "seconds"))


size <- c(1000, 10000, 20000, 22000)
G <- transform_counts(data, by = "mapped_reads")
G <- G[data@rowRanges@elementMetadata@listData$gene_type == "protein_coding",]
G <- G[-which(rowSums(G) <= 1),] # Filtering: remove genes with no counts
countMat=SummarizedExperiment::assay(DESeqDataSetFromMatrix(G, data.frame(row.names=seq_len(ncol(G))), ~1), 1)
vsd <- vst(countMat, blind=FALSE)

for(i in c(1, 2, 3, 4)){
  print(glue('SIZE: {size[i]} genes'))
  genes <- sample(seq(1, dim(vsd)[1]), size[i], replace=T)
  E <- G[genes,]
  
  start.time <- Sys.time()
  c <- cobra(X, E)
  Sigma_D <- abs(c$Q %*% diag(c$psi[2,]) %*% t(c$Q))
  print(paste("COBRA in ",round(as.numeric(difftime(Sys.time(), start.time,units = "secs")),1), "seconds"))
  
  start.time <- Sys.time()
  C <- cor(t(E))
  end.time <- Sys.time()
  print(paste("Co-expression in ",round(as.numeric(difftime(Sys.time(), start.time,units = "secs")),1), "seconds"))
}