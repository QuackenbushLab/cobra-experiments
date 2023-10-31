setwd("/home/soel/cobra/")
library(AnnotationDbi)
library(cluster)
library(clusterProfiler)
library(corpcor)
library(DESeq2)
library(fastDummies)
library(flashClust)
library(ggplot2)
library(gplots)
library(GOstats)
library(igraph)
library(KEGGandMetacoreDzPathwaysGEO)
library(MASS)
library(NetworkToolbox)
library(org.Hs.eg.db)
library(psych)
library(rARPACK)
library(recount3)
library(r2r)
library(WGCNA)

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
cobra_res <- cobra(X, vsd, method = "pcorsh")
Sigma_D <- abs(cobra_res$Q %*% diag(cobra_res$psi[2,]) %*% t(cobra_res$Q))
rownames(Sigma_D) <- rownames(vsd)
colnames(Sigma_D) <- rownames(vsd)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(Sigma_D, powerVector = powers, verbose = 5)
png("pre-processing.png")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 3;
adjacency = adjacency(Sigma_D, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = 200
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
write.csv(dynamicColors, "./data/thca/dragon/module_assignment.txt")

pdf("./figures/thca/dragon/modules_size.pdf")
#dynamicColors <- read.csv("./data/thca/dragon/module_assignment.txt", row.names = 1)
M <- hashmap()
M[['grey']] <- 'unassigned'
cnt <- 0
for(i in seq_len(length(table(dynamicColors)))){
  color <- rownames(table(dynamicColors))[i]
  #print(color)
  if(color != 'grey'){
    M[[color]] <- as.integer(cnt)
    cnt <- cnt + 1
  }
}
for(i in seq_len(dim(dynamicColors)[1])){
  dynamicColors[i, 2] <- M[[dynamicColors[i, 1]]]
}
barplot(table(dynamicColors[,2]), names.arg = c(as.character(seq_len(7)), "unassigned"), ylab = 'Number of genes', main = 'Module sizes with DRAGON after COBRA correction', las = 2, cex.names = .8)
dev.off()

png("cluster_wgcna.png")
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

enrichmentGO <- function(module, path){
  ans.go <- enrichGO(gene = unlist(mget(as.character(lapply(module, substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)), ont = "BP",
                     OrgDb ="org.Hs.eg.db",
                     universe = unlist(mget(as.character(lapply(rownames(vsd), substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                     readable=TRUE,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH"
  )
  if(dim(ans.go)[1] == 0){
    return(c())
  }
  dotplot(ans.go, showCategory=20)
  ggsave(file = path)
  if(dim(ans.go)[1] == 1){
    return(c(ans.go$Description[1]))
  }
  c(ans.go$Description[1], ans.go$Description[2])
}

enrichmentKEGG <- function(module, path){
  ans.kegg <- enrichKEGG(gene = unlist(mget(as.character(lapply(module, substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                         universe = unlist(mget(as.character(lapply(rownames(vsd), substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH"
  )
  if(dim(ans.kegg)[1] == 0){
    return(c())
  }
  dotplot(ans.kegg, showCategory=20)
  ggsave(file = path)
  if(dim(ans.kegg)[1] == 1){
    return(c(ans.kegg$Description[1]))
  }
  c(ans.kegg$Description[1], ans.kegg$Description[2])
}

pathways_go <- c()
for(color in rownames(table(dynamicColors[,1]))){
  if(color != 'grey'){
    pathways_go <- c(pathways_go, 
                         enrichmentGO(rownames(vsd)[dynamicColors[,1] == color], paste("./figures/thca/dragon/GO_", color, ".pdf", sep = "")))
  }
}

pathways_kegg <- c()
for(color in rownames(table(dynamicColors[,1]))){
  if(color != 'grey'){
    pathways_kegg <- c(pathways_kegg, 
                         enrichmentKEGG(rownames(vsd)[dynamicColors[,1] == color], paste("./figures/thca/dragon/KEGG_", color, ".pdf", sep = "")))
  }
}

data <- matrix(0, nrow = length(rownames(table(dynamicColors[,1]))) - 1, ncol = length(unique(pathways_go)))
idx <- 1
for(i in seq_len(length(rownames(table(dynamicColors[,1]))))){
  print(i)
  color <- rownames(table(dynamicColors[,1]))[i]
  if(color != 'grey'){
    ans.go <- enrichGO(gene = unlist(mget(as.character(lapply(rownames(vsd)[dynamicColors[,1] == color], substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)), ont = "BP",
                       OrgDb ="org.Hs.eg.db",
                       universe = unlist(mget(as.character(lapply(rownames(vsd), substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                       readable=TRUE,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
    v <- cbind(c(ans.go$Description), c(ans.go$p.adjust))
    for(j in seq_len(length(unique(pathways_go)))){
      pathway <- unique(pathways_go)[j]
      if(length(v[v[,1] == pathway,])[1] > 0){
        p_val <- v[v[,1] == pathway,][2]
        data[idx, j] <- -log10(as.numeric(p_val))
      }
    }
    idx <- idx + 1
  }
}
cols <- unique(pathways_go)
colnames(data) <- unique(cols)
rownames(data) <- paste0("Module ", 1:(length(rownames(table(dynamicColors[,1]))) - 1))
my_colors <- colorRampPalette(c("lightgrey", "blue"))
data[data == 0] <- NA
pdf("./figures/thca/dragon/GO_summary.pdf", width = 45, height = 45)
gplots::heatmap.2(t(data), scale = "none", Colv=NA, Rowv=NA, col = my_colors(100), revC = TRUE, cexRow = 3, cexCol = 2, margins =c(75,70),
                  trace = "none", density.info = "none", dendrogram='none')
dev.off()

data <- matrix(0, nrow = length(rownames(table(dynamicColors[,1]))) - 1, ncol = length(unique(pathways_kegg)))
idx <- 1
for(i in seq_len(length(rownames(table(dynamicColors[,1]))))){
  print(i)
  color <- rownames(table(dynamicColors[,1]))[i]
  if(color != 'grey'){
    ans.kegg <- enrichKEGG(gene = unlist(mget(as.character(lapply(rownames(vsd)[dynamicColors == color], substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                           universe = unlist(mget(as.character(lapply(rownames(vsd), substr, 1, 15)), envir=org.Hs.egENSEMBL2EG, ifnotfound = NA)),
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH"
    )
    v <- cbind(c(ans.kegg$Description), c(ans.kegg$p.adjust))
    for(j in seq_len(length(unique(pathways_kegg)))){
      pathway <- unique(pathways_kegg)[j]
      if(length(v[v[,1] == pathway,])[1] > 0){
        p_val <- v[v[,1] == pathway,][2]
        data[idx, j] <- -log10(as.numeric(p_val))
      }
    }
    idx <- idx + 1
  }
}
colnames(data) <- unique(pathways_kegg)
rownames(data) <- paste0("Module ", 1:(length(rownames(table(dynamicColors[,1])))-1))
my_colors <- colorRampPalette(c("lightgrey", "blue"))
data[data == 0] <- NA
pdf("./figures/thca/dragon/KEGG_summary.pdf", width = 45, height = 45)
gplots::heatmap.2(t(data), scale = "none", Colv=NA, Rowv=NA, col = my_colors(100), revC = TRUE, cexRow = 3, cexCol = 3, margins =c(70,70),
                  trace = "none", density.info = "none", dendrogram='none')
dev.off()