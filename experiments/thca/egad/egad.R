#
setwd("/Users/soel/Documents/cobra-experiments/experiments/thca/egad")
options(stringsAsFactors = F)
library(MASS)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(flavin)
library(ontologyIndex)
library(EGAD)
library(dplyr)
library(data.table)
library(scales)
library(dismay)
library(corpcor)
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
library(WGCNA)
library(recount3)
library(DESeq2)
library(fastDummies)

# General-purpose functions

map_genes = function(dat, keys, from = "SYMBOL", to = "ENSEMBL", db = NULL) {
  if (is.null(db) | !"OrgDb" %in% class(db))
    stop("argument `db` must be an object of class OrgDb")
  if (is.null(from) | !from %in% keytypes(db))
    stop("invalid `from` argument: ", from)
  if (is.null(to) | !to %in% keytypes(db))
    stop("invalid `to` argument: ", to)
  map = suppressMessages(AnnotationDbi::select(db, keys = keys, columns = to, 
                                               keytype = from))
  values = map[[to]][match(keys, map[[from]])]
  aggregate(dat, by = list(values), FUN = sum)
}

write_and_gzip = function(data, filepath) {
  write.table(data, filepath, sep = "\t", quote = F, row.names = F)
  system(paste("gzip --force", filepath))
}

clean_metric = function(vec) {
  as.character(fct_recode(vec,
                          "Manhattan distance" = "manhattan",
                          "Euclidean distance" = "euclidean",
                          "Canberra distance" = "canberra",
                          "Zero-inflated Kendall correlation" = "zi_kendall",
                          "Kendall correlation" = "kendall",
                          "Spearman correlation" = "spearman",
                          "Pearson correlation" = "pearson",
                          "Jaccard index" = "jaccard",
                          "Dice coefficient" = "dice",
                          "Hamming distance" = "hamming",
                          "Co-dependency index" = "binomial",
                          "Biweight midcorrelation" = "bicor",
                          "Cosine distance" = "cosine",
                          "Weighted rank correlation" = "weighted_rank",
                          "ϕs" = "phi_s",
                          "ρp" = "rho_p"
  ))
}


coding = read.delim("data/protein_coding_genes.txt.gz")
get_species = function(vector) {
  subset = coding %>% dplyr::filter(gene %in% vector)
  names(which.max(table(subset$source)))
}

get_score <- function(coexpr){
  if (any(is.na(coexpr)))
    coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
  if (any(is.infinite(coexpr)))
    coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)
  coexpr = scales::rescale(coexpr, to = c(-1, 1))
  
  ontology = get_ontology("data/go-basic.obo.gz")
  human = read_gpa("data/goa_human.gpa.gz", accession = "ENSEMBL", 
                   database = org.Hs.eg.db, ontology = ontology, propagate = T)
  human = filter_breadth(human, min = 1, max = 1e4)
  go = get(get_species(colnames(coexpr)))
  
  # make EGAD coexpression network
  n = nrow(coexpr)
  genes = rownames(coexpr)
  net = matrix(rank(coexpr, na.last = "keep", ties.method = "average"), 
               nrow = n, ncol = n)
  rownames(net) = colnames(net) = genes
  net = net / max(net, na.rm = T)
  diag(net) = 1
  
  go_subset = go %>%
    dplyr::filter(ENSEMBL %in% genes) %>%
    dplyr::select(ENSEMBL, GO.ID) %>%
    as.matrix()
  
  go_terms = unique(go_subset[, "GO.ID"])
  annotations = make_annotations(go_subset, genes, go_terms)
  start.time <- Sys.time()
  gba = run_GBA(net, annotations, min = 100, max = 1000)
  end.time <- Sys.time()
  print(end.time - start.time)
  aurocs = gba[[1]][, "auc"]
  return(aurocs)
}

get_expr <- function(project, home){
  data <- recount3::create_rse_manual(
    project = project,
    project_home = paste0("data_sources/", home),
    organism = "human",
    annotation = "gencode_v26",
    type = "gene",
    recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
  )
  
  # Gene expression pre-processing, for more information see 
  # https://github.com/netZoo/netbooks/blob/main/netbooks/netZooR/gene_expression_for_coexpression_nets.ipynb
  G <- transform_counts(data, by = "mapped_reads")
  G <- G[data@rowRanges@elementMetadata@listData$gene_type == "protein_coding",]
  G <- G[-which(rowSums(G) <= 1),] # Filtering: remove genes with no counts
  countMat=SummarizedExperiment::assay(DESeqDataSetFromMatrix(G, data.frame(row.names=seq_len(ncol(G))), ~1), 1)
  vsd <- vst(countMat, blind=FALSE)
  rownames(vsd) <- substr(rownames(vsd), 1, 15)
  metadata_url <- locate_url(
    project,
    paste0("data_sources/", home),
    recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
  metadata <- read_metadata(file_retrieve(url = metadata_url))
  list(vsd=vsd, metadata=metadata)
}

pipeline <- function(project, home, control){
  data <- get_expr(project, home)
  vsd <- data$vsd
  metadata <- data$metadata
  print("Extracted data!")
  
  coexpr <- dismay::dismay(as.matrix(t(vsd)), metric = "pearson")
  print("Computed coexpression!")
  aurocs <- get_score(coexpr)
  print("Computed coexpression score!")
  
  X <- cbind(rep(1, dim(vsd)[2]))
  for(i in control){
    type <- i$type
    cov <- i$cov
    
    if(any(is.na(metadata[[cov]])))
      metadata[[cov]][is.na(metadata[[cov]])] = median(metadata[[cov]], na.rm = T)
    
    if(type == "cont"){
      X <- cbind(X, metadata[[cov]])
    }
    else{
      X <- cbind(X, as.matrix(dummy_cols(metadata[[cov]])[-1]))
    }
  }
  print("Computed COBRA design matrix!")
  method <- cobra(X, vsd)
  C <- method$Q %*% diag(method$psi[1,]) %*% t(method$Q)
  print("Computed COBRA!")
  rownames(C) <- rownames(vsd)
  colnames(C) <- rownames(vsd)
  aurocs_cobra <- get_score(C)
  
  print("Coexpression mean score...")
  print(mean(aurocs))
  print("COBRA corrected mean score...")
  print(mean(aurocs_cobra))
}

pipeline("THCA", "tcga", list(list(type="cat", cov="tcga.gdc_cases.demographic.gender"), list(type="cont", cov="tcga.cgc_case_age_at_diagnosis"), list(type="cat", cov="tcga.gdc_cases.demographic.race"), list(type="cat", cov="tcga.gdc_cases.diagnoses.tumor_stage"), list(type="cat", cov="tcga.cgc_case_batch_number"))) #0.647, 0.650

#LUNG: 0.655, 0.659
#pipeline("LUNG", "gtex", list(list(type="cont", cov="gtex.smtsisch"), list(type="cat", cov="gtex.smgebtch"), list(type="cat", cov="gtex.smcenter")))

# MUSCLE: 0.656, 0.665
#pipeline("MUSCLE", "gtex", list(list(type="cont", cov="gtex.smtsisch"), list(type="cat", cov="gtex.smgebtch"), list(type="cat", cov="gtex.smcenter")))