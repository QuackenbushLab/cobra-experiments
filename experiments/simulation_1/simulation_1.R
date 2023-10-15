library(ggplot2)
library(limma)
library(MASS)
library(netZooR)
library(reshape)
library(reshape2)
library(sva)

setwd("/home/soel/cobra/")
path_fig <- "figures/simulation_1/"
path_data <- "data/simulation_1/"
numGenes <- 1000
numBatch1 <- 100
numBatch2 <- 100
nperm <- 10
plotSubset <- c(T,rep(F,10))

makeDiffCoexDF <- function(exprData, batches){
  coex <- cor(t(exprData[,batches=="Batch1"])) - cor(t(exprData[,batches=="Batch2"]))
  coexvalues <- sort(abs(coex[row(coex)>col(coex)]))
  
  nullMatrix <- replicate(nperm,{
    cat(".")
    batchPerm <- sample(c(rep(T,numBatch2),rep(F,numBatch2)))
    coex_perm <- abs(cor(t(data[,batchPerm])) - cor(t(data[,!batchPerm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
  })
  colnames(nullMatrix) <- paste0("null",seq_len(ncol(nullMatrix)))
  diffCoex <- melt(cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues=coexvalues))
  as.data.frame(diffCoex)[plotSubset,]
}

cobra_experiment <- function(exprData){
  X <- cbind(c(rep(1, 200)), c(rep(0, numBatch1), rep(1, numBatch2)), c(rep(0, numBatch1), rep(1, numBatch2)))
  res <- cobra(X, exprData)
  coex <- res$Q%*%diag(res$psi[3,])%*%t(res$Q)
  coexvalues <- sort(abs(coex[row(coex)>col(coex)]))
  
  nullMatrix <- replicate(10,{
    cat(".")
    batchPerm <- sample(c(rep(1,numBatch1),rep(0,numBatch2)))
    X <- cbind(c(rep(1, 200)), c(rep(0, numBatch1), rep(1, numBatch2)), batchPerm)
    res <- cobra(X, exprData)
    coex_perm <- abs(res$Q%*%diag(res$psi[3,])%*%t(res$Q))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
  })
  colnames(nullMatrix) <- paste0("null",seq_len(ncol(nullMatrix)))
  diffCoex <- melt(cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues=coexvalues))
  as.data.frame(diffCoex)[plotSubset,]
}

plot_coexpression <- function(diffCoexSubsetRaw, diffCoexSubsetComBat, diffCoexSubsetLimma, diffCoexSubsetCobra){
  diffCoexSubset <- rbind(diffCoexSubsetRaw,diffCoexSubsetComBat, diffCoexSubsetLimma, diffCoexSubsetCobra)
  pdf(paste0(path_fig, "demo_diff_coex_density.pdf"), width=12, height=6)
  densplot <- ggplot(diffCoexSubset) + 
    geom_density(aes(x=value, group=Var2),color="grey",alpha=.05,fill="lightgrey") +
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="coexvalues",], aes(x=value,color="blue"),fill="blue", alpha=.2) + 
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="expected",], aes(x=value,color="black")) + 
    coord_cartesian(xlim=c(0, 1), ylim = c(0, 33)) +
    facet_wrap(~type) +
    ggtitle("Differential Coexpression") + xlab("Absolute Differential Coexpression") + ylab("Density") +
    theme_bw(base_size = 20) + 
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,0),
          legend.position=c(.99,.65),
          legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"))+
    scale_colour_manual(name = 'Partition',
                        values =c('blue'='blue','black'='black'), labels = c('Mean Random','Across Batches')) 
  print(densplot)
  dev.off()
}

plot_expression <- function(data, expr_combat, batches){
  design <- model.matrix(~batches)
  colnames(design) <- c("Intercept","Batch2")
  fit <- lmFit(data, design)
  fit_combat <- lmFit(expr_combat, design)
  fit2 <- eBayes(fit)
  fit2_combat <- eBayes(fit_combat)
  output <- topTable(fit2, coef = 2,n=Inf)
  output_combat <- topTable(fit2_combat, coef = 2,n=Inf)
  sum(output$adj.P.Val<.001)
  sum(output_combat$adj.P.Val<.001)
  rawDF <- data.frame(pvalues=output$P.Value, 
                      logFC=output$logFC, type="Uncorrected")
  combatDF <- data.frame(pvalues=output_combat$P.Value,
                         logFC=output_combat$logFC, type="After Batch Correction") 
  df <- rbind(rawDF,combatDF)
  pdf(paste0(path_fig, "demo_diffexpress.pdf"), width=12, height=6)
  diffexpressPlot <- ggplot(df) + 
    geom_point(aes(x=logFC,y=-log10(pvalues)), shape=3) + 
    facet_wrap(~type) + 
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle("Differential Expression") + xlab("Log Fold Change") + ylab("-Log p-values")
  plot(diffexpressPlot)
  dev.off()
}

experiment <- function(data){
  batches <- c(rep("Batch1",numBatch1),rep("Batch2",numBatch2))
  diffCoexSubsetRaw    <- cbind(makeDiffCoexDF(data, batches),type="Uncorrected")
  expr_combat = ComBat(dat=data, batch=batches, par.prior=TRUE, prior.plots=FALSE)
  diffCoexSubsetComBat <- cbind(makeDiffCoexDF(expr_combat, batches),type="After Combat Batch Correction")
  diffCoexSubsetComBat[diffCoexSubsetComBat$Var2 == "null1",]$value = diffCoexSubsetRaw[diffCoexSubsetRaw$Var2 == "null1",]$value
  diffCoexSubsetComBat[diffCoexSubsetComBat$Var2 == "expected",]$value = diffCoexSubsetRaw[diffCoexSubsetRaw$Var2 == "expected",]$value
  diffCoexSubsetCobra <- cbind(cobra_experiment(data), type = "After COBRA batch correction")
  diffCoexSubsetLimma <- cbind(makeDiffCoexDF(removeBatchEffect(data, batches), batches), type="After limma correction")
  plot_expression(data, removeBatchEffect(data, batches), batches)
  plot_coexpression(diffCoexSubsetRaw, diffCoexSubsetComBat, diffCoexSubsetLimma, diffCoexSubsetCobra)
}

lungData      <- read.table(paste0(path_data, "/expression_lung.csv"),sep=',', header=T,row.names=1)
data <- as.matrix(cbind(lungData[sample(nrow(lungData),numGenes),1:numBatch1],lungData[sample(nrow(lungData),numGenes),1:numBatch2]))

experiment(data)
