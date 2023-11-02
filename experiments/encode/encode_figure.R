setwd("/home/soel/cobra/")
local <- "./data/encode/"
local_fig <- "./figures/encode/"
library(DESeq2)
library(dplyr)
library(ggplot2)
library(limma)
library(netZooR)
library(readr)
library(reshape2)
library(rlang)
library(Rsubread)
library(RUVcorr)
library(sva)

expr <- read.csv(paste0(local,"rnaseq_filtered_med0genes.txt"), sep = "\t")
expr <- expr[,!grepl(pattern = "_2", colnames(expr))]
yaleBatch <- grepl(pattern = "yale",colnames(expr))
yale <- expr[, yaleBatch]
argonne <- expr[,!yaleBatch]
expr <- cbind(yale, argonne)
yaleBatch = c(rep(T, 63), rep(F, 63))

expr_limma <- removeBatchEffect(expr, yaleBatch)
design <- model.matrix(~yaleBatch)
colnames(design) <- c("Intercept","Batch2-Batch5")
fit <- lmFit(expr, design)
fit <- eBayes(fit)
fit_combat <- lmFit(expr_limma, design)
fit_combat <- eBayes(fit_combat)
output <- topTable(fit, coef = 2,n=Inf)
output_combat <- topTable(fit_combat, coef = 2,n=Inf)
sum(output$adj.P.Val<.01)
sum(output_combat$adj.P.Val<.01)

topNgenes <- 1000
expr <- expr[order(-apply(expr,1,var))[seq_len(topNgenes)],]
expr_limma <- removeBatchEffect(expr, yaleBatch)

negative_controls <- c()
for(i in seq_len(30)){
  idx <- i * 2
  tmp <- rownames(which(abs(cor(t(expr))) == as.numeric(sort(unlist(as.list(abs(cor(t(expr))))))[idx]), arr.ind = TRUE))
  negative_controls <- c(negative_controls, tmp[1], tmp[2])
}
RUV <- t(RUVNaiveRidge(t(expr), center=FALSE, negative_controls, nu = 5, kW = 50))
expr_combat = ComBat(dat=expr, batch=yaleBatch, par.prior=TRUE, prior.plots=FALSE)

B1 <- 10
B2 <- 1000
compute <- function(E){
  res <- c()
  for(rep in seq_len(B1)){
    Sigma1 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
    Sigma2 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
    for(B in seq_len(B2)){
      group1 <- sample(c(rep(T, 32), rep(F, 31)))
      group1 <- c(cbind(group1, !group1))
      coex1 <- cor(t(E[,group1])) - cor(t(E[,!group1]))
      group2 <- sample(c(rep(T, 48), rep(F, 15)))
      group2 <- c(cbind(group2, !group2))
      coex2 <- cor(t(E[,group2])) - cor(t(E[,!group2]))
      Sigma1 <- Sigma1 + coex1
      Sigma2 <- Sigma2 + coex2
    }
    res <- c(res, ks.test(Sigma1, Sigma2)$statistic)
    print(ks.test(Sigma1, Sigma2))
  }
  return(res)
}

stat_cobra <- c()
for(rep in seq_len(B1)){
  Sigma1 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
  Sigma2 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
  for(B in seq_len(B2)){
    group1 <- sample(c(rep(T, 32), rep(F, 31)))
    group1 <- c(cbind(group1, !group1))
    X <- cbind(rep(1, dim(expr)[2]), group1, yaleBatch)
    res <- cobra(X, expr)
    coex1 <- res$Q %*% diag(res$psi[2,]) %*% t(res$Q)
    
    group2 <- sample(c(rep(T, 48), rep(F, 15)))
    group2 <- c(cbind(group2, !group2))
    X <- cbind(rep(1, dim(expr)[2]), group2, yaleBatch)
    res <- cobra(X, expr)
    coex2 <- res$Q %*% diag(res$psi[2,]) %*% t(res$Q)
    Sigma1 <- Sigma1 + coex1
    Sigma2 <- Sigma2 + coex2
  }
  stat_cobra <- c(stat_cobra, ks.test(Sigma1, Sigma2)$statistic)
  print(ks.test(Sigma1, Sigma2))
}

stat_sva <- c()
for(rep in seq_len(B1)){
  Sigma1 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
  Sigma2 <- matrix(0, nrow = topNgenes, ncol = topNgenes)
  for(B in seq_len(B2)){
    group1 <- sample(c(rep(T, 32), rep(F, 31)))
    group1 <- c(cbind(group1, !group1))
    X <- cbind(rep(1, dim(expr)[2]), group1, yaleBatch)
    nsv=num.sv(expr,X, method = "be")
    pc_corrected = t(sva_network(t(expr), nsv))
    pc_diff1 <- cor(t(pc_corrected[,group1])) - cor(t(pc_corrected[,!group1]))
    
    group2 <- sample(c(rep(T, 48), rep(F, 15)))
    group2 <- c(cbind(group2, !group2))
    X <- cbind(rep(1, dim(expr)[2]), group2, yaleBatch)
    nsv=num.sv(expr,X, method = "be")
    pc_corrected = t(sva_network(t(expr), nsv))
    pc_diff2 <- cor(t(pc_corrected[,group2])) - cor(t(pc_corrected[,!group2]))
    Sigma1 <- Sigma1 + pc_diff1
    Sigma2 <- Sigma2 + pc_diff2
  }
  stat_sva <- c(stat_sva, ks.test(Sigma1, Sigma2)$statistic)
  print(ks.test(Sigma1, Sigma2))
}

stat <- t(compute(expr))
stat_limma <- t(compute(expr_limma))
stat_combat <- t(compute(expr_combat))
stat_RUV <- t(compute(RUV))

data <- cbind(t(cbind(stat, stat_limma, t(stat_cobra), stat_combat, stat_RUV, t(stat_sva))), t(cbind(t(rep("Naive", B1)), t(rep("Limma", B1)), t(rep("COBRA", B1)), t(rep("Combat", B1)), t(rep("RUVCorr", B1)), t(rep("SVA", B1)))))
colnames(data) <- c("val", "group")
rownames(data) <- seq_len(B1 * 6)
data <- as.data.frame(data)
data$val <- as.numeric(data$val)
pdf(paste0(local_fig, "stability.pdf"))
boxplot(val ~ group,
        data=data,
        main="Stability for different batch proportions",
        xlab="Method",
        ylab="D statistic",
        col="orange",
        border="brown"
)
dev.off()

expr <- read.csv(paste0(local,"rnaseq_filtered_med0genes.txt"), sep = "\t")
expr <- expr[,!grepl(pattern = "_2", colnames(expr))]
yaleBatch <- grepl(pattern = "yale",colnames(expr))
yale <- expr[, yaleBatch]
argonne <- expr[,!yaleBatch]
expr <- cbind(yale, argonne)
yaleBatch = c(rep(T, 63), rep(F, 63))

expr_limma <- removeBatchEffect(expr, yaleBatch)
negative_controls <- c()
for(i in seq_len(30)){
  idx <- i * 2
  tmp <- rownames(which(abs(cor(t(expr))) == as.numeric(sort(unlist(as.list(abs(cor(t(expr))))))[idx]), arr.ind = TRUE))
  negative_controls <- c(negative_controls, tmp[1], tmp[2])
}
RUV <- t(RUVNaiveRidge(t(expr), center=FALSE, negative_controls, nu = 5, kW = 50))
expr_combat = ComBat(dat=expr, batch=yaleBatch, par.prior=TRUE, prior.plots=FALSE)

experiment <- function(group){
  G <- c(cbind(group, !group))
  naive_diff <- cor(t(expr[,G])) - cor(t(expr[,!G]))
  limma_diff <- cor(t(expr_limma[,G])) - cor(t(expr_limma[,!G]))
  X <- cbind(rep(1, dim(expr)[2]), G, yaleBatch)
  res <- cobra(X, expr)
  cobra_diff <- res$Q %*% diag(res$psi[2,]) %*% t(res$Q)
  ruv_diff <- cor(t(RUV[,G])) - cor(t(expr[,!G]))
  pc_corrected = t(sva_network(t(expr), nsv))
  pc_diff <- cor(t(pc_corrected[,G])) - cor(t(pc_corrected[,!G]))
  combat_diff <- cor(t(expr_combat[,G])) - cor(t(expr_combat[,!G]))
  return(c(mean(naive_diff^2), mean(limma_diff^2), mean(cobra_diff^2), mean(ruv_diff^2), mean(pc_diff^2), mean(combat_diff^2)))
}

B <- 10

naive_res <- matrix(0, nrow = B, ncol = 9)
limma_res <- matrix(0, nrow = B, ncol = 9)
cobra_res <- matrix(0, nrow = B, ncol = 9)
ruv_res <- matrix(0, nrow = B, ncol = 9)
pc_res <- matrix(0, nrow = B, ncol = 9)
combat_res <- matrix(0, nrow = B, ncol = 9)
P <- c(0, 8, 16, 24, 32, 39, 47, 55, 63)
for(p in seq_len(length(P))){
  print(p)
  for(i in seq_len(B)){
    print(i)
    group <- sample(c(rep(T, P[p]), rep(F, 63-P[p])))
    res <- experiment(group)
    naive_res[i, p] <- res[1]
    limma_res[i, p] <- res[2]
    cobra_res[i, p] <- res[3]
    ruv_res[i, p] <- res[4]
    pc_res[i,p] <- res[5]
    combat_res[i,p] <- res[6]
  }
}
naive_summ <- colMeans(naive_res)
cobra_summ <- colMeans(cobra_res)
limma_summ <- colMeans(limma_res)
ruv_summ <- colMeans(ruv_res)
pc_summ <- colMeans(pc_res)
combat_summ <- colMeans(combat_res)
p <- c(0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1)

DATA <- data.frame(p = p, naive = naive_summ, cobra = cobra_summ, limma = limma_summ, ruv = ruv_summ, pc = pc_summ, combat = combat_summ)

pdf(paste0(local_fig, "encodeFigure.pdf"))
plot(p, DATA$naive, type="b", pch=19, col="green", ylim = c(0, 0.05), xlab="Ratio of Yale patients in group 1", ylab="Mean residual co-epxression", main="Residual differential co-expression between groups in ENCODE")
lines(p, DATA$limma, pch=19, col="blue", type="b")
lines(p, DATA$cobra, pch=19, col="red", type="b")
lines(p, DATA$ruv, pch=19, col="turquoise", type="b")
lines(p, DATA$combat, pch=19, col="yellow", type="b")
lines(p, DATA$pc, pch=19, col="darkgreen", type="b")
legend(0.5, 0.05, legend=c("Naive", "Limma", "Cobra", "RUVCorr", "Combat", "sva"),
       col=c("green", "blue", "red", "turquoise", "yellow", "darkgreen"), lty=1, cex=0.8)
dev.off()