library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)
library(rARPACK)
library(limma)
library(netZooR)
library(RUVcorr)
library(sva)

## Imports

setwd("/home/soel/cobra/") # Put your local path here
source('./experiments/simulation_2/diagnosis_plots.R')
source('./experiments/simulation_2/simulateData.R')
source('./experiments/simulation_2/generateMasterDF.R')

## Basic setup


seed <- sample(10000,1)
set.seed(5620)
numGenes <- 4000 
numSamples <-400
addedError <- 8
batchEffectMultiplier <- 2
path <- "figures/simulation_2/"

mu <- rnorm(numGenes,mean = 9)

batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
caseControl <- c(rep(0,numSamples/4),rep(1,numSamples/2),rep(0,numSamples/4))
X <- cbind(rep(1,numSamples),batches, caseControl)
blockSeq <- sample(LETTERS[1:10],50,replace=T)

study <- simulateStudy(numGenes=numGenes, numSamples=numSamples, addedError=addedError, 
                       blockSeq=blockSeq, mu=mu, caseControl=caseControl, batches=batches, batchEffectMultiplier)

# Recreate the truth
batchMat <- tcrossprod(study$trueEffects$batch1Effect) - tcrossprod(study$trueEffects$batch2Effect) 
realMat <- tcrossprod(study$trueEffects$casesEffect) - tcrossprod(study$trueEffects$controlsEffect)

truePairwiseLabels <- rep("Background",choose(numGenes,2))
truePairwiseLabels[batchMat[row(batchMat) > col(batchMat)]!=0] <- "Batch effect"
truePairwiseLabels[realMat[row(realMat) > col(realMat)]!=0] <- "Real effect"

trueGeneLabels <- rep("Background", numGenes)
trueGeneLabels[study$batchEffectedGenes] <- "Batch"
trueGeneLabels[study$realEffectedGenes] <- "Real"

coex <- cor(t(study$data))
diag(coex) <- NA
png(paste0(path,'/coex_heatmap.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", 
          labRow=trueGeneLabels[c(T,F,F,F)], col="bluered", dendrogram = "none", RowSideColors = cbPalette[as.factor(trueGeneLabels[c(T,F,F,F)])])
dev.off()
insilico_result <- cobra(X, study$data)
differentialCorrelationNaive <- cor(t(insilico_result$G[,caseControl==1]))-cor(t(insilico_result$G[,caseControl==0]))
differentialCorrelationNaivewBatch <- 
  (cor(t(insilico_result$G[,caseControl==1&batches==1]))-cor(t(insilico_result$G[,caseControl==0&batches==1])) +
     cor(t(insilico_result$G[,caseControl==1&batches==0]))-cor(t(insilico_result$G[,caseControl==0&batches==0])))/2
expr_limma <- removeBatchEffect(insilico_result$G, batches==1)
limma <- cor(t(expr_limma[,caseControl == 1])) - cor(t(expr_limma[,caseControl == 0]))
RUV <- t(RUVNaiveRidge(t(insilico_result$G), center=FALSE, seq_len(numGenes)[trueGeneLabels == "Background"], nu = 5, kW = 50))
expr_combat = ComBat(dat=insilico_result$G, batch=(batches==1), par.prior=TRUE, prior.plots=FALSE)
combat <- cor(t(expr_combat[,caseControl == 1])) - cor(t(expr_combat[,caseControl == 0]))
nsv=num.sv(study$data,cbind(rep(1,numSamples),batches), method = "be")
pc_corrected = t(sva_network(t(study$data), 250))
sva <- cor(t(pc_corrected[,batches == 1])) - cor(t(pc_corrected[,batches == 0]))
insilico_MasterDF <- generateMasterDF(insilico_result, differentialCorrelationNaive, differentialCorrelationNaivewBatch, combat, limma, sva, RUV, truePairwiseLabels, maxPoints=800000)

plotEigenvectors(insilico_result, 
                 trueGeneLabels,
                 path, numEigenvectors=6)

mean(abs(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Real effect"]))
mean(abs(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Batch effect"]))

diagnosticPlots(insilico_MasterDF, path)