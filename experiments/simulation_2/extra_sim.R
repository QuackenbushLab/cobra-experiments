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

#setwd("/Users/soel/Documents/cobra-experiments/") # Put your local path here
## Basic setup

seed <- sample(10000,1)
set.seed(5620)
numGenes <- 4000 
numSamples <-400
addedError <- 8
batchEffectMultiplier <- 2
path <- "figures/simulation_2/extra_sim"

mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
X <- cbind(rep(1,numSamples),batches)
blockSeq <- sample(LETTERS[1:10],50,replace=T)

blocks <- sort(rep(blockSeq, length.out=numGenes))
randomUnif <- runif(numGenes,-1,1)
randomUnif <- randomUnif + sign(randomUnif)*.5
blockA <- as.numeric(blocks=="A")*randomUnif # Batch 1 only
blockB <- as.numeric(blocks=="B")*randomUnif # Batch 1 only
blockC <- as.numeric(blocks=="C")*randomUnif # Batch 2 only
blockD <- as.numeric(blocks=="D")*randomUnif # Batch 1 only
blockE <- as.numeric(blocks=="E")*randomUnif # Batch 2 only
blockF <- as.numeric(blocks=="F")*randomUnif # All Samples
blockG <- as.numeric(blocks=="G")*randomUnif # Batch2
blockH <- as.numeric(blocks=="H")*randomUnif # Batch2 only
blockI <- as.numeric(blocks=="I")*randomUnif # Batch1 only
blockJ <- as.numeric(blocks=="J")*randomUnif # All Samples

batch1Effect   <- cbind(blockA, blockB, blockD, blockI, blockJ, blockF)
batch2Effect   <- cbind(blockC, blockE, blockG, blockH, blockJ, blockF)

SigmaBatch1 <- batchEffectMultiplier*tcrossprod(batch1Effect)
SigmaBatch2 <- batchEffectMultiplier*tcrossprod(batch2Effect)

Sigmas <- list(Batch1=SigmaBatch1, 
               Batch2=SigmaBatch2)
Sigmas <- lapply(Sigmas, function(x){
  x<-x/addedError
  diag(x) <- 1 # adding a bit of random error
  x
})
counts <- list(Batch1=sum(batches==0),
               Batch2=sum(batches==1))
data <- cbind(t(mvrnorm(counts[['Batch1']],mu=mu, Sigma = Sigmas[['Batch1']])),
              t(mvrnorm(counts[['Batch2']],mu=mu, Sigma = Sigmas[['Batch2']])))
data[data<0] <-0

batchEffectedGenes <- rowSums(batch2Effect)!=0
realEffectedGenes <- (rowSums(batch1Effect))!=0
study <- list(data=data, blocks=blocks, realEffectedGenes=realEffectedGenes, batchEffectedGenes=batchEffectedGenes, 
     trueEffects=list(
       batch1Effect=batch1Effect,
       batch2Effect=batch2Effect))

# Recreate the truth
batchMat <- tcrossprod(study$trueEffects$batch2Effect)
realMat <- tcrossprod(study$trueEffects$batch1Effect)

truePairwiseLabels <- rep("Background",choose(numGenes,2))
truePairwiseLabels[batchMat[row(batchMat) > col(batchMat)]!=0] <- "Batch effect"
truePairwiseLabels[realMat[row(realMat) > col(realMat)]!=0] <- "Real effect"

trueGeneLabels <- rep("Background", numGenes)
trueGeneLabels[study$batchEffectedGenes] <- "Batch"
trueGeneLabels[study$realEffectedGenes] <- "Real"

insilico_result <- cobra(X, data)
correlationNaive <- cor(t(insilico_result$G))
correlationNaivewBatch <- 
  .5 * (cor(t(insilico_result$G[,batches==1]))+cor(t(insilico_result$G[,batches==0])))
expr_limma <- removeBatchEffect(insilico_result$G, batches==1)
limma <- cor(t(expr_limma))
RUV <- cor((RUVNaiveRidge(t(insilico_result$G), nu = 5, kW = 50)))
expr_combat = ComBat(dat=insilico_result$G, batch=(batches==1), par.prior=TRUE, prior.plots=FALSE)
combat <- cor(t(expr_combat))
nsv=num.sv(study$data,cbind(rep(1,numSamples),batches), method = "be")
pc_corrected = t(sva_network(t(study$data), nsv))
sva <- cor(t(pc_corrected))

cobra_corrected <- insilico_result$Q%*%diag(insilico_result$psi[1,])%*%t(insilico_result$Q)

insilico_MasterDF <- data.frame(newMeth=cobra_corrected[row(cobra_corrected) > col(cobra_corrected)],
                              naiveMeth=correlationNaive[row(correlationNaive) > col(correlationNaive)],
                              naiveWBatch=correlationNaivewBatch[row(correlationNaivewBatch) > col(correlationNaivewBatch)],
                              combat=combat[row(combat) > col(combat)],
                              limma=limma[row(limma) > col(limma)],
                              sva=sva[row(sva) > col(sva)],
                              ruv=RUV[row(RUV) > col(RUV)],
                              labels=truePairwiseLabels)

mean(abs(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Real effect"]))
mean(abs(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Background"]))
mean(abs(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Batch effect"]))

mean(abs(insilico_MasterDF$naiveMeth[insilico_MasterDF$labels=="Real effect"]))
mean(abs(insilico_MasterDF$naiveMeth[insilico_MasterDF$labels=="Background"]))
mean(abs(insilico_MasterDF$naiveMeth[insilico_MasterDF$labels=="Batch effect"]))

differentialCorrelationsDF <- insilico_MasterDF

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


onlyEffects <- differentialCorrelationsDF[differentialCorrelationsDF$labels!="Background",]

########

plotROC <- function(corrDF, positive, plottitle="Title"){
  library(ROCR)
  corrDF$newMeth <- abs(corrDF$newMeth)
  corrDF$naiveMeth <- abs(corrDF$naiveMeth)
  corrDF$naiveWBatch <- abs(corrDF$naiveWBatch)
  corrDF$limma <- abs(corrDF$limma)
  corrDF$combat <- abs(corrDF$combat)
  corrDF$sva <- abs(corrDF$sva)
  corrDF$ruv <- abs(corrDF$ruv)
  
  methodPred  <- prediction(corrDF$newMeth, corrDF$labels%in%positive)
  roc.methodPred  <- performance(methodPred, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
  
  methodPredNaive  <- prediction(corrDF$naiveMeth, corrDF$labels%in%positive)
  roc.methodPred.naive  <- performance(methodPredNaive, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.naive  <- performance(methodPredNaive, "auc")@y.values[[1]]
  
  methodPredNaiveWBatch  <- prediction(corrDF$naiveWBatch, corrDF$labels%in%positive)
  roc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, "auc")@y.values[[1]]
  
  methodlimma  <- prediction(corrDF$limma, corrDF$labels%in%positive)
  roc.methodPred.limma  <- performance(methodlimma, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.limma  <- performance(methodlimma, "auc")@y.values[[1]]
  
  methodcombat  <- prediction(corrDF$combat, corrDF$labels%in%positive)
  roc.methodPred.combat  <- performance(methodcombat, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.combat  <- performance(methodcombat, "auc")@y.values[[1]]
  
  methodsva  <- prediction(corrDF$sva, corrDF$labels%in%positive)
  roc.methodPred.sva  <- performance(methodsva, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.sva  <- performance(methodsva, "auc")@y.values[[1]]
  
  corrDF$ruv[is.na(corrDF$ruv)] <- 0
  methodruv  <- prediction(corrDF$ruv, corrDF$labels%in%positive)
  roc.methodPred.ruv  <- performance(methodruv, measure = c("tpr"), x.measure = "fpr")
  auc.methodPred.ruv  <- performance(methodruv, "auc")@y.values[[1]]
  
  plot(roc.methodPred, main=plottitle, col = "red", lwd=3)
  lines(roc.methodPred.naive@x.values[[1]], roc.methodPred.naive@y.values[[1]], col = "green", lwd=3)
  lines(roc.methodPred.naive.w.batch@x.values[[1]], roc.methodPred.naive.w.batch@y.values[[1]], col = "blue", lwd=3)
  lines(roc.methodPred.limma@x.values[[1]], roc.methodPred.limma@y.values[[1]], col = "pink", lwd=3)
  lines(roc.methodPred.combat@x.values[[1]], roc.methodPred.combat@y.values[[1]], col = "purple", lwd=3)
  lines(roc.methodPred.sva@x.values[[1]], roc.methodPred.sva@y.values[[1]], col = "darkgreen", lwd=3)
  lines(roc.methodPred.ruv@x.values[[1]], roc.methodPred.ruv@y.values[[1]], col = "gray", lwd=3)
  
  legend("bottomright", c(paste("COBRA",round(auc.methodPred,4)), 
                          paste("Naive",round(auc.methodPred.naive,4)), 
                          paste("Naive Batch",round(auc.methodPred.naive.w.batch,4)),
                          paste("Limma",round(auc.methodPred.limma,4)),
                          paste("ComBat",round(auc.methodPred.combat,4)),
                          paste("SVA",round(auc.methodPred.sva,4))),
                          # paste("RUVCorr",round(auc.methodPred.ruv,4))),
         lty=1,lwd=1,col=c("red", "green", "blue", "pink", "purple", "darkgreen"),title="Area under ROC curve")
  abline(0,1)
}
unique(differentialCorrelationsDF$labels)
png(paste0(path,'/OursVsOthersROC.png'), width = 800, height = 400)
par(mfrow=c(1,1))
plotROC(onlyEffects, positive=c("Real effect","Negative effect"), "Real Effects vs Batch Effects")
dev.off()

# OTHER PLOT

cobra_corrected <- insilico_result$Q%*%diag(insilico_result$psi[2,])%*%t(insilico_result$Q)

insilico_MasterDF <- data.frame(newMeth=cobra_corrected[row(cobra_corrected) > col(cobra_corrected)],
                                naiveMeth=correlationNaive[row(correlationNaive) > col(correlationNaive)],
                                naiveWBatch=correlationNaivewBatch[row(correlationNaivewBatch) > col(correlationNaivewBatch)],
                                combat=combat[row(combat) > col(combat)],
                                limma=limma[row(limma) > col(limma)],
                                sva=sva[row(sva) > col(sva)],
                                ruv=RUV[row(RUV) > col(RUV)],
                                labels=truePairwiseLabels)
differentialCorrelationsDF <- insilico_MasterDF

methodPred  <- prediction(abs(differentialCorrelationsDF$newMeth), differentialCorrelationsDF$labels%in%c("Batch effect"))
roc.methodPred  <- performance(methodPred, measure = c("tpr"), x.measure = "fpr")
auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
print(auc.methodPred)