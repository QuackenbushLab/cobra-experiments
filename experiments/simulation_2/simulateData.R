
simulateStudy <- function(numGenes, numSamples, addedError, blockSeq, mu, caseControl, batches, batchEffectMultiplier=1){
    start <- Sys.time()
    blocks <- sort(rep(blockSeq, length.out=numGenes))
    randomUnif <- runif(numGenes,-1,1)
    randomUnif <- randomUnif + sign(randomUnif)*.5
    blockA <- as.numeric(blocks=="A")*randomUnif # All Samples
    blockB <- as.numeric(blocks=="B")*randomUnif # Batch 1 only
    blockC <- as.numeric(blocks=="C")*randomUnif # Batch 2 only
    blockD <- as.numeric(blocks=="D")*randomUnif # Cases only
    blockE <- as.numeric(blocks=="E")*randomUnif # Controls only
    blockF <- as.numeric(blocks=="F")*randomUnif # Cases and negative with D 
    blockG <- as.numeric(blocks=="G")*randomUnif # No Samples
    blockH <- as.numeric(blocks=="H")*randomUnif # All Samples
    blockI <- as.numeric(blocks=="I")*randomUnif # No Samples
    blockJ <- as.numeric(blocks=="J")*randomUnif # No Samples
    
    batch1Effect   <- cbind(blockA, blockB, blockH)
    batch2Effect   <- cbind(blockA, blockC, blockH)
    casesEffect    <- cbind(blockD - blockF)
    controlsEffect <- cbind(blockE)
    
    SigmaBatch1 <- batchEffectMultiplier*tcrossprod(batch1Effect)
    SigmaBatch2 <- batchEffectMultiplier*tcrossprod(batch2Effect)
    SigmaCase <- tcrossprod(casesEffect)
    SigmaControl <- tcrossprod(controlsEffect)
    
    
    Sigmas <- list(Batch1Control=SigmaBatch1+SigmaControl, 
                   Batch1Case=SigmaBatch1+SigmaCase, 
                   Batch2Case=SigmaBatch2+SigmaCase, 
                   Batch2Control=SigmaBatch2+SigmaControl)
    Sigmas <- lapply(Sigmas, function(x){
        x<-x/addedError
        diag(x) <- 1 # adding a bit of random error
        x
    })
    counts <- list(Batch1Control=sum(batches==0&caseControl==0),
                   Batch1Case=sum(batches==0&caseControl==1),
                   Batch2Case=sum(batches==1&caseControl==1),
                   Batch2Control=sum(batches==1&caseControl==0))
    data <- cbind(t(mvrnorm(counts[['Batch1Control']],mu=mu, Sigma = Sigmas[['Batch1Control']])),
                  t(mvrnorm(counts[['Batch1Case']],mu=mu, Sigma = Sigmas[['Batch1Case']])),
                  t(mvrnorm(counts[['Batch2Case']],mu=mu, Sigma = Sigmas[['Batch2Case']])),
                  t(mvrnorm(counts[['Batch2Control']],mu=mu, Sigma = Sigmas[['Batch2Control']])))
    data[data<0] <-0
    
    print(paste("Data generation in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    batchEffectedGenes <- (rowSums(batch1Effect)-rowSums(batch2Effect))!=0
    realEffectedGenes <- (rowSums(casesEffect)-rowSums(controlsEffect))!=0
    list(data=data, blocks=blocks, realEffectedGenes=realEffectedGenes, batchEffectedGenes=batchEffectedGenes, 
         trueEffects=list(
        batch1Effect=batch1Effect,
        batch2Effect=batch2Effect,
        casesEffect=casesEffect,
        controlsEffect=controlsEffect))
}