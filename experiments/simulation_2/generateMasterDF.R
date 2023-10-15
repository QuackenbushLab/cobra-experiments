generateMasterDF <- function(ourMethodResult, correlationResult, correlationWBatchResult, combat, limma, sva, ruv, trueLabels, maxPoints=1000000){

    X_bar <- matrix(colMeans(X),nrow=1)
    fitValues <- c(X_bar%*%ourMethodResult$psi)
    correctedCorrelation <- ourMethodResult$Q%*%diag(fitValues)%*%t(ourMethodResult$Q)
    b <- ourMethodResult$Q%*%diag(ourMethodResult$D)%*%t(ourMethodResult$Q)
    print(sum(correctedCorrelation - cor(t(ourMethodResult$G))))
    differentialCorrelation <- ourMethodResult$Q%*%diag(ourMethodResult$psi[3,])%*%t(ourMethodResult$Q)

    
    allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
                                  naiveMeth=correlationResult[row(correlationResult) > col(correlationResult)],
                                  naiveWBatch=correlationWBatchResult[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  combat=combat[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  limma=limma[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  sva=sva[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  ruv=ruv[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  labels=trueLabels)
    if(maxPoints){
        numberOfPointsToPlot <- min(maxPoints, nrow(allCorrelations))
        allCorrelations <- allCorrelations[sample(nrow(allCorrelations),numberOfPointsToPlot),]
    }
    allCorrelations
}