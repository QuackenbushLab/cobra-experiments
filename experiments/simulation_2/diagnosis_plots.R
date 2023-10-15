library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)
library(rARPACK)
library(grid)
library(gridExtra)
library(reshape2)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Creates all the diagnosis plots for a result object
diagnosticPlots <- function(differentialCorrelationsDF, path){
    
    plotOursVsNaive <- ggplot(differentialCorrelationsDF[differentialCorrelationsDF$labels!="Background",]) + 
        geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.5, size=2) +
        ggtitle("Pairwise Differential Coexpression Estimates") + ylab("Our Method") + xlab("Naive Approach") + theme_bw(base_size = 60) +
        guides(color=guide_legend(title="True Interaction", override.aes = list(size=10), keyheight=1, keywidth=1, default.unit="inch")) + 
        scale_colour_manual(values=cbPalette[-1]) 
    
    png(paste0(path,'/OursVsNaive.png'), width = 1600, height = 1200)
    print(plotOursVsNaive)
    dev.off()
    
    plotOursVsNaive <- ggplot(differentialCorrelationsDF[differentialCorrelationsDF$labels!="Background",]) + 
        geom_point(aes(x=naiveWBatch,y=newMeth, color=factor(labels)), alpha=.5, size=2) +
        ggtitle("Pairwise Differential Coexpression Estimates") + ylab("Our Method") + xlab("Naive with Batch Approach") + theme_bw(base_size = 60) +
        guides(color=guide_legend(title="True Interaction", override.aes = list(size=10), keyheight=1, keywidth=1, default.unit="inch")) + 
        scale_colour_manual(values=cbPalette[-1]) 
    
    png(paste0(path,'/OursVsNaiveWBatch.png'), width = 1600, height = 1200)
    print(plotOursVsNaive)
    dev.off()
    
    naiveDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(naiveMeth, group=factor(labels), color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
        ggtitle("Naive Method") + xlab("Absolute Pairwise Differential Coexpression") +
        guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
        scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    naiveWBatchDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(naiveWBatch, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
        ggtitle("Naive Method with Batch") + xlab("Pairwise Differential Coexpression") +
        guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
        scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    ourMethodDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(newMeth, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
        ggtitle("COBRA") + xlab("Pairwise Differential Coexpression") +
        guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
        scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    limmaDensity <- ggplot(differentialCorrelationsDF) + 
      geom_density(aes(limma, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
      theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
      ggtitle("Limma") + xlab("Pairwise Differential Coexpression") +
      guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
      scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    combatDensity <- ggplot(differentialCorrelationsDF) + 
      geom_density(aes(combat, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
      theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
      ggtitle("Combat") + xlab("Pairwise Differential Coexpression") +
      guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
      scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    svaDensity <- ggplot(differentialCorrelationsDF) + 
      geom_density(aes(sva, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
      theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
      ggtitle("Sva") + xlab("Pairwise Differential Coexpression") +
      guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
      scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    ruvDensity <- ggplot(differentialCorrelationsDF) + 
      geom_density(aes(combat, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
      theme_bw(base_size = 60) + theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
      ggtitle("RUVCorr") + xlab("Pairwise Differential Coexpression") +
      guides(fill=guide_legend(title="True Interaction", keyheight=1, keywidth=1, default.unit="inch"), color=F) + 
      scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    png(paste0(path,'/NaiveDensity.png'), width = 1600, height = 1200)
    print(naiveDensity)
    dev.off()
    
    png(paste0(path,'/NaiveWBatchDensity.png'), width = 1600, height = 1200)
    print(naiveWBatchDensity)
    dev.off()
    
    png(paste0(path,'/ourMethodDensity.png'), width = 1600, height = 1200)
    print(ourMethodDensity)
    dev.off()
    
    png(paste0(path,'/limmaDensity.png'), width = 1600, height = 1200)
    print(limmaDensity)
    dev.off()
    
    png(paste0(path,'/combatDensity.png'), width = 1600, height = 1200)
    print(combatDensity)
    dev.off()
    
    png(paste0(path,'/svaDensity.png'), width = 1600, height = 1200)
    print(svaDensity)
    dev.off()
    
    png(paste0(path,'/ruvcorrDensity.png'), width = 1600, height = 1200)
    print(ruvDensity)
    dev.off()
    
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
        
        legend("bottomright", c(paste("Our Method",round(auc.methodPred,4)), 
                                paste("Naive",round(auc.methodPred.naive,4)), 
                                paste("Naive Batch",round(auc.methodPred.naive.w.batch,4)),
                                paste("Limma",round(auc.methodPred.limma,4)),
                                paste("ComBat",round(auc.methodPred.combat,4)),
                                paste("SVA",round(auc.methodPred.sva,4)),
                                paste("RUVCorr",round(auc.methodPred.ruv,4))),
               lty=1,lwd=1,col=c("red", "green", "blue", "pink", "purple", "darkgreen", "gray"),title="Area under ROC curve")
        abline(0,1)
    }
    unique(differentialCorrelationsDF$labels)
    png(paste0(path,'/OursVsOthersROC.png'), width = 800, height = 400)
    par(mfrow=c(1,2))
    plotROC(differentialCorrelationsDF, positive=c("Real effect","Negative effect"), "Real Effects vs All Pairs")
    plotROC(onlyEffects, positive=c("Real effect","Negative effect"), "Real Effects vs Batch Effects")
    dev.off()
    
}

plotEigenvectors <- function(ourMethodResult, trueLabels, path, numEigenvectors=6){
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    numVectors <- min(ncol(ourMethodResult$psi),20)
    numGenes <- nrow(ourMethodResult$G)
    nrows<-nrow(ourMethodResult$psi)
    
    plottingDF <- data.frame(ourMethodResult$Q[,seq_len(numEigenvectors)], trueLabels)
    names(plottingDF)[seq_len(numEigenvectors)] <- paste("Eigenvector",seq_len(numEigenvectors))
    plottingDFMelt <- melt(plottingDF, id.vars="trueLabels")
    plottingDFMelt$trueLabels <- trueLabels
    plottingDFMelt$Gene <- seq_len(numGenes)
    
    eigenvectorPlots <- ggplot(plottingDFMelt) + geom_point(aes(x=Gene,y=value,color=trueLabels)) +
        facet_wrap(~variable) + theme_bw(base_size = 30) + theme(legend.position="right") +
        guides(color=guide_legend(title="Gene Group:",override.aes = list(size=10)), title.position="top", title.hjust =0.5) + 
        scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette) + xlab("Genes")
    
    est <- t(ourMethodResult$psi[,1:numVectors])
    colnames(est)[1:3] <- c("Intercept","Batch","Case-Control") 
    eigenvalueDF <- melt(est,value.name = "Eigenvalue")
    eigenvaluePlots <- ggplot(eigenvalueDF) + geom_point(aes(x=Var1,y=Eigenvalue),size=4) + 
        facet_wrap(~Var2) + theme_bw(base_size = 30) +
        ylab("Eigenvalue Coefficent") + xlab("Top 20 Eigenvectors") + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
    
    png(paste0(path,'/EigenvectorPlots.png'), width = 1600, height = 1200)
    multiplot(eigenvectorPlots,eigenvaluePlots)
    dev.off()
}




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}