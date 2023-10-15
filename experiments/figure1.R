library(MASS)
library(ggplot2)
library(limma)
library(gridExtra)

setwd("/home/soel/cobra/figures/")
numSamples <- 1000
gene1 <- c(rnorm(numSamples,7),rnorm(numSamples,9))
gene2 <- c(rnorm(numSamples,7),rnorm(numSamples,9,2))
Batch <- c(rep("A",numSamples),rep("B",numSamples))

df11 <- data.frame(gene1,gene2,Batch, Correction="1. Raw Data", type="Example 1")

correct <- c(rep(0,numSamples),rep(2,numSamples))
scalecorrect <- c(rep(1,numSamples),rep(sqrt(2),numSamples))
gene1_ls <- gene1-correct
gene2_ls <- (gene2-correct-7)/scalecorrect+7
df12 <- data.frame(gene1=gene1_ls,gene2=gene2_ls,Batch, Correction="2. Batch Corrected",type="Example 1")

genes <- mvrnorm(n = numSamples, c(7,7), matrix(c(1,.95,.95,1),nrow=2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
gene1 <- c(genes[,1],rnorm(numSamples,9))
gene2 <- c(genes[,2],rnorm(numSamples,9))


df21 <- data.frame(gene1=gene1,gene2=gene2,Batch, Correction="1. Raw Data", type="Example 2")
df22 <- data.frame(gene1=gene1-correct,gene2=gene2-correct,Batch, Correction="2. Batch Corrected", type="Example 2")

df<- rbind(df11,df12,df21,df22)
ggplot(df, aes(x=gene1,y=gene2)) +geom_point(aes(col=Batch),size=3, alpha=.5) +
    ggtitle("Standard Batch Effect Correction") + xlab("Gene 1") + ylab("Gene 2") + 
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_bw(base_size = 40) + theme(plot.title = element_text(hjust = 0.5)) +  facet_grid(type ~Correction)
ggsave("figure1.pdf", plot = last_plot(), width = 6000, height = 6000, units = "px")

