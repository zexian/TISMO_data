suppressMessages(library(ggplot2))
suppressMessages(library(ggfortify))
suppressMessages(library(RColorBrewer))
suppressMessages(library(optparse))
suppressMessages(library(ggpubr))

#label the cell line
Batch <- "Cell_Line"

#gene expression 
final.exp <- readRDS("")

#metasheet 
batch.dat <- read.csv("")

#PCA figure 
pca_p <- function(exprTable, annot,title, Batch) {
  batch_n <- length(unique(as.character(annot[colnames(exprTable),Batch])))
  df <- cbind.data.frame(t(exprTable),phenotype = as.character(annot[colnames(exprTable),Batch]))
  pca_plot <- autoplot(prcomp(t(exprTable)), data = df, col = 'phenotype', size = 1, frame = FALSE, frame.type = 'norm')+
    labs(title=title)+
    #   scale_color_manual(values = brewer.pal(name = "Set1", n = 9)[1:batch_n])+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"))
  return(pca_plot)
}


p2 <- pca_p(exprTable = final.exp, annot = batch.dat, title = "PCA plot After Batch Removal",Batch=Batch)
