#https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01b_DGE_setup_and_overview.html  also has link to info on using raw counts rather than pseudocounts
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)

#data from https://www.ebi.ac.uk/gxa/experiments/E-GEOD-81089/Experiment%20Design


## List all directories containing data  
samples <- list.files(path = "./salmon_results/salmon_results", full.names = T,pattern="SRR*")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./salmon_results/salmon_results/", "") %>% 
  str_replace(".fastq", "")
files

tx2gene <- read.delim("tx2gene_grch38_ens94.txt")
tx2gene

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM",ignoreTxVersion=T)
attributes(txi)

#txi$counts %>% View()


names_to_status=c("tumor","tumor","tumor","healthy","healthy","healthy") #tumour are non small cell lung cancer
names(names_to_status)=c("SRR3474721","SRR3474733","SRR3474738","SRR3475324","SRR3475325","SRR3475335")
ages=c(78,74,48,62,62,55)
sex<-c("M","F","M","F","F","M")
smoker<-c("non-smoker","ex-smoker","smoker","non-smoker","non-smoker","non-smoker")

meta <- data.frame(names_to_status)
meta["sex"]<-sex
meta["age"]<-ages
meta["smoker"]<-smoker

colnames(meta)<-c("sampletype","sex","age","smoker")
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/02_DGE_count_normalization.html

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sex + sampletype)

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/03_DGE_QC_analysis.html

### Transform counts with rlog - for data visualization only
rld <- rlog(dds, blind=TRUE)#blind means it doesn't take sample groups into account - unbiased

for (variable in c("sampletype","sex","smoker")){
  plotPCA(rld, intgroup=variable)+ggtitle(paste0("Non-Small Cell Lung Cancer PCA coloured by",variable))
  ggsave(paste0("plots/PCA_by_",variable,".png"))
}



rld_mat <- assay(rld)    
## "assay()" is part of the "SummarizedExperiment" package which is a DESeq2 dependency and is loaded with the DESeq2 library
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function### Load pheatmap package
### Plot heatmap using the correlation matrix and the metadata object
#display.brewer.all() will give all color maps if you want to try others
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
plot <- pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20,main="Non-Small Cell Lung Cancer RNA-Seq \nHierarchical Clustering Heatmap")
ggsave("plots/hierarchical_clustering.png",plot)
