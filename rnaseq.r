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
txi$counts %>% View()

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01c_RNAseq_count_distribution.html
#distribution is similar to binomial or poisson, however in reality variance does not equal mean like in poisson
#negative binomial is a better fit when mean<variance which is generally the case
#evalute mean vs variance for control samples - this is for demonstration and is probably rarely done in practice
#hist(txi$counts,breaks=200,xlab="Raw Expression counts",ylab="Number of genes",main="RNA-seq salmon pseudocount distribution")
#data<-txi$counts
#mean_counts <- apply(data[,1:3], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
#variance_counts <- apply(data[,1:3], 1, var)
#df <- data.frame(mean_counts, variance_counts)

#ggplot(df) +
#  geom_point(aes(x=mean_counts, y=variance_counts)) + 
#  scale_y_log10(limits = c(1,1e9)) +
#  scale_x_log10(limits = c(1,1e9)) +
#  geom_abline(intercept = 0, slope = 1, color="red")+
#  title("Tumor samples variance vs mean")

#mean_counts <- apply(data[,4:6], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
#variance_counts <- apply(data[,4:6], 1, var)
#df <- data.frame(mean_counts, variance_counts)

#ggplot(df) +
#  geom_point(aes(x=mean_counts, y=variance_counts)) + 
#  scale_y_log10(limits = c(1,1e9)) +
#  scale_x_log10(limits = c(1,1e9)) +
#  geom_abline(intercept = 0, slope = 1, color="red")+
#  title("Control samples variance vs mean")

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
#below is also demostration of individual steps normally done automatically be DeSeq
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sex + sampletype)
#View(counts(dds))
#dds <- estimateSizeFactors(dds)
#sizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)
#write.table(normalized_counts, file="salmon_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/03_DGE_QC_analysis.html

### Transform counts with rlog - for data visualization only
rld <- rlog(dds, blind=TRUE)#blind means it doesn't take sample groups into account - unbiased
plotPCA(rld, intgroup="sampletype")
plotPCA(rld,intgroup="sex")
plotPCA(rld,intgroup="smoker")

rld_mat <- assay(rld)    
## "assay()" is part of the "SummarizedExperiment" package which is a DESeq2 dependency and is loaded with the DESeq2 library
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function### Load pheatmap package

### Plot heatmap using the correlation matrix and the metadata object
#display.brewer.all() will give all color maps if you want to try others
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20,main="Non-Small Cell Lung Cancer \nBulk RNA-Seq Hierarchical Clustering Heatmap")

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/04a_design_formulas.html
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05a_hypothesis_testing.html
# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ sex)
# or Wald test
dds_wald <- DESeq(dds, test="Wald")
plotDispEsts(dds_lrt,main="Non-Small Cell Lung Cancer \nBulk RNA-Seq LRT Dispersion Estimates")
plotDispEsts(dds_wald,main="Non-Small Cell Lung Cancer \nBulk RNA-Seq Wald Dispersion Estimates")

#Wald test results
#syntax: contrast <- c("condition", "level_to_compare", "base_level")
contrast<-c("sampletype","tumor","healthy")
res_table <- results(dds_wald, contrast=contrast, alpha = 0.05)#note this does padj, defaults to BH
res_table %>% 
  data.frame() %>% 
  View()#note padj in the results

#genes with all 0 counts - were excluded automatically before analysis as they had very low chance of being detected as differentially expressed due to all 0 counts
res_table[which(res_table$baseMean == 0),] %>% 
  data.frame() %>% 
  View()

#genes with an extreme outlier, apparently none
res_table[which(is.na(res_table$pvalue) & 
                    is.na(res_table$padj) &
                    res_table$baseMean > 0),] %>% 
  data.frame() %>% 
  View()

## Save the unshrunken results to compare
res_table_unshrunken <- res_table
# Apply fold change shrinkage - shrinks fold changes for genes with low expression
res_table <- lfcShrink(dds_wald, coef="sampletype_tumor_vs_healthy", type="apeglm")
# MA plot using unshrunken fold changes
plotMA(res_table_unshrunken, ylim=c(-2,2),main="Unshrunken fold changes MA plot")
# MA plot using shrunken fold changes
plotMA(res_table, ylim=c(-2,2),main="Shrunken fold changes MA plot")

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05c_summarizing_results.html
summary(res_table,alpha=0.05)#alpha 0.05 actually refers to FDR even though output says p value

### now extract differentially expressed genes
padj.cutoff <- 0.05
# Create a tibble of results
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the tibble to keep only significant genes
sig <- res_table_tb %>%
  filter(padj < padj.cutoff)
new_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

#nonw join with gene name info
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(tx2gene, by=c("gene" = "ensgene")) %>%
  left_join(res_table_tb,by=c("gene"="gene"))

write_csv(normalized_counts,"gene_results.csv")
#examine some significant genes of interest with good p values in more detail
#tx2gene[tx2gene$symbol=="HLA-A","ensgene"]
plot<-plotCounts(dds, gene="ENSG00000206503", intgroup="sampletype",returnData=TRUE) #can also not return data to plot plotcounts directly, but appearance isn't super flexible
#plot
# Plot the  normalized counts, using the samplenames (rownames(d) as labels)
ggplot(plot, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(plot))) + 
  theme_bw() +
  ggtitle("HLA-A Non-Small Cell Lung Cancer Bulk RNA-Seq results") +
  theme(plot.title = element_text(hjust = 0.5))

# now make heatmap
### Extract normalized expression for significant genes from the tumor and control samples
norm_sig <- normalized_counts%>% 
  filter(gene %in% sig$gene)  

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_sig[2:7], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta[colnames(meta)!="age"], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,main="Non-Small Cell Lung Cancer Bulk RNA-Seq \nSignificant Gene Hierarchical Clustering Z scores")

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction, for volcano plot
res_table_tb <- res_table_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

###get top 10 most significant genes for labelling volcano plot
## Add all the gene symbols as a column from the grch38 table using bind_cols()
res_table_tb <- bind_cols(res_table_tb, symbol=tx2gene$symbol[match(res_table_tb$gene, tx2gene$ensgene)])
## Create an empty column to indicate which genes to label
res_table_tb <- res_table_tb %>% mutate(genelabels = "")
## Sort by padj values 
res_table_tb <- res_table_tb %>% arrange(padj)
## Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
res_table_tb$genelabels[1:10] <- as.character(res_table_tb$symbol[1:10])

## Volcano plot
ggplot(res_table_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("NSCLC vs healthy tissue Volcano plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

#now get significant genes from lrt
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/08a_DGE_LRT_results.html
res_LRT <- results(dds_lrt)#padj here is BH correction
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)
# Compare to numbers we had from Wald test
nrow(sig)

# Obtain rlog values for those significant genes
#cluster_rlog <- rld_mat[sigLRT_genes$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
#clusters <- degPatterns(cluster_rlog, metadata = meta[,"sampletype",drop=FALSE], time = "sampletype", col=NULL)

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/genomic_annotation.html has list of lots of good annotation databases with R tools!
#annnotationhub lets you access multiple resources though
# Load libraries
library(AnnotationHub)
library(ensembldb)

# Connect to AnnotationHub
ah <- AnnotationHub()#another similar package often used is AnnotationDbi
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens#determine most recent version of annotations
human_ens<-human_ens[["AH104864"]]
# See gene-level information 
genes(human_ens, return.type = "data.frame") %>% View()
# See transcript-level information
transcripts(human_ens, return.type = "data.frame") %>% View()
# See exon-level information
exons(human_ens, return.type = "data.frame") %>% View()
#here we'll just use the gene-level data
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_id %in% res_table_tb$gene)
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()#reduce to 1-1 mappings  by taking first identifier for each

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)

# How many rows does annotations_ahb have?
annotations_ahb %>% nrow()

# Return only the non-duplicated genes using indices
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

# How many rows are we left with after removing?
annotations_ahb %>% nrow()
# Determine how many of the Entrez column entries are NA
which(is.na(annotations_ahb$entrezid)) %>%  length()


#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/10_FA_over-representation_analysis.html
## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_table_tb, annotations_ahb, by=c("gene"="gene_id"))   
## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(res_ids$gene)

## Extract significant results
sig <- dplyr::filter(res_ids, padj < 0.05)

sig_genes <- as.character(sig$gene)
## Run GO enrichment analysis 
ego <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "GO_enrichment_results.csv")

#you can also do GO analysis with just upregulated or just downregulated
## Extract upregulated genes
sig_up <- dplyr::filter(res_ids, padj < 0.05 & log2FoldChange > 0)
sig_up_genes <- as.character(sig_up$gene)
upregulated_go<-data.frame(enrichGO(gene = sig_up_genes, 
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE))

## Extract downregulated genes
sig_down <- dplyr::filter(res_ids, padj < 0.05 & log2FoldChange < 0)
sig_down_genes <- as.character(sig_down$gene)
downregulated_go<-data.frame(enrichGO(gene = sig_down_genes, 
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE))

#back to all GO results
## Dotplot 
dotplot(ego, showCategory=40,font.size=10)

## Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)
## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
foldchanges <- sig$log2FoldChange
names(foldchanges) <- sig$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
foldchanges <- ifelse(foldchanges > 1, 1, foldchanges)
foldchanges <- ifelse(foldchanges < -1, -1, foldchanges)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=foldchanges, 
         vertex.label.font=6)

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/11_FA_functional_class_scoring.html

## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid
## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
set.seed(123456)
library(stats)
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

# Write results to file
write.csv(gseaKEGG_results, "gsea_kegg.csv", quote=F)
gseaplot(gseaKEGG, geneSetID = 'hsa04110')#cell cycle set ID

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa04110",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))

## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = foldchanges, 
           pathway.id = gseaKEGG_results$ID[x], 
           species = "hsa",
           limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), 
           get_kegg_plots)
