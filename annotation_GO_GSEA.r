# gene identification, GO and GSEA analysis of NSCLC dataset after Wald test in DESeq2
# should be run after DESeq_Wald.r

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05c_summarizing_results.html
### extract differentially expressed genes
padj.cutoff <- 0.05
# Create a tibble of results
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


# extract just significant genes
sig <- res_table_tb %>%
  filter(padj < padj.cutoff)

new_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble() #tibble of metadata

#now join with gene name info from tx2gene
normalized_counts <- counts(dds_wald, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(tx2gene, by=c("gene" = "ensgene")) %>%
  left_join(res_table_tb,by=c("gene"="gene"))

dir.create("results")
write_csv(normalized_counts,"results/gene_results.csv")


#examine a significant gene of interest - note the human leukocyte antigen A is the
# second most significant and is very underexpressed in the tumours to aid in 
# immune evasion
#tx2gene[tx2gene$symbol=="HLA-A","ensgene"]
#

# Plot the  normalized counts, using the samplenames (rownames(d) as labels)
plot<-plotCounts(dds, gene="ENSG00000206503", intgroup="sampletype",returnData=TRUE) #can also not return data to plot plotcounts directly, but appearance isn't super flexible
ggplot(plot, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(plot))) + 
  theme_bw() +
  ggtitle(paste("HLA-A Non-Small Cell Lung Cancer Bulk RNA-Seq results, padj=",as.character(signif(min(normalized_counts[normalized_counts$symbol=="HLA-A","padj"]),4)))) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/HLA-A_normalized_counts.png")

# now make a heatmap with clustering of samples and genes
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
ggsave("plots/clustering_heatmap.png")


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
  ggtitle("NSCLC vs Healthy tissue Volcano plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("plots/Wald_volcanoplot.png")
# overexpression of SPP1 in lung cancer has been previously described, it seems 
# to promote invasion/ metastatis: see 
#https://cancerci.biomedcentral.com/articles/10.1186/s12935-022-02749-x


#now let's get useful information about differenntially regulated genes beyond
# just their names using AnnotationHub
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/genomic_annotation.html has list of lots of good annotation databases with R tools!
#annnotationhub lets you access multiple resources 
# Load libraries
library(AnnotationHub)
library(ensembldb)

# Connect to AnnotationHub
ah <- AnnotationHub()#another similar package often used is AnnotationDbi
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens#determine most recent version of annotations, right now AH104864
human_ens<-human_ens[["AH104864"]]

## how to extract info from AnnotationHub object:
# See gene-level information 
#genes(human_ens, return.type = "data.frame") %>% View()
# See transcript-level information
#transcripts(human_ens, return.type = "data.frame") %>% View()
# See exon-level information
#exons(human_ens, return.type = "data.frame") %>% View()

#here we'll just use the gene-level data
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_id %in% res_table_tb$gene)
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()#reduce to 1-1 mappings  by taking first identifier for each

#get just nonduplicated genes from above
annotations_ahb <- annotations_ahb[which(duplicated(annotations_ahb$gene_name) == FALSE), ]


#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/10_FA_over-representation_analysis.html
## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_table_tb, annotations_ahb, by=c("gene"="gene_id"))   
## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(res_ids$gene)

## Extract significant results
sig <- dplyr::filter(res_ids, padj < 0.05)
sig_genes <- as.character(sig$gene)


## Run GO overrepresentation analysis 
ego <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

write.csv(data.frame(ego), "results/GO_enrichment_results.csv")
#you can overrepresentation analysis on just overexpressed or underexpressed by
# first running something like sig_up <- dplyr::filter(res_ids, padj < 0.05 & log2FoldChange > 0)

## Dotplot 
dotplot(ego, showCategory=30,font.size=10)+ggtitle("GO Enrichment Dotplot for \n healthy vs NSCLC RNA-Seq")+
  theme(plot.title=element_text(hjust=0.5,size=rel(1.5)))
ggsave("plots/GO_enrichment_dotplot.png")

ego <- enrichplot::pairwise_termsim(ego)
## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50,cex_label_category=0.7)+
  ggtitle("Enrichment map of 50 most significant genes,\n healthy vs NSCLC RNA-Seq")+
  theme(plot.title=element_text(hjust=0.7,size=rel(1.5)))
ggsave("plots/GO_enrichment_map.png")


## now Gene Set Enrichment Analysis
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/11_FA_functional_class_scoring.html
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
## Extract the foldchanges and sort
foldchanges <- res_entrez$log2FoldChange
names(foldchanges) <- res_entrez$entrezid
foldchanges <- sort(foldchanges, decreasing = TRUE)

set.seed(431435)
library(stats)
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)


write.csv(gseaKEGG@result, "results/GSEA_KEGG.csv", quote=F)
#GSEA plot of hsa04110 i.e. cell cycle, the most significantly enriched term
#very strong enrichment in large log2fold changes, suggesting upregulation of cell cycle
gseaplot(gseaKEGG, geneSetID = 'hsa04110')+
  ggtitle("Cell cycle GO term Enrichment Plot, Healthy vs NSCLC RNA-Seq")+
  theme(plot.title=element_text(hjust=0.5,size=rel(1.5)))
ggsave("plots/cell_cycle_GSEA.png")

# we can visualize significant KEGG pathway genes with over/underexpression
# indicated using pathview

dir.create("gsea_pathview")
## Output images for all significant KEGG pathways
get_kegg_plots <- function(id) {
  pathview(gene.data = foldchanges,
           kegg.dir="gsea_pathview",
           pathway.id = id, 
           species = "hsa",
           limit = list(gene = 2, cpd = 1))
}

temp <- gseaKEGG_results$ID[gseaKEGG_results$ID!="hsa05206"]
for (sig_result in temp){
  get_kegg_plots(sig_result)
}

#these's a lot to comb through here, but hsa04110/cell cycle is a good place to
# start - notice overexpression of many cell cycle proteins and underexpression
# of tumour suppressors like Rb
