# differential gene expression analysis of healthy vs NSCLC samples using DESeq2
# this script should be run after bulk_metadata_visualization.R

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/04a_design_formulas.html
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05a_hypothesis_testing.html

# Wald test is the default for two sample groups, could use likelihood ratio test (LRT) otherwise
dds_wald <- DESeq(dds, test="Wald")
#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/04b_DGE_DESeq2_analysis.html

plotDispEsts(dds_wald,main="Non-Small Cell Lung Cancer \nBulk RNA-Seq Wald Dispersion Estimates")
dev.copy2pdf(file="plots/Wald_test_dispersion.pdf")
#cloud of dispersion nvalues should follow the decreasing curve relative to mean,
# so this dataset appears to follow DESeq2 model assumptions

#Wald test results
#syntax: contrast <- c("condition", "level_to_compare", "base_level")
contrast<-c("sampletype","tumor","healthy")
res_table <- results(dds_wald, contrast=contrast, alpha = 0.05)#note this does padj, defaults to BH
res_table %>% data.frame() %>% View()#note padj in the results

# Save the unshrunken results to compare
res_table_unshrunken <- res_table
# Apply fold change shrinkage - shrinks fold changes for genes with low expression
res_table <- lfcShrink(dds_wald, coef="sampletype_tumor_vs_healthy", type="apeglm")
# MA plot using unshrunken fold changes
plotMA(res_table_unshrunken, ylim=c(-2,2),main="Unshrunken fold changes MA plot",returnData=TRUE)
dev.copy2pdf(file="plots/unshrunken_FC.pdf")
# MA plot using shrunken fold changes
plotMA(res_table, ylim=c(-2,2),main="Shrunken fold changes MA plot")
dev.copy2pdf(file="plots/clustering_heatmap.pdf")

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05c_summarizing_results.html
summary(res_table,alpha=0.05)#alpha 0.05 actually refers to FDR even though output says p value
#10% of genes have significantly increased expression in tumours, 9.6% have significantly
# decreased expression
# we'll identify genes annd perform GSEA/GO analysis next...