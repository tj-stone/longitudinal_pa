##PA RNASeq script for differential expression analysis

#Load packages, read sample metadata, load ensembl hg38 transcript database
library(tximport)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(EnhancedVolcano)
library(Rtsne)
library(ggplot2)
library(GSEABase)
library(dplyr)
library(msigdbr)
library(GSVA)
library(limma)
library(RColorBrewer)
library(gplots)
library(scales)
library(reshape)

samples <- read.csv("rec_pa_rnaseq_metadata_qc_trimmed.csv")
rownames(samples) <- samples$Sample
ensembl.hg38.tx2gene <- tx2gene <- read.csv("ensembl.hg38.tx2gene.csv")

#Read in all files (PA + Controls) with tximport, construct DESeq dataset object, set design to Pri_Long
#to allow construction of the dds object (won't be used for diff testing, only tSNE plot)
files <- file.path(getwd(), paste0("quants_PA/", samples$Sample, "_quant", "/", "quant.sf"))
txi <- tximport(files,
                type = "salmon",
                tx2gene = ensembl.hg38.tx2gene,
                ignoreTxVersion = T)
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~Pri_Long)

#Regularised log transformation prior to tSNE (minimises bias effect of genes with small counts, also
#normalises for library size)
rds <- rlog(dds)

#Extract rlog transformed count data and use median absolute deviation to find top 5,000 most 
#variable genes across cohort
counts <- assay(rds)
mads <- apply(counts, 1, mad)
top5k <- counts[rev(order(mads))[1:5000],]
top5k <- t(top5k)

#Calculate TSNE components prior to plotting. Seed for reproducibility.
set.seed(1)
tsne <- Rtsne(as.matrix(top5k),
              check_duplicates = FALSE,
              pca = TRUE,
              perplexity = 5,
              theta = 0.5,
              dims = 2)

#Plot TSNE and write to file, to reproduce publication figure change width 4, height 2.5 and remove geom_text
tiff("PA_longitudinal_plus_controls_TSNE.tiff",
     units = "in",
     width = 7.5,
     height = 5.5,
     res = 300)
qplot(tsne$Y[,1],
      tsne$Y[,2],
      xlim = c(-50,50),
      ylim = c(-75,110),
      size = I(2.5),
      color = as.factor(samples$Pri_Long)) +
  labs(title = "Longitudinal pairings by case",
       subtitle = "",
       caption = "",
       x = "TSNE-1",
       y = "TSNE-2",
       color = "Sample type") +
  scale_colour_manual(values = c('springgreen',  "indianred", "steelblue2")) +
  geom_text(aes(label = samples$Pair),
            hjust = -0.5,
            vjust = -0.5) +
  theme_classic()
dev.off()

#Subset to include only paired primary/longitudinal PA, or read in files fresh with tximport, construct
#DESeq dataset object with Pair as a blocking factor and Pri_Long (timepoint) as comparative.
samplesPairs <- samples[1:22, ]
filesPairs <- file.path(getwd(), paste0("quants_PA/", samplesPairs$Sample, "_quant", "/", "quant.sf"))
txiPairs <- tximport(filesPairs,
                     type = "salmon",
                     tx2gene = ensembl.hg38.tx2gene,
                     ignoreTxVersion = T)
ddsPairs <- DESeqDataSetFromTximport(txiPairs,
                                     colData = samplesPairs,
                                     design = ~ Pair + Pri_Long)

#Filter out genes with low counts (< 10 per sample) across cohort. Keep remainder
keepPairs <- rowSums(counts(ddsPairs)) >= 10
ddsPairs <- ddsPairs[keepPairs,]

#Relevel Pri_Long (longitudinal timepoint) variable so that the base value for comparison is Primary
ddsPairs$Pri_Long <- relevel(ddsPairs$Pri_Long, ref = "Primary")

#Perform differential expression analysis, extract results, count and subset total number of
#differentially expressed genes with a false discovery adjusted p-val < 0.1.
#count and display total number with negative (label: -1) and positive (label: 1) fold change
ddsPairs <- DESeq(ddsPairs)
resPairs <- results(ddsPairs)
sum(resPairs$padj < 0.1, na.rm=TRUE)
resPairsSig <- subset(resPairs, padj < 0.1)
table(sign(resPairsSig$log2FoldChange))

#Order by Wald statistic and write to file for pre-ranked gene set enrichment analysis
resPairsOrdered <- resPairs[order(resPairs$stat),]
write.csv(resPairsOrdered, "pa_longitudinal_genes_preranked_waldstat.csv")

#Plot volcano plot to visualise differential expression, use log2 fold change (value 1) and false 
#discovery adjusted p-val (value 0.1) visualiation for cut-offs. 
volcanoData <- resPairs[rev(order(resPairs$padj)),]

tiff("PA_Longitudinal_vs_Primary_volcano_padj_0.1.tiff",
     units = "in",
     width = 5,
     height = 5,
     res = 300)
EnhancedVolcano(volcanoData,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                ylim = c(0, 4.5),
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                ylab = bquote(~-Log[10] ~ italic(P.adj)),
                col=c("grey", "grey", "grey", "black"),
                colAlpha = 1,
                pCutoff = 0.1,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 6.0,
                border = "full")
dev.off()

#Convert significant DEGs to data frame and fetch HGNC symbols from biomaRt
resPairsSig <- as.data.frame(resPairsSig)
resPairsSig$id <- rownames(resPairsSig)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(resPairsSig)
HGNC_list <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   values = genes,
                   mart =  mart)

#Attach and match HGNC symbols to results data, mathching by ensembl gene id. Write to file.
resPairsSig <- left_join(resPairsSig, HGNC_list, by = c("id" = "ensembl_gene_id"))
rownames(resPairsSig) <- resPairsSig$id
resPairsSig$id <- NULL
write.csv("paired_DEGs_pa_longitudinal_padj_0.1.csv")


##GSVA from dds RNAseq expression object##

#Normalise deseq dataset by variance stabilising transformation
ddsNorm <- vst(ddsPairs)

#Retrieve normalised data from ddsNorm
vst_df <- assay(ddsNorm) %>%
  as.data.frame() %>% #Make into a data frame
  tibble::rownames_to_column("ensembl_id") #Move gene IDs into their own column

#Check for duplicated IDs (needs to be 0 to work)
sum(duplicated(vst_df$ensembl_id))

vst_matrix <- vst_df %>%
  #Move gene identifiers to row names
  tibble::column_to_rownames("ensembl_id") %>%
  #Convert object into a matrix
  as.matrix()

#Overwrite column names
colnames(vst_matrix) <- paste(samplesPairs$Pair, samplesPairs$Pri_Long, sep = "-")

#Load msigdbr and get C8 cell type category gene sets
c8_cellType_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C8"
)

#Convert gene sets into lists for compatibility with GSVA
c8_cellType_list <- split(
  c8_cellType_gene_sets$human_ensembl_gene, #Genes to be split into pathways
  c8_cellType_gene_sets$gs_name #Pathways to be made higher levels of the list
)

#Subset to isolate gene sets from the Descartes fetal single cell dataset
c8_list_descartes <- c8_cellType_list[grep("DESCARTES", names(c8_cellType_list))]

#Test gene set variability for msigdb c8 single cell gene sets across cohort.
#Method in gsva set to zscore to scale all samples together for comparison
gsva_results_c8_descartes <- gsva(
  vst_matrix,
  c8_list_descartes,
  method = "zscore",
  kcdf = "Gaussian",
  min.sz = 15,
  max.sz = 500,
  mx.diff = TRUE,
  verbose = FALSE
)

#Set up experimental design matrix prior to using limma to find pathways that are sig different
#between primary and longitudinal in the GSVA scores of their respective samples
design <- model.matrix(~ samplesPairs$Pair + factor(samplesPairs$Pri_Long, levels=c("Primary", "Longitudinal")))
colnames(design)[ncol(design)] <- "LongitudinalVsPrimary"
fit <- lmFit(gsva_results_c8_descartes, design)
fit <- eBayes(fit)

#Extract significant pathways 
sigPathways <- topTable(fit, 
                        coef = "LongitudinalVsPrimary",
                        number = Inf,
                        p.value = 0.05,
                        adjust = "BH")

#Subset GSVA results to pull out those with Primary-Longitudinal significant differences
topMatrixGSVA <- t(gsva_results_c8_descartes[rownames(sigPathways),])

#Scale z scores to re-center into -4/+4 range for cleaner plotting
topMatrixGSVA <- rescale(topMatrixGSVA, to = c(-4, 4))

#Subtract even numbered rows (Primary) from odd numbered rows (Longitudinal), to get Primary -> Longitudinal enrichment z-score changes
GSVA_diff_paired <- topMatrixGSVA[!c(TRUE,FALSE),] - topMatrixGSVA[c(TRUE,FALSE),]
rownames(GSVA_diff_paired) <- paste("Pair", gsub("\\-.*","",rownames(GSVA_diff_paired)))
GSVA_diff_paired <- as.data.frame(GSVA_diff_paired)
names(GSVA_diff_paired) <- c("fetal cerebrum microglia", "fetal eye microglia")
GSVA_diff_paired$pair <- rownames(GSVA_diff_paired)

#Reshape data for ggplot
GSVA_barplot_data <- melt(GSVA_diff_paired[,c("pair", "fetal cerebrum microglia", "fetal eye microglia")],id.vars = 1)

#Plot waterfall plot of longitudinal change in enrichment z-scores for each patient.
tiff("PA_GSVA_diff_pairs_waterfall.tiff",
     units = "in",
     width = 8,
     height = 3,
     res = 300)
ggplot(GSVA_barplot_data, 
       aes(x = reorder(pair, -value),
           y = value)) +
  geom_bar(aes(fill = variable),
           colour = "black",
           stat = "identity",
           position = "dodge") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5.5)) + 
  scale_fill_manual(values = c("steelblue4", "steelblue1")) + 
  labs(x = NULL,
       y = "Longitudinal enrichment\nz-score increase") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 14))
dev.off()
