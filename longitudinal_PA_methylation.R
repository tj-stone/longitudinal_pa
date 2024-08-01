##Methylation script clustering & differential methylation + KEGG enrichment
library(minfi)
library(DMRcate)
library(limma)
library(Rtsne)
library(missMethyl)
library(msigdbr)
library(ggplot2)

#Begin from RGSet containing longitudinal PA cases + 4x temporal lobe controls
RGSet

#Check proportion of failed probes on each array as a general QC metric for array performance
detP <- detectionP(RGSet)
failed <- detP > 0.01
failedProbes <- colMeans(failed)
failedProbes

#Extract SNP beta values for paired samples with built in getSnpBeta() command in minfi
#Note, control samples are not excluded here, they automatically drop out when correlation
#filters are applied later, and can be pulled out if desired as true 'unrelated' samples to
#give an idea of the range of correlation values unrelated samples result in.
snpBeta <- getSnpBeta(RGSet)

#Re label the columns with sample label pulled from RGSet metadata
colnames(snpBeta) <- paste(RGSet$Primary_Longitudinal, RGSet$Sample_Group, sep = "-")

#Get individual/patient/subject ID as factor
labelID <- as.factor(RGSet$Sample_Group)

#MDS plot to check broad similarity across samples
#Most cases cluster alongside the other member of their pair, exception is pair 13
mdsPlot(snpBeta,
        sampNames = labelID,
        xlim = c(-10, 10),
        ylim = c(-10, 10),
        pal = "black",
        main = "MDS projection on SNP probes, labelled by subject")

#Generate correlation matrix of SNP beta values to check similarity of profiles between samples
#Spearman correlation for non-parametric data. Alternatively Kendall, similar distributions; lower values
corr <- cor(snpBeta, method = "spearman")

#Convert to data frame, remove 1:1 corrs (samples correlated against themselves) and plot
#histogram to check distribution of correlation scores
#Comparing all samples vs all others, expect data to form two distributions:
#Comparisons from same pair/patient should have high corr and lay on the right of the histogram
#Comparisons from different pairs/poor quality samples should have lower corr and lay on the left.
#From histogram, Rho of 0.7 is suitable cut-off between matching pairs and unrelated samples
corr_df <- as.data.frame(as.table(corr))
corr_df <- subset(corr_df, corr_df$Var1 != corr_df$Var2)
hist(corr_df$Freq, breaks = 20, main = NULL, xlab = "Spearman's Rho")
abline(v = 0.7, col = "red")

#Subset to remove comparisons with Spearman corr lower than 0.7
#Also subset to remove samples compared against themselves
corr_df_high <- subset(corr_df, abs(Freq) > 0.7)
corr_df_high

#Sample pair 13 doesn't appear in high correlation subset, pull out pair 13 manually to confirm low correlation scores
subset(corr_df, corr_df$Var1 == "Primary-13" & corr_df$Var2 == "Longitudinal-13")

#Remove sample pair 13 from RGSet prior to further analysis
RGSet <- RGSet[,RGSet$Sample_Group != 13]

#Preprocess data with functional normalisation, remove CpGs on X and Y chromosomes, use built-in DMRcate
#rmSNPandCH function to remove cross-hybridising probes + SNPs + minor allele frequency (MAF).
MSet <- preprocessFunnorm(RGSet)
A <- getAnnotation(MSet)
remove <- rownames(A)[A$chr == "chrX" | A$chr == "chrY"]
MSet <- MSet[!rownames(MSet) %in% remove,]
MValues <- getM(MSet)
MValues.noSNPs <- rmSNPandCH(MValues, dist = 2, mafcut = 0.05)
MSet <- MSet[rownames(MValues.noSNPs),]

#Extract beta values and use median absolute deviation to find top 10,000 most variable probes across
#cohort.
beta <- getBeta(MSet)
beta <- na.omit(beta)
mads <- apply(beta, 1, mad)
top10k <- beta[rev(order(mads))[1:10000],]

#Calculate TSNE components prior to plotting. Seed for reproducibility.
set.seed(1)
tsne <- Rtsne(as.matrix(t(top10k)),
              check_duplicates = FALSE,
              pca = TRUE,
              perplexity = 10, # Set to 9 if repeating after removing outlier sample
              theta = 0.5,
              dims = 2)

#Plot TSNE and write to file, to reproduce publication figure change width 4, height 2.5 and remove geom_text
tiff("PA_control_plus_longitudinal_Meth_TSNE.tiff",
     units = "in",
     width = 7.5,
     height = 5.5,
     res = 300)
qplot(tsne$Y[,1],
      tsne$Y[,2],
      xlim = c(-10, 10),
      ylim = c(-10, 10),
      size = I(2.5),
      color = as.factor(MSet$Primary_Longitudinal)) +
  labs(title = "Longitudinal pairings by case",
       subtitle = "",
       caption = "",
       x = "TSNE-1",
       y = "TSNE-2",
       color = "Sample type") +
  scale_colour_manual(values = c('springgreen',  "indianred", "steelblue2")) +
  geom_text(aes(label = MSet$Sample_Group),
            hjust = -0.5,
            vjust = -0.5) +
  theme_classic()
dev.off()

#Primary sample from pair 9 clusters with controls, so remove pair from further analysis
MSet <- MSet[,MSet$Sample_Group != "9"]
#Can repeat tsne from above (repeat median absolute deviation and top10k selection, 
#also change perplexity to 9) to check clustering after removal. This will replicate
#the final methylation TSNE from manuscript

#Subset dataset to remove controls prior to  paired differential methylation testing. 
#Extract M (methyl) values.
MSetPA <- MSet[,MSet$Sample_Group != "Control"]
M <- getM(MSetPA)

#Convert experimental design variables to factor (patient, primary/longit) or numeric (time - continuous data).
#Relevel longitude factor variable to make primary the base value.
patient <- factor(MSetPA$Sample_Group)
time <- as.numeric(MSetPA$Lapsed_Time_Yr)
longitude <- as.factor(MSetPA$Primary_Longitudinal)
longitude <- relevel(longitude, ref = "Primary")

#Model matrix and differential methylation to find diff meth between primary and longitudinal timepoints
#(pairwise analysis of primary vs longitudinal used in manuscript)
designLongitude <- model.matrix(~patient + longitude)
fit <- lmFit(M, designLongitude)
fit <- eBayes(fit)
topTable(fit, coef = 'longitudeLongitudinal')
#Sum of top diff meth CpGs with FDR adjusted pvalue < 0.1 (none found)
sum(topTable(fit, coef = 'longitudeLongitudinal')$adj.P.Val < 0.1)
#Remove FDR adjustment and use dim to count num of CpGs that have pval < 0.05.
dim(topTable(fit, coef = 'longitudeLongitudinal', number = Inf, p.value = 0.05, adjust = "none"))
#Extract CpGs with pval < 0.05 and get corresponding gene names using Illumina annotation.
sigCpGs <- topTable(fit, coef = "longitudeLongitudinal", number = Inf, p.value = 0.05, adjust = "none")
cpgsLongitude <- rownames(sigCpGs)
GeneSymbols <- A$UCSC_RefGene_Name[A$Name %in% cpgsLongitude]

#Count number of CpGs with positive and negative fold change (relative to Primary samples)
dim(sigCpGs[which(sigCpGs$logFC > 0), ])
dim(sigCpGs[which(sigCpGs$logFC < 0), ])

#Split significantly diff meth probes into Up- and Down-methylated
cpgsLongitudeUp <- rownames(sigCpGs[which(sigCpGs$logFC > 0), ])
cpgsLongitudeDn <- rownames(sigCpGs[which(sigCpGs$logFC < 0), ])

#KEGG enrichment analysis of longitudinal diff meth pval < 0.05 CpGs using missMethyl
all.cpgs <- rownames(M)
KEGG_all <- gometh(cpgsLongitude, all.cpg = all.cpgs, collection = "KEGG", array.type = "EPIC")
KEGG_all <- KEGG_all[order(KEGG_all$FDR),]
KEGG_sig <- KEGG_all[KEGG_all$FDR < 0.1,]
KEGG_sig

#KEGG enrichment of longitudinal diff meth pval < 0.05, split into hyper (up) and hypo (dn) meth
KEGG_all_up <- gometh(cpgsLongitudeUp, all.cpg = all.cpgs, collection = "KEGG", array.type = "EPIC")
KEGG_all_up <- KEGG_all_up[order(KEGG_all_up$FDR),]
KEGG_sig_up <- KEGG_all_up[KEGG_all_up$FDR < 0.1,]
KEGG_sig_up

KEGG_all_dn <- gometh(cpgsLongitudeDn, all.cpg = all.cpgs, collection = "KEGG", array.type = "EPIC")
KEGG_all_dn <- KEGG_all_dn[order(KEGG_all_dn$FDR),]
KEGG_sig_dn <- KEGG_all_dn[KEGG_all_dn$FDR < 0.1,]
KEGG_sig_dn

#Additional enrichment analyses using gene sets from molecular signatures database
#Pull hallmark gene sets data from molecular signatures database
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H"
)

#Convert hallmark gene sets into lists
hallmarks_list <- split(
  hallmark_gene_sets$human_entrez_gene,
  hallmark_gene_sets$gs_name
)

#Overrepresentation analysis with gsameth function from missMethyl, subset by FDR pval 0.1 for significance
Hallmark_all <- gsameth(cpgsLongitude, all.cpg = all.cpgs, collection = hallmarks_list, array.type = "EPIC")
Hallmark_all <- Hallmark_all[order(Hallmark_all$FDR),]
Hallmark_sig <- Hallmark_all[Hallmark_all$FDR < 0.1,]
Hallmark_sig

#Pull C8 single cell gene sets data from molecular signatures database
c8_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C8"
)

#Convert genesets into lists
c8_list <- split(
  c8_gene_sets$human_entrez_gene,
  c8_gene_sets$gs_name
)

#Subset to isolate gene sets from the Descartes fetal single cell dataset
c8_list_descartes <- c8_list[grep("DESCARTES", names(c8_list))]

#Overrepresentation analysis with gsameth function from missMethyl, subset by FDR pval 0.1 for significance
#Note, can change 'c8_list_descartes' to 'c8_list' to assay all 700 C8 gene sets
c8_all <- gsameth(cpgsLongitude, all.cpg = all.cpgs, collection = c8_list_descartes, array.type = "EPIC")
c8_all <- c8_all[order(c8_all$FDR),]
c8_sig <- c8_all[c8_all$FDR < 0.1,]
c8_sig

#Overrepresentation using only CpGs with increased methylation in longitudinal compared to primary
c8_all_up<- gsameth(cpgsLongitudeUp, all.cpg = all.cpgs, collection = c8_list_descartes, array.type = "EPIC")
c8_all_up <- c8_all_up[order(c8_all_up$FDR),]
c8_sig_up <- c8_all_up[c8_all_up$FDR < 0.1,]
c8_sig_up

#Overrepresentation using only CpGs with reduced methylation in longitudinal compared to primary
c8_all_dn<- gsameth(cpgsLongitudeDn, all.cpg = all.cpgs, collection = c8_list_descartes, array.type = "EPIC")
c8_all_dn <- c8_all_dn[order(c8_all_dn$FDR),]
c8_sig_dn <- c8_all_dn[c8_all_dn$FDR < 0.1,]
c8_sig_dn
