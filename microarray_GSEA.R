library(affy)
library(limma)
#

Data <- ReadAffy(celfile.path = "GSE20051_RAW/")
eset <- rma(Data)
head(eset)
pData(eset)

strain <- c("DMSO","PLX","DMSO","PLX","DMSO","PLX","DMSO","PLX","PLX","DMSO","PLX","PLX")
design <- model.matrix(~factor(strain))
colnames(design) <- c("DMSO","PLX")
fit <- lmFit(eset, design)
fit <- eBayes(fit)
options(digits=2)
res<-topTable(fit, number=Inf, adjust.method="none", coef=1)# save your results
write.table(res,"dif_exp.txt",sep="\t")
##

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library("AnnotationDbi")
library(hgu133a2.db)
#
hs_hallmark_df <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "H"
)
hs_gocc_df <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C5", subcategory = "GO:CC"
)
hs_gobp_df <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C5", subcategory = "GO:BP"
)
#
hs_reatome_df <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C2", subcategory = "CP:REACTOME"
)
#
keytypes(org.Hs.eg.db)
#
ID     <- featureNames(eset)

out <- mapIds(hgu133a2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
#
# library("biomaRt")
# ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
# affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id")
# getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)

###
res$Gene <- rownames(res)
dge_mapped_df <- data.frame(
  entrez_id = mapIds(
    # Replace with annotation package for the organism relevant to your data
    hgu133a2.db,
    keys = res$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "PROBEID",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(res, by = c("Ensembl" = "Gene"))


dge_mapped_df2 <- dge_mapped_df[!duplicated(dge_mapped_df$entrez_id), ]

# #
# dup_entrez_ids <- dge_mapped_df %>%
#   dplyr::filter(duplicated(entrez_id)) %>%
#   dplyr::pull(entrez_id)
# 
# dup_entrez_ids
# dge_mapped_df %>%
#   dplyr::filter(entrez_id %in% dup_entrez_ids)

#
filtered_dge_mapped_df <- dge_mapped_df2 %>%
  # Sort so that the highest absolute values of the t-statistic are at the top
  dplyr::arrange(dplyr::desc(abs(t))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`-- this will keep
  # the first row with the duplicated value thus keeping the row with the
  # highest absolute value of the t-statistic
  dplyr::distinct(entrez_id, .keep_all = TRUE)

# Let's create a named vector ranked based on the t-statistic values
t_vector <- filtered_dge_mapped_df$t
names(t_vector) <- filtered_dge_mapped_df$entrez_id

# We need to sort the t-statistic values in descending order here
t_vector <- sort(t_vector, decreasing = TRUE)
head(t_vector)

gsea_results <- GSEA(
  geneList = t_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hallmark_df,
    gs_name,
    entrez_gene
  )
)
###
head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 3, order_by = NES)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  color.line = "#0d76ff"
)

most_positive_nes_plot
###

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(n = 3, order_by = NES)
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "HALLMARK_E2F_TARGETS",
  color.line = "#0d76ff"
)

most_negative_nes_plot
##########################################
######
gsea_results <- GSEA(
  geneList = t_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_gocc_df,
    gs_name,
    entrez_gene
  )
)
###
head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 3, order_by = NES)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOCC_CYTOSOLIC_RIBOSOME",
  title = "GOCC_CYTOSOLIC_RIBOSOME",
  color.line = "#0d76ff"
)

most_positive_nes_plot
###

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(n = 3, order_by = NES)
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOCC_ENDOCYTIC_VESICLE",
  title = "GOCC_ENDOCYTIC_VESICLE",
  color.line = "#0d76ff"
)

most_negative_nes_plot

# GOBP_CHAPERONE_MEDIATED_AUTOPHAGY


######
gsea_results <- GSEA(
  geneList = t_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_gobp_df,
    gs_name,
    entrez_gene
  )
)
###
head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 3, order_by = NES)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOBP_CYTOPLASMIC_TRANSLATION",
  title = "GOBP_CYTOPLASMIC_TRANSLATION",
  color.line = "#0d76ff"
)

most_positive_nes_plot
###

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(n = 3, order_by = NES)
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "GOBP_VIRAL_PROCESS",
  title = "GOBP_VIRAL_PROCESS",
  color.line = "#0d76ff"
)

most_negative_nes_plot


##
# filtered_and_selected_data <- gsea_results@result %>%
#   filter(ID == "GOBP_CHAPERONE_MEDIATED_AUTOPHAGY") 
# 
# 
# any(gsea_results@result$ID == "GOBP_CHAPERONE_MEDIATED_AUTOPHAGY")

#####E

######
gsea_results <- GSEA(
  geneList = t_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_reatome_df,
    gs_name,
    entrez_gene
  )
)
###
head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 3, order_by = NES)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  title = "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  color.line = "#0d76ff"
)

most_positive_nes_plot
###

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(n = 3, order_by = NES)
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "REACTOME_SELECTIVE_AUTOPHAGY",
  title = "REACTOME_SELECTIVE_AUTOPHAGY",
  color.line = "#0d76ff"
)

most_negative_nes_plot
any(gsea_results@result$ID == "REACTOME_SELECTIVE_AUTOPHAGY")

filtered_and_selected_data <- gsea_results@result %>%
  filter(ID == "REACTOME_SELECTIVE_AUTOPHAGY") 
