# Replace target with your variable of interest
# tse is a treesummarizedexperiment table - some tests works with that best
# ps is phyloseq
# some output require refLevel, which is based on the reference level of your target variable

suppressPackageStartupMessages({
  library(mia)
  library(patchwork)
  library(tidySummarizedExperiment)
  library(ALDEx2)
  library(Maaslin2)
  library(MicrobiomeStat)
  library(knitr)
  library(tidyverse)
  library(ANCOMBC)
  library(phyloseq)
  library(microbiome)
  library(tidyverse)
})

# ALDEx2 
x <- aldex.clr(
  reads = assay(tse),
  conds = colData(tse)$target, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 1000, 
  denom = "all",
  verbose = FALSE
)
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  hist.plot = TRUE,
  verbose = TRUE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
aldex_out <- data.frame(x_tt, x_effect)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps)), Genus = tax_table(ps)[, "Genus"])

# Merge genus names with the aldex_out data frame
aldex_out_with_names <- merge(aldex_out, genus_info, by.x = "row.names", by.y = "Genus_ID")

# Set the row names of the merged data frame to the actual genus names
rownames(aldex_out_with_names) <- aldex_out_with_names$Genus

# Remove the 'Row.names' and 'Genus' columns
aldex_out_with_names <- aldex_out_with_names[, !(names(aldex_out_with_names) %in% c("Row.names", "Genus"))]

library(dplyr)
aldex_out_with_names %>%
  filter(wi.eBH <= 0.05) %>%
  dplyr::select(we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()

# ANCOM-BC
out = ancombc(phyloseq = ps, formula = "target",
              p_adj_method = "holm", group = "target", alpha = 0.05)
res = out$res
significant_taxa_indices <- which(res$diff_abn$targetrefLevel == TRUE)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(ps)), Genus = tax_table(ps)[, "Genus"])
significant_taxa_info <- genus_info[genus_info$Genus_ID %in% rownames(otu_table(ps)[significant_taxa_indices,]), ]
print(significant_taxa_info)

# MaAslin2
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
fit_data <- Maaslin2(
  asv,
  meta_data,
  output = "target",
  transform = "NONE",
  fixed_effects = "target",
  reference = "target,refLevel",  
  normalization = "CLR",
  standardize = FALSE
)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(Genus_target)), Genus = tax_table(Genus_target)[, "Genus"])

# Filter significant taxa
significant_taxa <- filter(fit_data$results, qval <= 0.05)

# Merge significant_taxa and genus_info to display names
significant_taxa_with_names <- merge(significant_taxa, genus_info, by.x = "feature", by.y = "Genus_ID")

# Remove the 'feature' column, which contains the full sequences
significant_taxa_with_names$feature <- NULL

# Display the table
kable(significant_taxa_with_names)

# LinDA
otu.tab <- as.data.frame(assay(tse))
meta <- as.data.frame(colData(tse)) %>% select(target)
res_linda <- linda(
  otu.tab,
  meta,
  formula = '~target',
  alpha = 0.05,
  prev.filter = 0.10)

# Extract genus names from the tax_table
genus_info <- data.frame(Genus_ID = rownames(tax_table(Genus_target)), Genus = tax_table(Genus_target)[, "Genus"])

# Convert LinDA output to a data frame and add genus names
linda_output_df <- as.data.frame(res_linda$output)
linda_output_with_names <- merge(linda_output_df, genus_info, by.x = "row.names", by.y = "Genus_ID")

# Make the Genus names unique
unique_genus_names <- make.unique(as.character(linda_output_with_names$Genus))

# Set row names as the unique Genus names
rownames(linda_output_with_names) <- unique_genus_names

# Remove the 'row.names' and 'Genus' columns
linda_output_with_names$Row.names <- NULL
linda_output_with_names$Genus <- NULL

# Filter significant results and display the table using kable
library(knitr)
significant_linda <- linda_output_with_names[linda_output_with_names$targetrefLevel.reject == TRUE, ]
kable(significant_linda)