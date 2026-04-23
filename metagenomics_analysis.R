# shotgun_metagenomics_analysis.R - COMPLETE SELF-CONTAINED VERSION
# NO INTERNET DOWNLOAD REQUIRED - Works entirely offline

# ============================================
# STEP 1: Install and load packages
# ============================================

# Install packages if needed (only need to do once)
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("vegan")) install.packages("vegan")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("ggpubr")) install.packages("ggpubr")

# Load packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

# ============================================
# STEP 2: Create realistic simulated dataset
# ============================================

cat("\n========================================\n")
cat("SHOTGUN METAGENOMICS ANALYSIS PIPELINE\n")
cat("========================================\n\n")
cat("Creating simulated metagenomics dataset...\n")

set.seed(456)  # For reproducible results

# Dataset parameters
n_samples <- 30
n_taxa <- 40

# Create realistic taxonomic names (common gut microbiome taxa)
phyla <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", 
           "Verrucomicrobia", "Tenericutes", "Cyanobacteria", "Fusobacteria")

genera <- c("Bacteroides", "Prevotella", "Faecalibacterium", "Ruminococcus",
            "Escherichia", "Bifidobacterium", "Akkermansia", "Lactobacillus",
            "Clostridium", "Streptococcus", "Roseburia", "Alistipes",
            "Parabacteroides", "Coprococcus", "Dorea", "Blautia")

# Generate taxonomic names by combining phylum and genus
taxon_names <- paste0(sample(phyla, n_taxa, replace = TRUE), "_", 
                      sample(genera, n_taxa, replace = TRUE))
taxon_names <- make.unique(taxon_names)

# Create count matrix with realistic distribution
counts <- matrix(0, nrow = n_taxa, ncol = n_samples)

# Control samples (1-15) - baseline microbiome
for (i in 1:15) {
  # Negative binomial distribution typical for sequencing data
  counts[, i] <- rnbinom(n_taxa, size = 0.8, mu = 100 + rnorm(n_taxa, 0, 30))
}

# Treatment samples (16-30) - with perturbation
for (i in 16:30) {
  # Some taxa increase, some decrease
  fold_changes <- exp(rnorm(n_taxa, mean = 0, sd = 1.2))
  counts[, i] <- rnbinom(n_taxa, size = 0.8, 
                         mu = (100 + rnorm(n_taxa, 0, 30)) * fold_changes)
}

# Add zeros to simulate sparse data (typical in metagenomics)
counts[counts < 0] <- 0
zero_proportion <- 0.15  # 15% zeros
counts[sample(1:length(counts), size = length(counts) * zero_proportion)] <- 0

# Set row and column names
rownames(counts) <- taxon_names
colnames(counts) <- paste0("Sample_", 1:n_samples)

# Create metadata with clinical information
metadata <- data.frame(
  SampleID = colnames(counts),
  Group = factor(rep(c("Control", "Treatment"), each = 15),
                 levels = c("Control", "Treatment")),
  Patient = rep(paste0("P", 1:15), 2),
  Timepoint = rep(c("Baseline", "Follow-up"), each = 15),
  Age = round(runif(n_samples, 25, 70)),
  BMI = round(runif(n_samples, 18, 35), 1),
  Gender = sample(c("M", "F"), n_samples, replace = TRUE)
)

cat("✓ Dataset created successfully!\n")
cat(paste("  - Samples:", n_samples, "(15 Control, 15 Treatment)\n"))
cat(paste("  - Taxa:", n_taxa, "\n"))
cat(paste("  - Zero proportion:", round(mean(counts == 0) * 100, 1), "%\n\n"))

# ============================================
# STEP 3: Process and filter data
# ============================================

cat("STEP 1: Processing data...\n")

# Filter low-prevalence taxa (present in <10% of samples)
prevalence <- rowSums(counts > 0) / ncol(counts)
keep_taxa <- prevalence >= 0.1
filtered_counts <- counts[keep_taxa, ]
cat(paste("  - Kept", sum(keep_taxa), "out of", nrow(counts), "taxa\n"))

# Calculate relative abundance (percentage)
rel_abundance <- sweep(filtered_counts, 2, colSums(filtered_counts), "/") * 100

# Rarefaction to minimum library size for diversity analyses
min_lib_size <- min(colSums(filtered_counts))
rarefied <- rrarefy(t(filtered_counts), min_lib_size)
colnames(rarefied) <- rownames(filtered_counts)
rownames(rarefied) <- colnames(filtered_counts)
cat(paste("  - Rarefied to", min_lib_size, "reads/sample\n"))
cat(paste("  - Library sizes range:", 
          round(min(colSums(filtered_counts))), "-", 
          round(max(colSums(filtered_counts))), "\n\n"))

# ============================================
# STEP 4: Taxonomic composition plot
# ============================================

cat("STEP 2: Generating visualizations...\n")
cat("  - Taxonomic composition plot...\n")

# Get top 10 most abundant taxa
mean_abund <- rowMeans(rel_abundance)
top_taxa <- names(sort(mean_abund, decreasing = TRUE))[1:10]
top_abundance <- rel_abundance[top_taxa, ]
other_abundance <- colSums(rel_abundance[!rownames(rel_abundance) %in% top_taxa, ])

# Combine top taxa with "Other" category
plot_data <- rbind(top_abundance, Other = other_abundance)

# Prepare data for ggplot
plot_df <- as.data.frame(t(plot_data))
plot_df$Sample <- rownames(plot_df)
plot_df_long <- pivot_longer(plot_df, 
                             cols = -Sample, 
                             names_to = "Taxon", 
                             values_to = "Abundance")

# Add group information
plot_df_long <- left_join(plot_df_long, metadata[, c("SampleID", "Group")], 
                          by = c("Sample" = "SampleID"))

# Create stacked bar plot
p1 <- ggplot(plot_df_long, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Taxonomic Composition of Gut Microbiome",
       y = "Relative Abundance (%)",
       x = "Samples") +
  scale_fill_brewer(palette = "Set3") +
  facet_grid(~Group, scales = "free_x", space = "free")

print(p1)
ggsave("01_taxonomic_composition.png", p1, width = 14, height = 6, dpi = 300)
cat("    ✓ Saved: 01_taxonomic_composition.png\n")

# ============================================
# STEP 5: Alpha diversity analysis
# ============================================

cat("  - Alpha diversity analysis...\n")

# Calculate diversity metrics
shannon <- diversity(rarefied, index = "shannon")
simpson <- diversity(rarefied, index = "simpson")
richness <- specnumber(rarefied)
evenness <- shannon / log(richness)

diversity_df <- data.frame(
  Sample = names(shannon),
  Shannon = shannon,
  Simpson = simpson,
  Richness = richness,
  Evenness = evenness
)

diversity_df <- merge(diversity_df, metadata, by.x = "Sample", by.y = "SampleID")

# Shannon diversity boxplot with statistics
p2 <- ggplot(diversity_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Treatment Group",
       y = "Shannon Diversity (H')",
       x = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Treatment" = "#A23B72")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4)

print(p2)
ggsave("02_shannon_diversity.png", p2, width = 6, height = 5, dpi = 300)

# Richness plot
p3 <- ggplot(diversity_df, aes(x = Group, y = Richness, fill = Group)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  theme_minimal() +
  labs(title = "Species Richness by Treatment Group",
       y = "Number of Taxa (Observed Richness)",
       x = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Treatment" = "#A23B72")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4)

print(p3)
ggsave("03_richness.png", p3, width = 6, height = 5, dpi = 300)
cat("    ✓ Saved: 02_shannon_diversity.png, 03_richness.png\n")

# ============================================
# STEP 6: Beta diversity analysis
# ============================================

cat("  - Beta diversity analysis...\n")

# Calculate Bray-Curtis distance matrix
bc_dist <- vegdist(rarefied, method = "bray")

# Perform PCoA
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  PCo1 = pcoa$points[, 1],
  PCo2 = pcoa$points[, 2],
  Sample = rownames(pcoa$points)
)

# Add metadata
pcoa_df <- merge(pcoa_df, metadata, by.x = "Sample", by.y = "SampleID")

# Calculate variance explained
var_exp <- round(pcoa$eig[1:2] / sum(pcoa$eig) * 100, 1)

# Perform PERMANOVA
set.seed(123)
permanova <- adonis2(bc_dist ~ Group, data = metadata, permutations = 999)
p_value <- permanova$`Pr(>F)`[1]

# Create PCoA plot
p4 <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Beta Diversity Analysis - PCoA",
       subtitle = paste("Bray-Curtis distance, PERMANOVA p =", round(p_value, 3)),
       x = paste0("PCo1 (", var_exp[1], "%)"),
       y = paste0("PCo2 (", var_exp[2], "%)")) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Treatment" = "#A23B72")) +
  scale_shape_manual(values = c("Control" = 16, "Treatment" = 17))

print(p4)
ggsave("04_beta_diversity.png", p4, width = 8, height = 6, dpi = 300)
cat("    ✓ Saved: 04_beta_diversity.png\n")

# ============================================
# STEP 7: Heatmap of top variable taxa
# ============================================

cat("  - Generating heatmap...\n")

# Identify most variable taxa
taxon_var <- apply(rel_abundance, 1, var)
top_variable <- names(sort(taxon_var, decreasing = TRUE))[1:25]
heatmap_data <- rel_abundance[top_variable, ]
heatmap_data <- log10(heatmap_data + 0.01)  # Log transform for visualization

# Create annotation for samples
annotation_col <- data.frame(
  Group = metadata$Group,
  Age = metadata$Age,
  BMI = metadata$BMI
)
rownames(annotation_col) <- metadata$SampleID

# Define color palette
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

# Generate and save heatmap
png("05_abundance_heatmap.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(heatmap_data,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = list(
           Group = c(Control = "#2E86AB", Treatment = "#A23B72")
         ),
         main = "Top 25 Most Variable Taxa",
         fontsize_row = 8,
         fontsize_col = 8,
         color = heatmap_colors,
         border_color = NA,
         show_colnames = TRUE,
         angle_col = 45)
dev.off()

cat("    ✓ Saved: 05_abundance_heatmap.png\n")

# ============================================
# STEP 8: Differential abundance analysis
# ============================================

cat("  - Differential abundance analysis...\n")

# Split samples by group
control_samples <- metadata$SampleID[metadata$Group == "Control"]
treatment_samples <- metadata$SampleID[metadata$Group == "Treatment"]

# Calculate statistics for each taxon
diff_results <- data.frame(
  Taxon = rownames(rel_abundance),
  mean_control = NA,
  mean_treatment = NA,
  sd_control = NA,
  sd_treatment = NA,
  log2FC = NA,
  p_value = NA
)

for (i in 1:nrow(rel_abundance)) {
  control_abund <- rel_abundance[i, control_samples]
  treatment_abund <- rel_abundance[i, treatment_samples]
  
  # Calculate means and standard deviations
  diff_results$mean_control[i] <- mean(control_abund)
  diff_results$mean_treatment[i] <- mean(treatment_abund)
  diff_results$sd_control[i] <- sd(control_abund)
  diff_results$sd_treatment[i] <- sd(treatment_abund)
  
  # Wilcoxon rank-sum test (non-parametric)
  test <- wilcox.test(control_abund, treatment_abund)
  diff_results$p_value[i] <- test$p.value
  
  # Log2 fold change with pseudocount to avoid division by zero
  fc <- (mean(treatment_abund) + 0.01) / (mean(control_abund) + 0.01)
  diff_results$log2FC[i] <- log2(fc)
}

# Adjust p-values for multiple testing (FDR)
diff_results$adj_p_value <- p.adjust(diff_results$p_value, method = "fdr")
diff_results$significant <- diff_results$adj_p_value < 0.05

# Classify direction of change
diff_results$direction <- ifelse(diff_results$log2FC > 1 & diff_results$significant, 
                                 "Enriched in Treatment",
                                 ifelse(diff_results$log2FC < -1 & diff_results$significant,
                                        "Enriched in Control",
                                        "Not significant"))

# Sort by significance
diff_results <- diff_results[order(diff_results$adj_p_value), ]

# Save results to CSV
write.csv(diff_results, "06_differential_abundance.csv", row.names = FALSE)
cat("    ✓ Saved: 06_differential_abundance.csv\n")

# Create volcano plot
p5 <- ggplot(diff_results, aes(x = log2FC, y = -log10(adj_p_value))) +
  geom_point(aes(color = significant, size = significant), alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#D62828")) +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Abundance",
       subtitle = paste("FDR-adjusted p-value < 0.05 (", 
                        sum(diff_results$significant), 
                        " significant taxa)"),
       x = expression("Log"[2] ~ "Fold Change (Treatment vs Control)"),
       y = expression("-Log"[10] ~ "Adjusted P-value")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "blue", alpha = 0.5, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.5) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  guides(color = guide_legend("Significant (FDR < 0.05)"),
         size = "none")

# Add labels for top significant taxa
sig_taxa <- diff_results[diff_results$significant & abs(diff_results$log2FC) > 1.5, ]
if (nrow(sig_taxa) > 0) {
  p5 <- p5 + geom_text(data = head(sig_taxa, 8),
                       aes(label = Taxon), 
                       vjust = -0.5, 
                       hjust = 0.5,
                       size = 3,
                       check_overlap = TRUE)
}

print(p5)
ggsave("07_volcano_plot.png", p5, width = 8, height = 6, dpi = 300)
cat("    ✓ Saved: 07_volcano_plot.png\n")

# ============================================
# STEP 9: Summary statistics and report
# ============================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

# Calculate summary statistics
cat("SUMMARY STATISTICS:\n")
cat(paste("  Total samples analyzed:", n_samples, "\n"))
cat(paste("  Control samples:", sum(metadata$Group == "Control"), "\n"))
cat(paste("  Treatment samples:", sum(metadata$Group == "Treatment"), "\n"))
cat(paste("  Total taxa analyzed:", nrow(rel_abundance), "\n"))
cat(paste("  Significant taxa (adj. p < 0.05):", sum(diff_results$significant), "\n"))
cat(paste("  - Enriched in Treatment:", 
          sum(diff_results$direction == "Enriched in Treatment"), "\n"))
cat(paste("  - Enriched in Control:", 
          sum(diff_results$direction == "Enriched in Control"), "\n\n"))

# Diversity statistics
cat("DIVERSITY METRICS:\n")
cat(paste("  Mean Shannon diversity (Control):", 
          round(mean(diversity_df$Shannon[diversity_df$Group == "Control"]), 3), "\n"))
cat(paste("  Mean Shannon diversity (Treatment):", 
          round(mean(diversity_df$Shannon[diversity_df$Group == "Treatment"]), 3), "\n"))
cat(paste("  Mean Richness (Control):", 
          round(mean(diversity_df$Richness[diversity_df$Group == "Control"]), 1), "\n"))
cat(paste("  Mean Richness (Treatment):", 
          round(mean(diversity_df$Richness[diversity_df$Group == "Treatment"]), 1), "\n\n"))

# Top significant taxa
top5 <- head(diff_results[diff_results$significant, ], 5)
if (nrow(top5) > 0) {
  cat("TOP 5 DIFFERENTIALLY ABUNDANT TAXA:\n")
  for (i in 1:nrow(top5)) {
    cat(paste("  ", i, ". ", top5$Taxon[i], " - Log2FC: ", 
              round(top5$log2FC[i], 2), ", Adj p: ", 
              format(top5$adj_p_value[i], scientific = TRUE), "\n", sep = ""))
  }
  cat("\n")
}

# Files generated
cat("OUTPUT FILES GENERATED:\n")
cat("  ✓ 01_taxonomic_composition.png - Stacked bar plot\n")
cat("  ✓ 02_shannon_diversity.png - Alpha diversity boxplot\n")
cat("  ✓ 03_richness.png - Species richness violin plot\n")
cat("  ✓ 04_beta_diversity.png - PCoA ordination plot\n")
cat("  ✓ 05_abundance_heatmap.png - Clustered heatmap\n")
cat("  ✓ 06_differential_abundance.csv - Complete statistical results\n")
cat("  ✓ 07_volcano_plot.png - Volcano plot\n\n")

# Save complete results as R object
results_list <- list(
  counts = filtered_counts,
  relative_abundance = rel_abundance,
  rarefied = rarefied,
  metadata = metadata,
  diversity = diversity_df,
  differential_results = diff_results,
  pcoa = pcoa,
  permanova = permanova
)

saveRDS(results_list, "metagenomics_results.rds")
cat("Complete results saved to: metagenomics_results.rds\n")
cat("To reload in R: results <- readRDS('metagenomics_results.rds')\n\n")

cat("========================================\n")
cat("NEXT STEPS:\n")
cat("========================================\n")
cat("1. View plots in current directory\n")
cat("2. Open '06_differential_abundance.csv' in Excel\n")
cat("3. Run 'results <- readRDS('metagenomics_results.rds')' to explore data\n")
cat("4. Modify parameters and re-run for different scenarios\n\n")

cat("Analysis pipeline completed successfully!\n")