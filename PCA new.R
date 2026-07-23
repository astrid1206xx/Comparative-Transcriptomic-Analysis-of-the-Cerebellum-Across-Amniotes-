install.packages(c("rotl", "ape", "phytools", "dplyr", "stringr"))

library(rotl)
library(ape)
library(phytools)
library(dplyr)
library(stringr)

# Species list mapping common names to scientific names
common_to_sci <- c(
  "Human"     = "Homo sapiens",
  "Chimpanzee"= "Pan troglodytes",
  "Bonobo"    = "Pan paniscus",
  "Gorilla"   = "Gorilla gorilla",
  "Orangutan" = "Pongo abelii",
  "Macaque"   = "Macaca mulatta",
  "Mouse"     = "Mus musculus",
  "Opossum"   = "Monodelphis domestica",
  "Platypus"  = "Ornithorhynchus anatinus",
  "Chicken"   = "Gallus gallus"
)

# Retrieve OpenTree data
sci <- unname(common_to_sci)
tn  <- tnrs_match_names(sci)
if (any(is.na(tn$ott_id))) {
  stop("Trouble with: ",
       paste(sci[is.na(tn$ott_id)], collapse=", "))
}

tr <- tol_induced_subtree(ott_ids = tn$ott_id)
tr <- multi2di(tr)
if (is.null(tr$edge.length)) tr <- compute.brlen(tr)

# Process tip labels
tip_full  <- tr$tip.label
tip_plain <- gsub("_ott[0-9]+$", "", tip_full)     
tip_pretty <- gsub("_", " ", tip_plain)            

# Map tip labels to common names
sci_plain <- gsub("_", " ", tip_plain)             
sci_to_common <- setNames(names(common_to_sci), common_to_sci)
tip_to_common <- sci_to_common[sci_plain]
tip_label_nice <- ifelse(is.na(tip_to_common), tip_pretty, tip_to_common)

# Plot tree (phylogram style)
plot(tr,
     show.tip.label = FALSE,
     type = "phylogram",
     direction = "rightwards",
     cex = 1.5,
     x.lim = c(0, max(nodeHeights(tr)) * 1.2), 
     label.offset = 0.02,                      
     no.margin = TRUE)
title("Phylogeny")
tiplabels(tip_label_nice,
          frame = "none",
          adj = -0.05,
          cex = 1.2)  

# Plot tree (standard style)
plot(tr, show.tip.label = FALSE, cex=1.5); title("Phylogeny")
tiplabels(tip_label_nice, frame="none", adj=-0.05, cex=1.5)

# Plot tree (fan style)
phytools::plotTree(tr, type="fan", ftype="off"); title("Phylogeny (fan)")
phytools::tiplabels(tip_label_nice, frame="none", adj=0.5, cex=0.9)

par(mfrow=c(1,1))

# Save tree and plot
write.tree(tr, file = "ten_species_opentree.newick")
pdf("ten_species_phylogeny.pdf", width=9, height=6)
plot(tr, show.tip.label = FALSE, cex=0.9); title("Ten-species phylogeny")
tiplabels(tip_label_nice, frame="none", adj=-0.05, cex=0.9)
dev.off()

# Print reference table
mapping <- data.frame(
  tip_label = tip_full,
  scientific_name = sci_plain,
  common_name = tip_label_nice,
  stringsAsFactors = FALSE
)
print(mapping, row.names = FALSE)

# Additional species processing
raw_names <- c(
  "Human"     = "Homo sapiens",
  "Chimpanzee"= "Pan troglodytes",
  "Bonobo"    = "Pan paniscus",
  "Gorilla"   = "Gorilla gorilla",
  "Orangutan" = "Pongo abelii",
  "Macaque"   = "Macaca mulatta",
  "Mouse"     = "Mus musculus",
  "Opossum"   = "Monodelphis domestica",
  "Platypus"  = "Ornithorhynchus anatinus",
  "Chicken"   = "Gallus gallus"
)

# Standardize names
raw_names <- raw_names |>
  stringr::str_trim() |>
  stringr::str_replace("^\\?+", "")

common_to_sci <- c(
  "human"       = "Homo sapiens",
  "chimpanzee"  = "Pan troglodytes",
  "gorilla"     = "Gorilla gorilla"
)

scientific_targets <- ifelse(
  tolower(raw_names) %in% names(common_to_sci),
  common_to_sci[tolower(raw_names)],
  raw_names
) |> unname()

matches <- tnrs_match_names(scientific_targets)

# Check for mismatches
if (any(is.na(matches$ott_id))) {
  warning("Trouble with: \n",
          paste(scientific_targets[is.na(matches$ott_id)], collapse = ", "))
}

ott_ids <- unique(na.omit(matches$ott_id))
tr <- tol_induced_subtree(ott_ids = ott_ids)
tr$tip.label <- gsub(" ", "_", tr$tip.label)

# Ensure binary tree and compute branch lengths if missing
tr <- multi2di(tr)
if (is.null(tr$edge.length)) tr <- compute.brlen(tr)
sci_clean <- gsub(" ", "_", matches$unique_name)

# Create alias mapping
alias_map <- tibble(
  input_alias = raw_names,
  resolved_scientific = matches$unique_name,
  tip = sci_clean
)

extra_alias <- tibble(
  input_alias = names(common_to_sci),
  resolved_scientific = unname(common_to_sci),
  tip = gsub(" ", "_", unname(common_to_sci))
) |>
  filter(tip %in% tr$tip.label)

alias_map <- bind_rows(alias_map, extra_alias) |>
  distinct(input_alias, .keep_all = TRUE)

alias_map <- bind_rows(
  alias_map,
  alias_map %>%
    transmute(input_alias = gsub(" ", "_", resolved_scientific),
              resolved_scientific, tip)
) |>
  distinct(input_alias, .keep_all = TRUE)

# Function to find tips by name
find_tip <- function(name, tree = tr, aliases = alias_map) {
  key <- tolower(gsub(" ", "_", str_trim(name)))
  aliases_lower <- aliases
  aliases_lower$input_alias <- tolower(gsub(" ", "_", aliases_lower$input_alias))
  hit <- aliases_lower %>% filter(input_alias == key)
  if (nrow(hit) == 0) {
    stop(sprintf("No names found: ", name))
  }
  tip <- hit$tip[1]
  idx <- which(tree$tip.label == tip)
  list(tip_label = tip, index = idx)
}

# Test the function
example_queries <- c("human", "Homo sapiens", "chimpanzee", "Pan troglodytes",
                     "gorilla", "Gorilla gorilla",
                     "Peromyscus maniculatus", "Cavia porcellus")

hits <- lapply(example_queries, find_tip)

# Plot tree with labels
pdf("tree_labels.pdf", width=8, height=6)
plot(tr,
     show.tip.label = TRUE,   
     cex = 1.2,               
     label.offset = 0.01,     
     no.margin = TRUE)
plot(tr, show.tip.label = FALSE)  
tiplabels(text = gsub("_", " ", tr$tip.label),
          adj = -0.1,
          frame = "none",
          cex = 1.2)  

tip_indices <- sapply(hits, function(x) x$index)
tiplabels(pch = 19, tip = unique(tip_indices), cex = 0.8)

print(alias_map)

# Install and load additional packages
options(repos = c(CRAN = "https://cran.rstudio.com"))
pkgs <- c("tidyverse", "matrixStats", "pheatmap", "viridis", "ggrepel", "uwot", "scales",
          "rotl", "phytools", "ape", "Matrix", "irlba", "stringr")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(pheatmap)
  library(viridis); library(ggrepel); library(uwot); library(scales)
  library(rotl); library(phytools); library(ape); library(Matrix); library(irlba); library(stringr)
})

# Read count data
counts <- read.table("C:/Users/Astrid Xu/Desktop/cerebellum/merged_all.txt",
                     sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, quote = "", comment.char = "")
# Replace NA with 0
counts[is.na(counts)] <- 0
counts <- as.matrix(counts); storage.mode(counts) <- "numeric"

# Process sample information
stopifnot(!is.null(colnames(counts)))
cn <- colnames(counts)
split_cn <- strsplit(cn, "\\.")
get_part <- function(i) vapply(split_cn, function(v) if(length(v) >= i) v[i] else NA_character_, "")
Species <- get_part(1)
Tissue  <- toupper(get_part(2))
Rep     <- suppressWarnings(as.integer(get_part(3)))

sample_info <- tibble(
  Sample  = cn,
  Species = Species,
  Tissue  = Tissue,
  Rep     = Rep
)

# Normalization functions
median_ratio_normalize <- function(counts_mat){
  geom_means <- exp(rowMeans(log(pmax(counts_mat,1)), na.rm=TRUE))
  geom_means[apply(counts_mat, 1, function(x) all(x==0))] <- NA
  ratios <- sweep(counts_mat, 1, geom_means, "/")
  sf <- apply(ratios, 2, function(x) median(x[is.finite(x) & !is.na(x) & x>0], na.rm=TRUE))
  sf[!is.finite(sf) | is.na(sf) | sf==0] <- 1
  norm <- sweep(counts_mat, 2, sf, "/")
  list(norm_counts = norm, size_factors = sf)
}

cpm_normalize <- function(counts_mat){
  libsize <- colSums(counts_mat); libsize[libsize==0] <- 1
  cpm <- sweep(counts_mat, 2, libsize, "/") * 1e6
  list(norm_counts = cpm, size_factors = libsize/mean(libsize))
}

# Apply normalization
norm_method <- "MRN"   
norm_res    <- if(norm_method=="MRN") median_ratio_normalize(counts) else cpm_normalize(counts)
norm_counts <- norm_res$norm_counts
logMat      <- log2(norm_counts + 1)

# Quality control
libsize   <- colSums(counts)
zero_frac <- colMeans(counts==0)

theme_clean <- function(base=12){
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color="grey85", linewidth=0.3),
          axis.title = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust=0.02),
          legend.position = "right")
}

df_lib <- tibble(Sample = names(libsize),
                 LibrarySize = as.numeric(libsize),
                 ZeroFrac = zero_frac) %>%
  left_join(sample_info, by="Sample")

print(
  ggplot(df_lib, aes(x=reorder(Sample, LibrarySize), y=LibrarySize, fill=Species)) +
    geom_col() + coord_flip() +
    scale_fill_viridis(discrete = TRUE, option="D") +
    labs(title="Library Size per Sample", x=NULL, y="Total Raw Counts") +
    theme_clean()
)

print(
  ggplot(df_lib, aes(x=reorder(Sample, ZeroFrac), y=ZeroFrac, fill=Species)) +
    geom_col() + coord_flip() +
    scale_fill_viridis(discrete = TRUE, option="D") +
    scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
    labs(title="Zero-count Fraction per Sample", x=NULL, y="% zeros (genes)") +
    theme_clean()
)

df_dens <- as_tibble(logMat, rownames="Gene") |>
  pivot_longer(-Gene, names_to="Sample", values_to="logExpr") |>
  left_join(sample_info, by="Sample")

print(
  ggplot(df_dens, aes(x=logExpr, group=Sample, color=Species)) +
    geom_density(linewidth=0.8, alpha=0.8) +
    scale_color_viridis(discrete=TRUE, option="D") +
    labs(title=paste0("Expression Density (log2 + 1, ", norm_method, ")"),
         x="log2(normalized + 1)", y="Density") +
    theme_clean()
)

# Filter lowly expressed genes
keep <- rowSums(counts >= 10) >= max(2, round(ncol(counts)*0.25))
logF <- logMat[keep, , drop=FALSE]

# PhyloPCA analysis
species_tag <- sample_info$Species 

tag_to_sci <- c(
  "Human"              = "Homo sapiens",
  "Chimpanzee"         = "Pan troglodytes",
  "Gorilla"            = "Gorilla gorilla",
  "Gibbon"             = "Hylobates lar",
  "Cattle"             = "Bos taurus",
  "Pig"                = "Sus scrofa",
  "Longtailed_Macaque" = "Macaca fascicularis",
  "Guinea_pig"         = "Cavia porcellus",
  "Naked_molerat"      = "Heterocephalus glaber",
  "Deermice1"          = "Peromyscus maniculatus",
  "Deermice2"          = "Peromyscus maniculatus"
)

if(!all(unique(species_tag) %in% names(tag_to_sci))){
  missing_tags <- setdiff(unique(species_tag), names(tag_to_sci))
  stop("Error: \n", paste(missing_tags, collapse=", "))
}

sci_vec <- unname(tag_to_sci[species_tag])
unique_sci <- sort(unique(sci_vec))

tn <- rotl::tnrs_match_names(unique_sci)
if(any(is.na(tn$ott_id))){
  bad <- unique_sci[is.na(tn$ott_id)]
  stop("Error: Could not resolve the following scientific names in OpenTree: \n", paste(bad, collapse=", "))
}
tr <- rotl::tol_induced_subtree(ott_ids = tn$ott_id)
tr$tip.label <- gsub(" ", "_", tr$tip.label)
tr <- multi2di(tr)
if (is.null(tr$edge.length)) tr <- compute.brlen(tr)

# Create mapping from scientific names to tree tips
sci_to_tip_map <- setNames(rep(NA_character_, length(unique_sci)), unique_sci)
for(sci in unique_sci) {
  sci_key <- tolower(gsub(" ", "_", sci))
  tip_key <- tolower(gsub(" ", "_", tr$tip.label))
  exact_match <- tr$tip.label[sci_key == tip_key]
  
  if (length(exact_match) == 1) {
    sci_to_tip_map[sci] <- exact_match
  } else {
    pref <- paste0("^", gsub(" ", "_", sci))
    prefix_match <- grep(pref, tr$tip.label, value = TRUE, ignore.case = TRUE)
    if (length(prefix_match) >= 1) {
      sci_to_tip_map[sci] <- prefix_match[1]
    }
  }
}

if(any(is.na(sci_to_tip_map))){
  stop("Error: Could not find tips for the following scientific names: \n",
       paste(names(sci_to_tip_map)[is.na(sci_to_tip_map)], collapse=", "))
}

sample_info$Tip <- sci_to_tip_map[tag_to_sci[sample_info$Species]]

message("Number of samples per tree tip:")
print(table(sample_info$Tip, useNA = "ifany"))
if(sum(!is.na(sample_info$Tip)) < 2L) stop("Error: Less than 2 samples mapped to the tree.")

# Filter and aggregate data
ok_idx <- which(!is.na(sample_info$Tip) & sample_info$Tip %in% tr$tip.label)
stopifnot(length(ok_idx) >= 2L)
logF_ok <- logF[, ok_idx, drop=FALSE]
si_ok   <- sample_info[ok_idx, , drop=FALSE]

species_mat <- sapply(split(as.data.frame(t(logF_ok)), si_ok$Tip),
                      function(df) colMeans(df, na.rm = TRUE))

X <- t(species_mat)

# Align with tree
common_tips <- intersect(tr$tip.label, rownames(X))
if (length(common_tips) < 2L) {
  cat("Debug info: \n")
  cat("  Tree tip examples: ", paste(head(tr$tip.label, 5), collapse=", "), "\n")
  cat("  X row name examples: ", paste(head(rownames(X), 5), collapse=", "), "\n")
  stop("Error: Insufficient overlap with phylogenetic tree (<2 species).")
}

tr_use <- ape::keep.tip(tr, common_tips)
X      <- X[tr_use$tip.label, , drop=FALSE]

nzv <- matrixStats::colVars(X) > 1e-8
X <- X[, nzv, drop=FALSE]
message(sprintf("Retained %d non-zero variance genes for PhyloPCA.", ncol(X)))

# Low-memory PhyloPCA implementation
top_p <- 2000L
if (ncol(X) > top_p) {
  v <- matrixStats::colVars(X, na.rm = TRUE)
  keep <- order(v, decreasing = TRUE)[seq_len(top_p)]
  X <- X[, keep, drop=FALSE]
  message(sprintf("Filtered to top %d highly variable genes for PhyloPCA.", top_p))
}

Xc <- scale(X, center = TRUE, scale = TRUE)

C <- ape::vcv(tr_use, corr = TRUE)
lambda_hat <- 1
try({
  subp <- min(1000L, ncol(Xc))
  set.seed(1)
  cols <- sample(seq_len(ncol(Xc)), subp)
  pp0 <- phytools::phyl.pca(tr_use, Xc[, cols, drop=FALSE], method = "lambda", mode = "corr")
  if (!is.null(pp0$lambda) && is.finite(pp0$lambda)) lambda_hat <- pp0$lambda
}, silent = TRUE)
message(sprintf("Estimated Pagel's λ = %.3f", lambda_hat))

n <- nrow(C)
C_lambda <- (1 - lambda_hat) * diag(n) + lambda_hat * C

L <- chol(C_lambda)
Y <- backsolve(L, Xc, transpose = FALSE)

K <- max(2, min(8, nrow(Y) - 1))
set.seed(123)
rp <- irlba::prcomp_irlba(Y, n = K, center = FALSE, scale. = FALSE)

var_total <- sum(matrixStats::colVars(Y)^2)
var_expl <- round(100 * (rp$sdev^2 / var_total), 1)

scores <- as.data.frame(rp$x)
scores$Tip <- rownames(scores)
rownames(scores) <- NULL

scores$Display_Name <- gsub("_ott[0-9]+$", "", scores$Tip)
scores$Display_Name <- gsub("_", " ", scores$Display_Name)

ggplot(scores, aes(PC1, PC2, label = Display_Name)) +
  geom_point(aes(color = Display_Name), size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = 50, show.legend = FALSE,
                           box.padding = 0.5, point.padding = 0.5) +
  scale_color_viridis(discrete = TRUE, option = "D", name = "Species") +
  labs(title = paste0("PhyloPCA (GLS-λ=", round(lambda_hat, 3), ", IRLBA, corr-mode)"),
       x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)")) +
  theme_clean() +
  theme(legend.position = "bottom")

# Create mapping from Tip to common names
name_map <- sample_info %>%
  filter(!is.na(Tip)) %>%
  select(Tip, General_Name = Species) %>%
  distinct(Tip, .keep_all = TRUE)

scores_with_names <- scores %>%
  left_join(name_map, by = "Tip")

print(
  ggplot(scores_with_names, aes(x = PC1, y = PC2, label = General_Name, color = General_Name)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text_repel(
      size = 3.5, 
      max.overlaps = 50, 
      show.legend = FALSE,
      box.padding = 0.5, 
      point.padding = 0.5
    ) +
    scale_color_viridis(discrete = TRUE, option = "D", name = "Species") +
    labs(
      title = paste0("PhyloPCA (GLS-λ=", round(lambda_hat, 3), ", IRLBA, corr-mode)"),
      x = paste0("PC1 (", var_expl[1], "%)"),
      y = paste0("PC2 (", var_expl[2], "%)")
    ) +
    theme_clean() +
    theme(legend.position = "bottom")
)

# Sample correlation heatmap
cor_mat <- cor(logF, method="pearson")
ann <- sample_info %>% column_to_rownames("Sample") %>%
  select(where(~!is.numeric(.)))
pal <- viridis::viridis(50)

pdf("sample_correlation_heatmap.pdf", width = 12, height = 12)
pheatmap(cor_mat,
         color = colorRampPalette(pal)(100),
         annotation_col = ann, annotation_row = ann,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = paste0("Sample–Sample Correlation (", norm_method, " + log2)"),
         fontsize_row = 6,
         fontsize_col = 6)
dev.off()

# MDS analysis
d <- as.dist(1 - cor_mat)
mds <- cmdscale(d, k = 2)
mds_df <- as_tibble(mds, .name_repair="minimal") %>% setNames(c("MDS1","MDS2")) %>%
  mutate(Sample = colnames(logF)) %>% left_join(sample_info, by="Sample")
ggplot(mds_df, aes(MDS1, MDS2, color=Species, shape=Tissue, label=Sample)) +
  geom_point(size=3, alpha=0.95) +
  ggrepel::geom_text_repel(size=3, max.overlaps=30, show.legend=FALSE) +
  scale_color_viridis(discrete=TRUE, option="D") +
  labs(title=paste0("MDS on 1 - Pearson (", norm_method, " + log2)")) +
  theme_clean()

# UMAP analysis
logF_scaled <- scale(t(logF))
set.seed(123)
um <- umap(logF_scaled, n_neighbors = min(15, ncol(logF_scaled)-1), min_dist = 0.2, metric = "cosine")
um_df <- as_tibble(um, .name_repair="minimal") %>% setNames(c("UMAP1","UMAP2")) %>%
  mutate(Sample = colnames(logF)) %>% left_join(sample_info, by="Sample")
ggplot(um_df, aes(UMAP1, UMAP2, color=Species, shape=Tissue, label=Sample)) +
  geom_point(size=3, alpha=0.95) +
  ggrepel::geom_text_repel(size=3, max.overlaps=30, show.legend=FALSE) +
  scale_color_viridis(discrete=TRUE, option="D") +
  labs(title=paste0("UMAP of Samples (", norm_method, " + log2)")) +
  theme_clean()