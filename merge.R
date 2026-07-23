# Install and load required packages
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
library(biomaRt)
library(dplyr)
library(tidyr)

# Read input file
infile <- "C:/Users/Astrid Xu/Desktop/cerebellum/hg_cb_raw_counts.txt"
outfile <- "C:/Users/Astrid Xu/Desktop/cerebellum/o_hg_cb_raw_counts.txt"
df <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- names(df)[1]

# Connect to Ensembl mart
host_used <- "https://may2024.archive.ensembl.org"
mart_hs <- useMart("ENSEMBL_MART_ENSEMBL", host = host_used, dataset = "hsapiens_gene_ensembl")

# Get available attributes
ats <- listAttributes(mart_hs)$name
name_attrs <- ats[grepl("^hglaber.*_homolog_associated_gene_name$", ats)]
ensg_attrs <- ats[grepl("^hglaber.*_homolog_ensembl_gene$", ats)]
type_attrs <- ats[grepl("^hglaber.*_homolog_orthology_type$", ats)]

# Get mapping data
attr_list <- unique(c("ensembl_gene_id", name_attrs, ensg_attrs, type_attrs))
raw <- getBM(attributes = attr_list, mart = mart_hs)

# Standardize gene symbols using HGNChelper
if (!requireNamespace("HGNChelper", quietly = TRUE)) install.packages("HGNChelper")
library(HGNChelper)

x <- unique(na.omit(df[[1]]))
fix <- HGNChelper::checkGeneSymbols(x, species = "human")
sym_std <- ifelse(is.na(fix$Suggested.Symbol), fix$x, fix$Suggested.Symbol)
df[[1]] <- sym_std

# Process chimpanzee data
infile <- "C:/Users/Astrid Xu/Desktop/cerebellum/cp_cb_raw_counts.txt"
if (!requireNamespace("stringi", quietly = TRUE)) install.packages("stringi")
library(stringi)

# Create mapping table
map_unique <- raw %>%
  select(ensembl_gene_id, hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# Add ENSG column to output
key <- names(df)[1]
df_out <- df %>%
  left_join(map_unique, by = setNames("hgnc_symbol", key)) %>%
  relocate(ensembl_gene_id, .before = all_of(key)) %>%
  select(ensembl_gene_id, all_of(names(df)))

write.table(df_out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Process mouse to human mapping
infile <- "C:/Users/Astrid Xu/Desktop/cerebellum/pm-po_CB_rawcounts.txt"
outfile <- "C:/Users/Astrid Xu/Desktop/cerebellum/o_pmpo_cb_raw_counts.txt"

# Connect to mouse dataset
mart_mm <- useMart("ENSEMBL_MART_ENSEMBL",
                   host = "https://jul2024.archive.ensembl.org",
                   dataset = "mmusculus_gene_ensembl")

# Read data
df <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
key <- names(df)[1]
sym_raw <- trimws(as.character(df[[1]]))

# Map mouse symbols to ensembl IDs
map_mouse1 <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "mgi_symbol",
  values = unique(sym_raw),
  mart = mart_mm
)

leftover <- setdiff(unique(sym_raw), unique(map_mouse1$mgi_symbol))
map_mouse2 <- if (length(leftover)) {
  getBM(
    attributes = c("ensembl_gene_id", "external_synonym"),
    filters = "external_synonym",
    values = leftover,
    mart = mart_mm
  )
} else {
  data.frame(ensembl_gene_id = character(), external_synonym = character())
}

# Combine mappings
map_mouse <- bind_rows(
  map_mouse1 %>% transmute(query = mgi_symbol, ensembl_gene_id),
  map_mouse2 %>% transmute(query = external_synonym, ensembl_gene_id)
) %>%
  filter(!is.na(query), query != "", !is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
  distinct(query, .keep_all = TRUE)

# Map to human homologs
map_hom <- getBM(
  attributes = c("ensembl_gene_id",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_perc_id"),
  filters = "ensembl_gene_id",
  values = unique(map_mouse$ensembl_gene_id),
  mart = mart_mm
)

# Combine and select best matches
map_all <- map_mouse %>%
  left_join(map_hom, by = "ensembl_gene_id")

map_best <- map_all %>%
  filter(!is.na(hsapiens_homolog_ensembl_gene),
         hsapiens_homolog_ensembl_gene != "") %>%
  mutate(.prio = ifelse(hsapiens_homolog_orthology_type == "ortholog_one2one", 1L, 2L),
         .pid = suppressWarnings(as.numeric(hsapiens_homolog_perc_id))) %>%
  arrange(query, .prio, desc(.pid), hsapiens_homolog_ensembl_gene) %>%
  distinct(query, .keep_all = TRUE) %>%
  transmute(query,
            human_ensembl = hsapiens_homolog_ensembl_gene,
            human_symbol = hsapiens_homolog_associated_gene_name)

# Add human ENSG column to output
df_out <- df %>%
  mutate(.query = sym_raw) %>%
  left_join(map_best, by = c(".query" = "query")) %>%
  relocate(human_ensembl, .before = all_of(key)) %>%
  select(human_ensembl, all_of(names(df))) %>%
  select(-.query)

write.table(df_out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Merge all count files
if (!requireNamespace("tools", quietly = TRUE)) install.packages("tools")
library(tools)

files <- c(
  "C:/Users/Astrid Xu/Desktop/cerebellum/primates_CB_rawcounts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_bt_cb_rawcounts_merged.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_ss_cb_rawcounts_merged.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_mf_cb_rawcounts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_cp_cb_raw_counts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_hg_cb_raw_counts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_pmpo_cb_raw_counts.txt"
)

# Function to read and clean each file
read_one <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  key <- names(df)[1]
  df <- df |> filter(!is.na(.data[[key]]), .data[[key]] != "")
  cnt_cols <- setdiff(names(df), key)
  df[cnt_cols] <- lapply(df[cnt_cols], function(x) suppressWarnings(as.numeric(x)))
  df <- df |>
    group_by(.data[[key]]) |>
    summarise(across(all_of(cnt_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  names(df)[1] <- "gene"
  df
}

# Read all files
list_of_dfs <- lapply(files, read_one)

# Merge all files
merged <- Reduce(function(x, y) full_join(x, y, by = "gene"), list_of_dfs)
merged <- merged |>
  mutate(across(-gene, ~ replace(., is.na(.), 0))) |>
  arrange(gene)

# Write merged output
out_file <- "C:/Users/Astrid Xu/Desktop/cerebellum/merged_counts_all.txt"
write.table(merged, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Read and check merged file
infile <- "C:/Users/Astrid Xu/Desktop/cerebellum/merged_all.txt"
merged <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
nrow(merged)
ncol(merged)
dim(merged)

# Function to calculate alignment statistics
align_stats_one <- function(path) {
  df <- read.table(path, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE,
                   quote = "", comment.char = "")
  
  col_human <- names(df)[1]
  col_species <- names(df)[2]
  
  df_clean <- df %>%
    filter(!is.na(.data[[col_species]]),
           .data[[col_species]] != "",
           !grepl("^__", .data[[col_species]]))
  
  total_genes <- nrow(df_clean)
  na_rows <- sum(is.na(df_clean[[col_human]]) | df_clean[[col_human]] == "")
  aligned <- total_genes - na_rows
  aligned_pct <- if (total_genes > 0) round(100 * aligned / total_genes, 2) else NA_real_
  
  data.frame(
    file = basename(path),
    total_genes = total_genes,
    na_rows = na_rows,
    aligned_rows = aligned,
    aligned_pct = aligned_pct,
    stringsAsFactors = FALSE
  )
}

# Calculate statistics for all files
files <- c(
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_bt_cb_rawcounts_merged.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_ss_cb_rawcounts_merged.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_mf_cb_rawcounts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_cp_cb_raw_counts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_hg_cb_raw_counts.txt",
  "C:/Users/Astrid Xu/Desktop/cerebellum/o_pmpo_cb_raw_counts.txt"
)

stats_all <- do.call(rbind, lapply(files, align_stats_one))
print(stats_all)
out_file <- "C:/Users/Astrid Xu/Desktop/cerebellum/stats_of_all.txt"
write.table(stats_all, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Remove gene columns from merged file
fp <- "C:/Users/Astrid Xu/Desktop/cerebellum/merged_counts_all.txt"
df <- read.table(fp, sep = "\t", header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE, quote = "", comment.char = "")

keep <- seq_along(df) == 1 | !grepl("gene", names(df), ignore.case = TRUE)
df2 <- df[, keep, drop = FALSE]

out <- "C:/Users/Astrid Xu/Desktop/cerebellum/merged_all.txt"
write.table(df2, out, sep = "\t", quote = FALSE, row.names = FALSE)