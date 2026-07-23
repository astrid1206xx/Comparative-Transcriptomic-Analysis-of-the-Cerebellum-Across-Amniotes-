
infile  <- "C:/Users/Astrid Xu/Desktop/cerebellum/4 primates_8brain regions_rawcounts.txt"
outfile <- "C:/Users/Astrid Xu/Desktop/cerebellum/primates_CB_rawcounts.txt"
df <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
head(df)
names(df)


cb_idx <- grep("CB", names(df))  

show(cb_idx)
keep_idx <- unique(c(1, cb_idx))
dt_cb <- df[, ..keep_idx]

df <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
head(df)
names(df)


df_cb <- df[, c(1, 20:31)]
colnames(df_cb)[1] <- "Gene.ID"
write.table(df_cb,
            "C:/Users/Astrid Xu/Desktop/cerebellum/pm_cb_rawcounts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cb_idx <- grep("CB", names(df), ignore.case = TRUE)
print(cb_idx)
print(names(df)[cb_idx])

keep_idx <- c(1, cb_idx)
df_cb <- df[, keep_idx, drop = FALSE]

colnames(df_cb)[1] <- "Gene.ID"

write.table(df_cb,
            "C:/Users/Astrid Xu/Desktop/cerebellum/primates_CB_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(dt_cb, outfile, sep = "\t", quote = FALSE)
"C:\Users\Astrid Xu\Desktop\cerebellum\4 primates_8brain regions_rawcounts.txt"
"C:\Users\Astrid Xu\Desktop\cerebellum\le_cb_rawcounts.txt"
"C:\Users\Astrid Xu\Desktop\cerebellum\ss_cb_rawcounts.txt"
"C:\Users\Astrid Xu\Desktop\cerebellum\bt_cb_rawcounts 1.txt"
"C:\Users\Astrid Xu\Desktop\cerebellum\ss 2.txt"
"C:\Users\Astrid Xu\Desktop\cerebellum\pm_all_rawcounts.tsv"


"C:\Users\Astrid Xu\Desktop\cerebellum\cp_cb_raw_counts.txt"

if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(biomaRt); library(dplyr)

infile  <- "C:/Users/Astrid Xu/Desktop/cerebellum/mf_CB_rawcounts.txt"     
outfile <- "C:/Users/Astrid Xu/Desktop/cerebellum/o_mf_cb_rawcounts.txt"

df <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
head(df)

mart_mf <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")  
orth_map <- getBM(
  attributes = c("ensembl_gene_id", 
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type"),
  filters = "ensembl_gene_id",
  values  = unique(df[[1]]),
  mart    = mart_mf
)

#priotize one2one ortholog
orth_best <- orth_map %>%
  filter(!is.na(hsapiens_homolog_ensembl_gene), hsapiens_homolog_ensembl_gene != "") %>%
  mutate(prio = ifelse(hsapiens_homolog_orthology_type=="ortholog_one2one",1,2)) %>%
  arrange(ensembl_gene_id, prio) %>%
  distinct(ensembl_gene_id, .keep_all=TRUE) %>%
  select(ensembl_gene_id, human_ensembl=hsapiens_homolog_ensembl_gene)

library(dplyr)

key <- names(df)[1]  
df_out <- df %>%
  left_join(orth_best, by = setNames("ensembl_gene_id", key)) %>%  
  relocate(human_ensembl, .before = all_of(key))                   

write.table(df_out, outfile, sep="\t", quote=FALSE, row.names=FALSE)

cat("Done. Wrote:", outfile, "\n",
    "Mapped:", sum(!is.na(df_out$human_ensembl)), "/", nrow(df_out), "genes.\n")
