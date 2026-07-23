##TF selection
library(httr)
library(jsonlite)
library(dplyr)
install.packages("readxl")  
library(readxl)
excel_sheets("C:/Users/Astrid Xu/Downloads/1-s2.0-S1874939921000833-mmc1.xlsx")


human_tf <- read_excel("C:/Users/Astrid Xu/Downloads/1-s2.0-S1874939921000833-mmc1.xlsx", sheet = "A. Human dbTF catalogue")
head(human_tf)
TF_symbols<- unique (human_tf$`HGNC approved gene symbol`)
head(TF_symbols)

loading_df <- as.data.frame(pca_res$x)  # gene-wise PCA
loading_df$gene <- rownames(loading_df)

top_pc1_genes <- loading_df %>%
  arrange(desc(PC3)) %>%
  slice(1:200) %>%
  pull(gene)

bottom_pc1_genes <- loading_df %>%
  arrange(PC3) %>%
  slice(1:200) %>%
  pull(gene)

#merge
tb_genes <- unique(c(top_pc1_genes, bottom_pc1_genes))


PC2_tfs <- intersect(tb_genes, TF_symbols)
show(PC2_tfs)
PC2_TF <- data.frame(PC2_TF= PC2_tfs)
write.csv(PC2_TF, "PC2_TF.csv",row.names = FALSE)
getwd()




#GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

install.packages("Matrix")
install.packages("Seurat")

gene_df <- bitr(tb_genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)


missing_genes <- setdiff(tb_genes, gene_df$SYMBOL)
if (length(missing_genes) > 0) {
  message("The following genes could not be converted and were skipped: ", paste(missing_genes, collapse = ", "))
}

#GO
ego_tb <- enrichGO(gene         = gene_df$ENTREZID,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
print(head(ego_tb))


dotplot(ego_tb, showCategory = 15) + ggtitle("GO Enrichment of PC3 Top/Bottom Genes")


fisher_compare <- corr_by_group %>%
  select(PC, TF, group, rho, n) %>%
  pivot_wider(names_from = group,
              values_from = c(rho, n),
              names_sep = "_") %>%
  mutate(
    # Fisher
    z_pri   = atanh(rho_Primate),
    z_other = atanh(rho_Other),
    se      = sqrt(1/(pmax(n_Primate, 4) - 3) + 1/(pmax(n_Other, 4) - 3)),
    z_diff  = (z_pri - z_other) / se,
    p_diff  = 2 * pnorm(-abs(z_diff))
  )

#FDR
fisher_compare <- fisher_compare %>%
  group_by(PC) %>%
  mutate(FDR_diff = p.adjust(p_diff, method = "BH")) %>%
  ungroup()

top_diff <- fisher_compare %>%
  arrange(FDR_diff, desc(abs(z_diff))) %>%
  slice_head(n = 30)
top_diff

write.csv(corr_by_group, "TF_PC_spearman_by_group.csv", row.names = FALSE)
write.csv(fisher_compare, "TF_PC_group_correlation_difference.csv", row.names = FALSE)


expr_long_grp <- expr_mat_cs[, tf_cols, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("organism") %>%
  pivot_longer(-organism, names_to = "TF", values_to = "expr") %>%
  left_join(species_group_df, by = "organism")

expr_diff <- expr_long_grp %>%
  group_by(TF, group) %>%
  summarise(med = median(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = med) %>%
  mutate(log2FC_Primate_vs_Other = log2((Primate + 1e-9)/(Other + 1e-9)))

#Wilcoxon for every TF
expr_wcx <- expr_long_grp %>%
  group_by(TF) %>%
  summarise(
    p_wilcox = tryCatch(wilcox.test(expr ~ group, exact = FALSE)$p.value,
                        error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(FDR_wilcox = p.adjust(p_wilcox, method = "BH"))

expr_diff_full <- expr_diff %>% left_join(expr_wcx, by = "TF")
write.csv(expr_diff_full, "TF_expression_Primate_vs_Other.csv", row.names = FALSE)

##visualization
plot_tf_vs_pc_grouped <- function(tf, pc_name = "PC1_mean") {
  stopifnot(tf %in% colnames(expr_mat_cs))
  spp <- rownames(expr_mat_cs)
  df <- tibble(
    organism = spp,
    TF_expr = expr_mat_cs[, tf],
    PC = pc_mat[, pc_name]
  ) %>%
    left_join(species_group_df, by = "organism")
  
  #Spearman
  ann <- df %>%
    group_by(group) %>%
    summarise(
      rho = {
        out <- safe_spearman(TF_expr, PC)
        as.numeric(out["rho"])
      },
      n = sum(is.finite(TF_expr) & is.finite(PC)),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(group, ": ρ=", round(rho, 3), " (n=", n, ")"))
  
  ggplot(df, aes(x = TF_expr, y = PC, label = organism, color = group, shape = group)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste0(tf, " vs ", pc_name, " by group"),
         x = paste0(tf, " mean expression (species-level)"),
         y = pc_name) +
    theme_classic() +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 1)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = paste(ann$label, collapse = "\n"))
}