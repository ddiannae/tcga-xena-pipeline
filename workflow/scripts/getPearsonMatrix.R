library(data.table)
library(dplyr)
library(magick)
library(ComplexHeatmap)

cat("Reading matrix \n")
expr_file <- "/datos/ot/diana/regulacion-trans-take1/prostate/rdata/gc-full_length-no_tmm_norm_cpm10_arsyn_cancer.tsv"
annot_file <- "/datos/ot/diana/regulacion-trans-take1/prostate/rdata/annot.RData"

cat("Getting annotations \n")
chr_pal <- c("#D909D1", "#0492EE", "#D2656C", "#106F35", "#5BD2AE", "#199F41", 
             "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
             "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#651532", "#B925AE",
             "#86C91F", "#CB97E8", "#130B9E", "#EF782B", "#79A5B9", "#F7AA7C")

names(chr_pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                    "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21", "Y")

load(annot_file)
expr_matrix <- fread(expr_file)
annot <- annot %>% 
  mutate(chr = factor(chr, levels=c(as.character(1:22), "X", "Y"))) %>%
  arrange(chr, start)
genes <- expr_matrix %>% pull(gene)
expr_matrix <- expr_matrix %>% select(-gene) %>% as.matrix()
rownames(expr_matrix) <- genes

annot <- annot %>% filter(gene_id %in% genes)
expr_matrix <- expr_matrix[annot$gene_id, ]

expr_matrix <- t(expr_matrix)
cat("Calculating correlation \n")
corr <- cor(expr_matrix, expr_matrix, method="pearson")

chrs <- annot %>% dplyr::select(chr) %>% dplyr::rename(Chr = chr)

cat("Building heatmap \n")
ha <- rowAnnotation(df = chrs, 
                    name = "Chr", show_annotation_name = F,
                    col = list(Chr = chr_pal[as.character(chrs$Chr)]), 
                    width = unit(0.5, "cm"),
                    annotation_legend_param = list(
                      title = "Chr", title_gp = gpar(fontsize = 14), grid_height = unit(0.5, "cm")))


ht <- Heatmap(corr, cluster_rows = F, cluster_columns = F, show_row_names = F, 
        show_column_names = F,  heatmap_legend_param = list( title_gp = gpar(fontsize = 14),
                                                             title = "Pearson Correlation",
                                                             legend_height = unit(4, "cm")), 
        right_annotation = ha)

cat("Saving image \n")
png("/datos/ot/diana/regulacion-trans-take1/prostate/plots/cancer_heatmap.png", width = 1200, height = 1200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
