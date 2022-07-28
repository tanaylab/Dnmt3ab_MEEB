library(metacell)

scdb_init("scrna_db", force_reinit = TRUE)

scfigs_init("figs/")
tgconfig::override_params("config/eseb.yaml", package = "metacell")

source("scripts/generic_mc.r")

if (1) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_wt.txt", sep = "\t")$x))
    bad_genes07 <- union(bad_genes07, c(as.character(read.table("data/filter_gene_line_batch_wt.txt", h = T)$x)))
    generate_mc("eseb_07_wt",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        # 			Knn = 100, Knn_core=18,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (0) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_3a.txt", sep = "\t")$x))
    bad_genes07 <- union(bad_genes07, c(as.character(read.table("data/filter_gene_line_batch_3a.txt", h = T)$x)))
    generate_mc("eseb_07_3a",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (1) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_3b.txt", sep = "\t")$x))
    bad_genes07 <- union(bad_genes07, c(as.character(read.table("data/filter_gene_line_batch_3b.txt", h = T)$x)))
    generate_mc("eseb_07_3b",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (0) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_3ab.txt", sep = "\t")$x))
    generate_mc("eseb_07_3ab",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (0) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_D1.txt", sep = "\t")$x))
    generate_mc("eseb_07_D1",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (0) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_tko.txt", sep = "\t")$x))
    generate_mc("eseb_07_TKO",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
if (0) {
    bad_genes07 <- c(as.character(read.table("data/filtered_genes_eseb07_wt.txt", sep = "\t")$x))
    generate_mc("eseb_sEB8_wt",
        color_key = "config/wt05_mc_colorize.txt",
        recompute = T,
        Knn = 80, Knn_core = 15,
        min_mc_sz = 15,
        add_bad_genes = bad_genes07
    )
}
