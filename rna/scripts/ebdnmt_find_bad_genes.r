library(metacell)
library(Matrix)
scdb_init("scrna_db", force_reinit = T)

by_line <- 1
if (0) {
    mat_nm <- "eseb_07_wt"
    rep_id1 <- 4
    rep_id2 <- 5
}

if (1) {
    mat_nm <- "eseb_07_3a"
    rep_id1 <- 3
    rep_id2 <- 4
    rep_id1 <- "J1"
    rep_id2 <- "N15"
}

if (0) {
    mat_nm <- "eseb_07_3b"
    rep_id1 <- 3
    rep_id2 <- 4
}
if (0) {
    mat_nm <- "eseb_07_3ab"
    rep_id1 <- 2
    rep_id2 <- 3
}
if (0) {
    mat_nm <- "eseb_07_TKO"
    rep_id1 <- 1
    rep_id2 <- 2
}
# mat_nm = "eseb_07_D1"

T_vm <- 0.1
T_top3 <- 3
T_tot <- 50
K <- 80

mat <- scdb_mat(mat_nm)
md <- mat@cell_metadata

all_c <- colnames(mat@mat)
if (by_line) {
    rep1 <- intersect(all_c, rownames(md)[as.character(md$line) == rep_id1])
    rep2 <- intersect(all_c, rownames(md)[as.character(md$line) == rep_id2])
} else {
    rep1 <- intersect(all_c, rownames(md)[md$rep == rep_id1])
    rep2 <- intersect(all_c, rownames(md)[md$rep == rep_id2])
}

bat1 <- c()
bat2 <- c()
for (day in 0:7) {
    c_day <- rownames(md)[md$EB_day == day]
    cd1 <- intersect(rep1, c_day)
    cd2 <- intersect(rep2, c_day)

    n <- min(length(cd1), length(cd2))
    bat1 <- c(bat1, sample(cd1, n))
    bat2 <- c(bat2, sample(cd2, n))
}

g_tot1 <- rowSums(mat@mat[, bat1])
g_tot2 <- rowSums(mat@mat[, bat2])
e_tot1 <- g_tot1 / sum(g_tot1)
e_tot2 <- g_tot2 / sum(g_tot2)

bat_fold <- log2(1e-5 + e_tot1) - log2(1e-5 + e_tot2)

mcell_gset_filter_varmean(mat_nm, mat_nm, T_vm = T_vm, force_new = T)
mcell_gset_filter_cov(mat_nm, mat_nm, T_tot = T_tot, T_top3 = T_top3)

mcell_gset_split_by_dsmat(mat_nm, mat_nm, K)

gset <- scdb_gset(mat_nm)

gset_batch <- tapply(bat_fold[names(gset@gene_set)], gset@gene_set, mean)

dir.create(sprintf("figs/%s.gset_cors/", mat_nm), showWarnings = FALSE)
png(sprintf("figs/%s.gset_cors/batch_gset_bars.png", mat_nm), w = 400, h = 800)
barplot(gset_batch, horiz = T)
dev.off()

message("potentially batchy gsets:")
bat_ids <- which(abs(gset_batch) > 0.4)
ids <- bat_ids
for (i in ids) {
    message("Gset ", i)
    print(which(gset@gene_set == i))
}
message("Containing RPs:")
rps <- grep("^Rp[l|s]", names(gset@gene_set), v = T)
rp_ids <- setdiff(unique(gset@gene_set[rps]), ids)
for (i in rp_ids) {
    message("Rp containung gset ", i)
    print(which(gset@gene_set == i))
}
mrps <- grep("^Mrp", names(gset@gene_set), v = T)
hsps <- grep("^Hsp", names(gset@gene_set), v = T)
nduf <- grep("^Nduf", names(gset@gene_set), v = T)
nduf <- grep("^Psmb", names(gset@gene_set), v = T)
stress_ids <- setdiff(unique(gset@gene_set[c(mrps, hsps, nduf)]), ids)
for (i in stress_ids) {
    message("Stress containung gset ", i)
    print(which(gset@gene_set == i))
}
ids <- c(rp_ids, ids)
ids <- c(rp_ids, ids)
cc_anchors <- c("Top2a", "Ube2c", "Pcna", "Mki67", "Cenpf", "Hist1h1a", "Hist1h1b")
cc_ids <- setdiff(unique(gset@gene_set[cc_anchors]), ids)
for (i in cc_ids) {
    message("CC containung gset ", i)
    print(which(gset@gene_set == i))
}


# EDIT THIS
if (0) {
    message("Make sure you editted manually the list of gset that are filtered!!!")
    bad_genes_base <- c(as.character(read.table("data/filtered_genes_eseb07_wt.txt", sep = "\t")$x))
    gset_old_bad <- table(gset@gene_set, names(gset@gene_set) %in% bad_genes_base)
    bad_gsets_wt <- c(52, 43, 33, 68, 72, 35, 66, 79, 69, 80, 3, 9)
    bad_gene_wt <- names(gset@gene_set)[which(gset@gene_set %in% bad_gsets_wt)]
    write.table(bad_genes_wt, "data/filtered_genes_eseb07_wt.txt", quote = F, sep = "\t")
}

if (0) {
    message("Make sure you editted manually the list of gset that are filtered!!!")
    bad_genes_base <- c(as.character(read.table("data/filtered_genes_eseb07_3a.txt", sep = "\t")$x))
    gset_old_bad <- table(gset@gene_set, names(gset@gene_set) %in% bad_genes_base)
    bad_gsets_3a <- c(56, 58, 75, 60, 77, 51, 76, 78, 65)
    bad_genes_3a <- names(gset@gene_set)[which(gset@gene_set %in% bad_gsets_3a)]
    write.table(bad_genes_3a, "data/filtered_genes_eseb07_3a.txt", quote = F, sep = "\t")
}
if (0) {
    message("Make sure you editted manually the list of gset that are filtered!!!")
    bad_genes_base <- c(as.character(read.table("data/filtered_genes_eseb07_3b.txt", sep = "\t")$x))
    gset_old_bad <- table(gset@gene_set, names(gset@gene_set) %in% bad_genes_base)
    bad_gsets_3b <- c(46, 50, 72, 78, 30, 75, 70, 65, 58, 76)
    bad_genes_3b <- names(gset@gene_set)[which(gset@gene_set %in% bad_gsets_3b)]
    write.table(bad_genes_3b, "data/filtered_genes_eseb07_3b.txt", quote = F, sep = "\t")
}
if (0) {
    message("Make sure you editted manually the list of gset that are filtered!!!")
    bad_genes_base <- c(as.character(read.table("data/filtered_genes_eseb07_3ab.txt", sep = "\t")$x))
    gset_old_bad <- table(gset@gene_set, names(gset@gene_set) %in% bad_genes_base)
    bad_gsets_3ab <- c(56, 61, 77, 59, 72, 74, 68, 3)
    bad_genes_3ab <- names(gset@gene_set)[which(gset@gene_set %in% bad_gsets_3ab)]
    write.table(bad_genes_3ab, "data/filtered_genes_eseb07_3ab.txt", quote = F, sep = "\t")
}
if (0) {
    message("Make sure you editted manually the list of gset that are filtered!!!")
    bad_genes_base <- c(as.character(read.table("data/filtered_genes_eseb07_tko.txt", sep = "\t")$x))
    gset_old_bad <- table(gset@gene_set, names(gset@gene_set) %in% bad_genes_base)
    bad_gsets_tko <- c(44, 31, 52, 76, 75, 62, 71, 32, 33, 64, 56)
    bad_genes_tko <- names(gset@gene_set)[which(gset@gene_set %in% bad_gsets_tko)]
    write.table(bad_genes_tko, "data/filtered_genes_eseb07_tko.txt", quote = F, sep = "\t")
}
