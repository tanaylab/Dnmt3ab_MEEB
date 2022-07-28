
library("Matrix")

library(metacell)
scdb_init("scrna_db", force_reinit = T)
scfigs_init("figs/")

bats <- read.table("config/eseb_select_batches_Aug20.txt", h = T, stringsAsFactors = F)
rownames(bats) <- bats$Amp.Batch.ID

key <- read.table("data/sc_amp_batches_es_eb.txt", h = T, sep = "\t", stringsAsFactors = F)
fkey <- key[key$Amp.Batch.ID %in% bats$Amp.Batch.ID, ]
fkey$exp_state_r <- bats[fkey$Amp.Batch.ID, "Type"]
fkey$exp_state <- sub("_r\\d", "", fkey$exp_state)
fkey$rep <- bats[fkey$Amp.Batch.ID, "rep"]
fkey <- fkey[!is.na(fkey$rep), ]
fkey$perturb <- sub("_\\d+[_g]*", "", fkey$exp_state)
is_dko <- fkey$perturb == "DKO"
is_dko16 <- fkey$perturb == "DKO16"
fkey$perturb[is_dko16] <- "DKO"
fkey$line <- ifelse(grepl("N15S", fkey$perturb), "N15S", ifelse(grepl("N15", fkey$perturb), "N15", "J1"))
fkey$line[is_dko] <- "DKO"
fkey$line[is_dko16] <- "DKO16"
fkey$perturb <- sub("J1_", "", fkey$perturb)
fkey$perturb <- sub("N15_", "", fkey$perturb)
fkey$perturb <- sub("N15S_", "", fkey$perturb)
fkey$EB_day <- as.integer(as.character(fkey$EB_day))

import_set <- function(prefix, mat_nm, from_day, to_day) {
    fkey_early <- fkey[as.numeric(fkey$EB_day) <= to_day &
        as.numeric(fkey$EB_day) >= from_day &
        fkey$perturb == prefix, ]
    fkey_early <- fkey_early[grep("serum", fkey_early$exp_state, invert = T), ]
    fkey_early <- fkey_early[grep("TODO", fkey_early$exp_state, invert = T), ]
    fkey_early <- fkey_early[grep("mix", fkey_early$exp_state, invert = T), ]
    fkey_early <- fkey_early[grep("V6", fkey_early$exp_state, invert = T), ]
    fkey_early <- fkey_early[grep("vitA", fkey_early$exp_state, invert = T), ]
    # remove gated experiments?
    fkey_early <- fkey_early[grep("P1", fkey_early$exp_state, invert = T), ]
    key_fn <- sprintf("config/eseb_key_07_%s.txt", prefix)
    write.table(fkey_early, key_fn, quote = F, sep = "\t")

    message("Will import ", mat_nm)

    mcell_import_multi_mars(
        mat_nm = mat_nm,
        dataset_table_fn = key_fn,
        base_dir = "data/raw/",
        patch_cell_name = T,
        force = TRUE
    )
    mat <- scdb_mat(mat_nm)
    message("got ", ncol(mat@mat), " cells")
    nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
    bad_genes <- c(grep("^Mtrnr", nms, v = T), "Neat1", grep("ERCC", nms, v = T), "Atpase6", "Xist", "Malat1", "Tmsb4x")
    mcell_mat_ignore_genes(mat_nm, mat_nm,
        bad_genes,
        reverse = F
    )
    mcell_mat_ignore_small_cells(mat_nm, mat_nm, 1000)
}

import_set("WT", "eseb_07_wt", from_day = 0, to_day = 7)
# import_set("DNMT1", "eseb_07_D1", from_day=0, to_day=7)
import_set("DNMT3A", "eseb_07_3a", from_day = 0, to_day = 7)
import_set("DNMT3B", "eseb_07_3b", from_day = 0, to_day = 7)
# import_set("DKO", "eseb_07_3ab", from_day=0, to_day=7)
# import_set("TKO", "eseb_07_TKO", from_day=0, to_day=7)
