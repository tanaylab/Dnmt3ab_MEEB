library(tidyverse)
library(tgppt)
library(tgutil)
library(glue)
library(misha)
library(here)
library(patchwork)
library(misha.ext)
library(tgstat)

set.seed(17)

theme_set(tgppt::theme_arial(8))
gsetroot(here("db"))
options(gmax.data.size = 1e9, gmultitasking = FALSE, tgutil.verbose = FALSE)

source(here("code/utils.R"))

scripts <- list.files(here::here("code"), pattern = "\\.(r|R)$", full.names = TRUE, recursive = TRUE)
scripts <- scripts[!(scripts %in% c(here("code/init.R"), here("code/utils.R")))]

walk(scripts, source)

tracks_key <- fread(here("data/tracks_key.tsv"))
tracks_key <- tracks_key %>%
    group_by(day, line, sort, capture, group) %>%
    mutate(replicate = 1:n()) %>%
    mutate(name = glue("{day}_{line}_{sort}_{capture}_r{replicate}_{group}"), name = gsub("\\.", "_", name), name = gsub("/", "_", name), name = gsub("\\?", "", name)) %>%
    ungroup()

track_stats <- calc_track_stats() %cache_df% here("data/track_stats.csv") %>% filter(CHH <= 0.02)
tracks_key <- tracks_key %>% filter(track_name %in% track_stats$track)

ext_tracks_key <- fread(here("data/ext_tracks_key.tsv")) %>%
    mutate(name = glue("{day}_{line}_{sort}_{capture}_r{replicate}_{group}")) %>%
    mutate(sort = ifelse(sort == "Epi", "epi", sort)) %>%
    mutate(sort_fraction = NA)
chip_tracks_key <- fread(here("data/chip_tracks.tsv"))


meth_md <- tracks_key %>%
    mutate(meth_track = glue("{track_name}.meth"), cov_track = glue("{track_name}.cov"), meth_vtrack = glue("{name}.meth"), cov_vtrack = glue("{name}.cov")) %>%
    mutate(expr = glue("{meth_vtrack} / {cov_vtrack}"))

ext_meth_md <- ext_tracks_key %>%
    mutate(meth_track = glue("{track_name}.meth"), cov_track = glue("{track_name}.cov"), meth_vtrack = glue("{name}.meth"), cov_vtrack = glue("{name}.cov")) %>%
    mutate(expr = glue("{meth_vtrack} / {cov_vtrack}"))

init_vtracks <- function() {
    gvtrack.create("d_exon", "intervs.global.exons", "distance")
    gvtrack.create("d_tss", "intervs.global.tss", "distance")
    gvtrack.create("tor", "Encode.esd3.replichip.rep2", "avg")
    gvtrack.create("ab_score", "DNMT.ab_score", "avg")
    gvtrack.create("a_score", "DNMT.a_score", "avg")
    gvtrack.create("b_score", "DNMT.b_score", "avg")
    gvtrack.iterator("tor", sshift = -15000, eshift = 15000)
    gvtrack.create("cg_cont", "seq.CG_500_mean")
    gvtrack.create("gc_cont", "seq.GC_500_mean")

    walk2(chip_tracks_key$name, chip_tracks_key$track_name, ~ gvtrack.create(.x, .y, "global.percentile"))
    walk2(meth_md$meth_vtrack, meth_md$meth_track, ~ gvtrack.create(.x, .y, "sum"))
    walk2(meth_md$cov_vtrack, meth_md$cov_track, ~ gvtrack.create(.x, .y, "sum"))
}

init_vtracks()

md <- rbind(tracks_key[, intersect(colnames(tracks_key), colnames(ext_tracks_key))], ext_tracks_key[, intersect(colnames(tracks_key), colnames(ext_tracks_key))])

germ_layer_colors <<- c("epi" = "lightblue", "ps" = "orange", "ve" = "pink", "ecto" = "#00468BFF", "epi/ecto" = "#00468BFF", "endo" = "#ED0000FF", "meso" = "#1B1919FF", "ExE" = "#42B540FF", "blood?" = "#925E9FFF", "other" = "gray")
