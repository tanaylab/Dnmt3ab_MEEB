
# calc_track_on_target_stats <- function(track, regions){
#     print(track)
#     stats <- gpatterns::gpatterns.apply_tidy_cpgs(track, function(x){
#         poss <- x  %>% mutate(end = as.integer(ifelse(end == '-', start + 1, end)), new_start = pmin(start, end), new_end=pmax(start, end), new_end = ifelse(new_end == new_start, new_start + 1, new_end)) %>% select(chrom, start=new_start, end=new_end, read_id) %>% gpatterns::.gpatterns.force_chromosomes() %>% distinct(chrom, start, end)
#         tibble(n = nrow(poss), n_ontar = poss %>% gpatterns::gintervals.filter(regions, max_distance=200) %>% nrow())
#     })

#     stats <- stats %>% summarise(n = sum(n), n_ontar = sum(n_ontar)) %>% mutate(frac_ontar = n_ontar / n)
#     return(stats)
# }

calc_track_on_target_stats <- function(track, regions) {
    print(track)
    stats <- tibble(n = gsummary(paste0(track, ".cov"), intervals = gintervals.all())[5], n_ontar = gsummary(paste0(track, ".cov"), intervals = regions)[5])

    return(stats)
}

calc_on_target_stats <- function() {
    V1_stats <- map_dfr(tracks_key %>% filter(capture == "V1") %>% pull(track_name), function(track) calc_track_on_target_stats(track, regions = "intervs.captPBAT_probes.ES_EB_V1") %>% mutate(track = track))
    V2_stats <- map_dfr(tracks_key %>% filter(capture == "V2") %>% pull(track_name), function(track) calc_track_on_target_stats(track, regions = "intervs.captPBAT_probes.ES_EB_V2") %>% mutate(track = track))
    wg_stats <- map_dfr(tracks_key %>% filter(capture == "wgbs") %>% pull(track_name), function(track) calc_track_on_target_stats(track, regions = gintervals.union("intervs.captPBAT_probes.ES_EB_V1", "intervs.captPBAT_probes.ES_EB_V1")) %>% mutate(track = track))

    stats <- bind_rows(V1_stats, V2_stats, wg_stats) %>%
        organize_samp_id() %>%
        as_tibble()
    return(stats)
}

organize_samp_id <- function(df) {
    df <- df %>%
        left_join(tracks_key %>% select(track = track_name, mut = line, day) %>% unite("samp_id", day, mut, remove = TRUE)) %>%
        mutate(line = ifelse(grepl("N15", samp_id), "N15", "J1")) %>%
        mutate(
            samp_id = gsub(".cov$", "", samp_id),
            samp_id = gsub("_tko", "_TKO", samp_id),
            samp_id = gsub("_N15", "", samp_id),
            samp_id = gsub("_ko3a", "_3a", samp_id),
            samp_id = gsub("_ko3b", "_3b", samp_id),
            samp_id = gsub("d0S_", "S_", samp_id)
        ) %>%
        separate(samp_id, c("day", "mut"), sep = "_", remove = FALSE) %>%
        mutate(samp_id = gsub("e7.5_mouse", "e7.5", samp_id)) %>%
        mutate(samp_id = ifelse(line == "N15", paste0(line, "_", samp_id), samp_id)) %>%
        as_tibble()
}

organize_sc_samp_id <- function(df) {
    df <- df %>%
        unite("samp_id", day, line, remove = FALSE) %>%
        mutate(
            samp_id = gsub("_ko3a", "_3a", samp_id),
            samp_id = gsub("_ko3b", "_3b", samp_id),
            samp_id = gsub("_mouse", "", samp_id)
        )
    df
}


calc_track_stats <- function(tracks = fread(here("data/tracks_key.tsv"))$track_name) {
    tracks <- grep("sc_bulk", tracks, invert = TRUE, value = TRUE)
    stats <- map_dfr(tracks, gpatterns::gpatterns.track_stats)
    stats <- stats %>%
        organize_samp_id() %>%
        as_tibble()
    return(stats)
}


add_total_cov <- function(fn = here("data/tracks_key.tsv")) {
    tracks_key <- fread(fn)

    get_track_tot_cov <- function(x) {
        if (gtrack.exists(glue("{x}.cov"))) {
            return(gsummary(glue("{x}.cov"))[5])
        } else {
            return(NA)
        }
    }

    total_covs <- plyr::laply(tracks_key$track_name, get_track_tot_cov, .parallel = TRUE)
    tracks_key$total <- total_covs
    n_cpgs <- nrow(gintervals.load("intervs.global.seq_CG"))
    tracks_key <- tracks_key %>% mutate(x_cov = total / n_cpgs)
    fwrite(tracks_key, fn, sep = "\t")
    invisible(tracks_key)
}

add_enh_cov <- function(fn = here("data/tracks_key.tsv")) {
    tracks_key <- fread(fn)
    enh_intervs_gl <- define_enhancer_intervs_germ_layers() %cache_df% here("output/enh_intervs_germ_layers.csv") %>% as_tibble()
    enh_intervs_encode <- define_enhancer_intervs() %cache_df% here("output/enh_intervs_encode.csv") %>% as_tibble()
    enh_intervs <- bind_rows(enh_intervs_gl %>% select(chrom, start, end), enh_intervs_encode %>% select(chrom, start, end)) %>% distinct(chrom, start, end)

    get_track_enh_cov <- function(x) {
        if (gtrack.exists(glue("{x}.cov"))) {
            gvtrack.create(vtrack = glue("vt_{x}.cov"), src = glue("{x}.cov"), func = "sum")
            res <- gextract(glue("vt_{x}.cov"), intervals = enh_intervs, iterator = enh_intervs, colnames = "cov") %>%
                replace_na(replace = list(cov = 0)) %>%
                summarise(tot_cov = sum(cov), mean_cov = mean(cov, na.rm = TRUE))
        } else {
            res <- tibble(tot_cov = NA, mean_cov = NA)
        }
        return(res)
    }

    enh_covs <- plyr::ldply(tracks_key$track_name, get_track_enh_cov, .parallel = TRUE)
    tracks_key$enh_total <- enh_covs$tot_cov
    tracks_key$enh_mean <- enh_covs$mean_cov
    fwrite(tracks_key, fn, sep = "\t")
    invisible(tracks_key)
}
