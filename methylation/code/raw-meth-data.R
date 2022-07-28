get_sample_annotation <- function() {
    cfg <- fread(here("data/meth_sample_annotation.csv"))
    cfg <- cfg %>%
        set_names("seq_id", "sample_id", "sample_description", "eb_day", "mut", "treatment", "illumina.index", "index1", "index1.seq", "index2", "index2.seq") %>%
        mutate(track = glue("EarlyDev{seq_id}.{sample_id}_{sample_description}_{eb_day}_{mut}_{treatment}"))
    cfg <- cfg %>% mutate(eb_day = ifelse(grepl("Serum", sample_description), "S", eb_day))
    cfg <- cfg %>%
        mutate(mut = forcats::fct_recode(mut, "wt" = "J1", "D1" = "Dnmt1KO", "3a" = "Dnmt3aKO", "3b" = "Dnmt3bKO"), eb_day = ifelse(eb_day == "S", eb_day, glue("d{eb_day}"))) %>%
        unite("samp_id", eb_day, mut, remove = FALSE)
    return(cfg)
}

calc_eb_cpg_meth <- function(from, to, min_cov, max_na, intervals, iterator, cache_fn, rm_meth_cov, tracks_key = NULL) {
    tracks_key <- tracks_key %||% get("tracks_key", .GlobalEnv)
    track_df <- tracks_key %>%
        filter(day %in% c("d0S", paste0("d", from:to)), !grepl("N15", line)) %>%
        group_by(day, line) %>%
        mutate(name1 = glue("{day}_{line}_{1:n()}")) %>%
        ungroup() %>%
        select(line, day, name = name1, track_name) %>%
        unite("grp", day, line, remove = FALSE)


    m_all <- gextract_meth(tracks = track_df$track_name, names = track_df$name, intervals = intervals, extract_meth_calls = TRUE, iterator = iterator) %cache_df% cache_fn %>% as_tibble()


    m_merged <- m_all %>% select(chrom, start, end)
    grps <- track_df %>%
        arrange(line, day) %>%
        pull(grp) %>%
        unique()
    for (g in grps) {
        nms <- track_df %>%
            filter(grp == g) %>%
            pull(name)

        cov_col <- paste0(g, ".cov")
        meth_col <- paste0(g, ".meth")

        m_merged[[cov_col]] <- rowSums(m_all[, paste0(nms, ".cov")], na.rm = TRUE)
        m_merged[[meth_col]] <- rowSums(m_all[, paste0(nms, ".meth")], na.rm = TRUE)

        m_merged[[g]] <- ifelse(m_merged[[cov_col]] >= min_cov, m_merged[[meth_col]] / m_merged[[cov_col]], NA)
    }

    if (rm_meth_cov) {
        m_merged <- m_merged %>% select(-ends_with(".meth"), -ends_with(".cov"))
    }

    if (!is.null(max_na)) {
        f_cov <- rowSums(is.na(m_merged[, -1:-3])) <= max_na
    } else {
        f_cov <- rowSums(!is.na(m_merged[, -1:-3])) >= 1
    }

    m_f <- m_merged %>% filter(f_cov)


    colnames(m_f) <- gsub("ko3a", "3a", colnames(m_f))
    colnames(m_f) <- gsub("ko3b", "3b", colnames(m_f))


    return(m_f)
}


calc_eb_day0_to_day4_cpg_meth <- function(min_cov = 10, max_na = 5, intervals = gintervals.union("intervs.captPBAT_probes.ES_EB_V1", "intervs.captPBAT_probes.ES_EB_V2"), iterator = "intervs.global.seq_CG", cache_fn = here("output/eb_day0_to_day4_cpg_meth.tsv"), rm_meth_cov = TRUE) {
    calc_eb_cpg_meth(from = 0, to = 4, min_cov = min_cov, max_na = max_na, intervals = intervals, iterator = iterator, cache_fn = cache_fn, rm_meth_cov = rm_meth_cov, tracks_key = tracks_key %>% filter(line %in% c("wt", "tko", "dko", "ko3a", "ko3b")))
}



calc_eb_day5_to_day6_cpg_meth <- function(min_cov = 10, max_na = 3, intervals = gintervals.union("intervs.captPBAT_probes.ES_EB_V1", "intervs.captPBAT_probes.ES_EB_V2"), iterator = "intervs.global.seq_CG", cache_fn = here("output/eb_day5_to_day6_cpg_meth.tsv"), rm_meth_cov = TRUE) {
    calc_eb_cpg_meth(from = 5, to = 6, min_cov = min_cov, max_na = max_na, intervals = intervals, iterator = iterator, cache_fn = cache_fn, rm_meth_cov = rm_meth_cov, tracks_key = tracks_key %>% filter(line %in% c("wt", "tko", "dko", "ko3a", "ko3b")))
}

# sperates also to germ layers (sorting)
calc_eb_day0_to_day6_cpg_meth <- function(min_cov = 10, max_na = 5, intervals = gintervals.union("intervs.captPBAT_probes.ES_EB_V1", "intervs.captPBAT_probes.ES_EB_V2"), iterator = "intervs.global.seq_CG", cache_fn = here("output/eb_day0_to_day6_cpg_meth.tsv"), use_sort = TRUE) {
    track_df <- tracks_key %>%
        filter(day %in% c("d0S", paste0("d", 0:6)), !grepl("N15", line)) %>%
        group_by(day, line) %>%
        mutate(name1 = glue("{day}_{line}_{sort}_{1:n()}")) %>%
        ungroup() %>%
        select(line, day, sort, name = name1, track_name)

    if (use_sort) {
        track_df <- track_df %>%
            unite("grp", day, line, sort, remove = FALSE)
    } else {
        track_df <- track_df %>%
            unite("grp", day, line, remove = FALSE)
    }

    m_all <- gextract_meth(tracks = track_df$track_name, names = track_df$name, intervals = intervals, iterator = iterator, extract_meth_calls = TRUE) %cache_df% cache_fn %>% as_tibble()

    m_merged <- m_all %>% select(chrom, start, end)
    grps <- track_df %>%
        arrange(line, day) %>%
        pull(grp) %>%
        unique()

    for (g in grps) {
        nms <- track_df %>%
            filter(grp == g) %>%
            pull(name)

        cov_col <- paste0(g, ".cov")
        meth_col <- paste0(g, ".meth")

        if (all(has_name(m_all, c(paste0(nms, ".cov"), paste0(nms, ".meth"))))) {
            m_merged[[cov_col]] <- rowSums(m_all[, paste0(nms, ".cov")], na.rm = TRUE)
            m_merged[[meth_col]] <- rowSums(m_all[, paste0(nms, ".meth")], na.rm = TRUE)

            m_merged[[g]] <- ifelse(m_merged[[cov_col]] >= min_cov, m_merged[[meth_col]] / m_merged[[cov_col]], NA)
        }
    }

    m_merged <- m_merged %>% select(-ends_with("ko1"))

    f_cov <- rowSums(is.na(m_merged[, -1:-3])) <= max_na

    m_f <- m_merged %>% filter(f_cov)

    colnames(m_f) <- gsub("ko3a", "3a", colnames(m_f))
    colnames(m_f) <- gsub("ko3b", "3b", colnames(m_f))

    return(m_f)
}
