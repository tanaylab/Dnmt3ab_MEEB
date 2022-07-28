plot_spatial_meth <- function(intervals, tracks, expand_s = 0, expand_e = 0, min_cov = 10, names = NULL, md = NULL, k = NULL, scope = 5e5, separate_plots = FALSE) {
    expanded_intervals <- intervals %>% mutate(start = start - expand_s, end = end + expand_e)
    d <- gpatterns.get_avg_meth(tracks = tracks, names = names, intervals = expanded_intervals, iterator = NULL, min_cov = 1, min_samples = 1)
    if (!is.null(k)) {
        d <- d %>%
            group_by(samp) %>%
            arrange(samp, start) %>%
            mutate(meth_s = zoo::rollsum(meth, k = k, fill = "extend"), cov_s = zoo::rollsum(cov, k = k, fill = "extend")) %>%
            mutate(avg_s = meth_s / cov_s) %>%
            select(chrom:end, cov = cov_s, samp, avg = avg_s)
    }

    d1 <- d %>%
        filter(cov >= min_cov) %>%
        select(chrom:end, samp, avg)
    if (!is.null(md)) {
        d1 <- d1 %>% left_join(md)
    }
    p_meth <- d1 %>% ggplot(aes(x = start, y = avg, color = samp, group = samp)) +
        geom_point() +
        geom_line() +
        scale_x_continuous(limits = c(expanded_intervals$start, expanded_intervals$end)) +
        geom_vline(xintercept = c(intervals$start, intervals$end), color = "gray", linetype = "dashed") +
        ggsci::scale_color_lancet(name = "") +
        ylab("Avg. Meth.") +
        xlab("")
    p_meth <- p_meth + ggtitle(glue::glue("{intervals$chrom}: {scales::comma(expanded_intervals$start)}-{scales::comma(expanded_intervals$end)}"))

    scope_intervals <- intervals %>% mutate(start = start - scope, end = end + scope)

    p_genome <- expanded_intervals %>%
        # distinct(start) %>%
        ggplot(aes(x = start, xend = start, y = 1, yend = 0)) +
        geom_segment(color = "red") +
        # geom_point(size = 0.5) +
        ylim(0, 5) +
        scale_x_continuous(limits = c(scope_intervals$start, scope_intervals$end)) +
        xlab("") +
        ylab("") +
        ggpubr::theme_pubr() +
        theme(
            text = element_text(family = "ArialMT", size = 6),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )

    genes <- gintervals.neighbors1("intervs.global.tss", scope_intervals) %>%
        filter(dist == 0) %>%
        distinct(geneSymbol, .keep_all = TRUE) %>%
        select(chrom, start, end, strand, name = geneSymbol)

    for (i in 1:nrow(genes)) {
        p_genome <- p_genome +
            annotate("segment",
                x = genes$start[i],
                xend = genes$start[i],
                y = 0,
                yend = 1,
                color = "blue"
            ) +
            annotate("segment",
                x = genes$start[i],
                xend = genes$start[i] + genes$strand[i] * (scope_intervals$end - scope_intervals$start) * 0.03,
                y = 1,
                yend = 1,
                color = "blue",
                arrow = arrow(length = unit(0.05, "inches"))
            ) +
            annotate("text", label = genes$name[i], x = genes$start[i], y = 1.5, size = 2)
    }

    if (separate_plots) {
        return(list(p_genome = p_genome, p_meth = p_meth))
    }

    margin <- margin(t = 0, b = 0, r = 3, l = 3, unit = "pt")
    p <- cowplot::plot_grid(cowplot::plot_grid(p_genome + theme(plot.margin = margin), p_meth + guides(color = "none") + theme(plot.margin = margin), ncol = 1, align = "v", rel_heights = c(0.3, 0.7)), cowplot::get_legend(p_meth), nrow = 1, rel_widths = c(0.8, 0.2))

    return(p)
}

names_to_tracks <- function(names, new_names = NULL) {
    md <- rbind(tracks_key, ext_tracks_key)
    md <- md %>% filter(name %in% names)
    if (!is.null(new_names)) {
        md <- md %>% left_join(tibble(name = names, new_name = new_names))
        md <- md %>% select(track_name, name = new_name)
    }

    md <- md %>% select(track_name, name)
    return(md)
}

get_smooth_meth_data_ebdnmt <- function(names, new_names = names, ...) {
    md <- names_to_tracks(names, new_names)
    get_smooth_meth_data(md$track_name, names = md$name, ...)
}



get_smooth_meth_data <- function(tracks, names = NULL, min_cov = 20, d_expand = 100, hits_iterator = d_expand / 2, min_samples = NULL, screen = TRUE, intervals = NULL) {
    cov_vtracks <- glue("{names}_smoo.cov")
    walk2(cov_vtracks, glue("{tracks}.cov"), ~ gvtrack.create(.x, .y, "sum"))
    walk(cov_vtracks, ~ gvtrack.iterator(.x, sshift = -d_expand, eshift = d_expand))

    meth_vtracks <- glue("{names}_smoo.meth")
    walk2(meth_vtracks, glue("{tracks}.meth"), ~ gvtrack.create(.x, .y, "sum"))
    walk(meth_vtracks, ~ gvtrack.iterator(.x, sshift = -d_expand, eshift = d_expand))

    on.exit(walk(c(meth_vtracks, cov_vtracks), gvtrack.rm))

    if (screen) {
        if (is.null(min_samples)) {
            screen_expr <- paste(glue("{cov_vtracks} >= {min_cov}"), collapse = " & ")
        } else if (min_samples == 1) {
            screen_expr <- paste(glue("{cov_vtracks} >= {min_cov}"), collapse = " | ")
        } else {
            screen_expr <- paste(glue("{cov_vtracks} >= {min_cov}"), collapse = ", ")
            screen_expr <- glue("psum({screen_expr}, na.rm=TRUE) >= {min_samples}")
        }

        hits <- gscreen(screen_expr, iterator = hits_iterator)
        hits_mid <- hits %>% mutate(start = (start + end) / 2, end = start + 1)
    } else {
        hits_mid <- "intervs.global.seq_CG"
    }


    if (is.null(intervals)) {
        intervals <- hits_mid
    }

    hits_data <- gextract(c(glue("{meth_vtracks} / {cov_vtracks}"), cov_vtracks), iterator = hits_mid, intervals = intervals, colnames = c(names, glue("{names}.cov"))) %>%
        arrange(intervalID) %>%
        select(-intervalID) %>%
        as_tibble()

    replace_list <- map(glue("{names}.cov"), ~0) %>% set_names(glue("{names}.cov"))
    hits_data <- hits_data %>% replace_na(replace = replace_list)

    return(hits_data)
}

gextract_meth_ebdnmt <- function(names, new_names = names, ...) {
    md <- names_to_tracks(names, new_names)
    gextract_meth(md$track_name, names = md$name, ...)
}
