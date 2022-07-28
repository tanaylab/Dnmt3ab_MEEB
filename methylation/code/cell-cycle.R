calc_cell_cycle_ord <- function(df) {
    mat <- as.matrix(df %>% select(cell_id, early_late_diff, early_late_cov) %>% na.omit() %>% column_to_rownames("cell_id"))
    mat_norm <- t(t(mat) / apply(mat, 2, sd))
    circ <- princurve:::start_circle(mat_norm)
    pc <- princurve::principal_curve(x = mat_norm, start = circ, stretch = 2, smoother = "periodic_lowess")
    cell_ord <- tibble(cell_id = rownames(mat)[pc$ord]) %>% mutate(ord = 1:n())

    df <- df %>% left_join(cell_ord)

    return(list(ord = cell_ord, pc = pc, mat_norm = mat_norm, df = df))
}

wrap_ord <- function(df, n_wrap = 1) {
    orig_df <- df
    for (i in 1:n_wrap) {
        df <- bind_rows(df, orig_df %>% mutate(ord2 = ord2 + i))
    }
    return(df)
}

plot_principal_curve <- function(fit, x) {
    curve <-
        as_data_frame(fit$s) %>%
        mutate(lambda = fit$lambda) %>%
        slice(fit$ord)


    ggplot() +
        geom_point(aes(x = early_late_diff, y = early_late_cov), as.data.frame(x), colour = "darkgray") +
        geom_path(aes(x = early_late_diff, y = early_late_cov), curve) +
        theme_bw()
}

plot_cc_circle <- function(df, point_size = 0.5, fit = NULL) {
    p <- ggplot(df, aes(x = avg_late - avg_early, y = early_late_cov, fill = ord2)) +
        geom_point(size = point_size, stroke = 0.1, shape = 21) +
        scale_fill_gradientn(colors = c("white", "blue", "red", "yellow")) +
        # guides(fill=FALSE) +
        theme(aspect.ratio = 1) +
        xlab("Late meth. - Early meth.") +
        ylab("log2(Early cov. / Late cov.)")

    # if (!is.null(fit)){
    #     curve <-
    #         as_data_frame(fit$s) %>%
    #         mutate(lambda = fit$lambda) %>%
    #         slice(fit$ord)
    #     p <- p +
    #         geom_path(inherit.aes=FALSE, aes(x=early_late_diff, y=early_late_cov), data=curve)
    # }

    return(p)
}

get_cc_segments <- function(df, n_breaks = 2, psi = c(0.2, 0.5), labels = c("S-start", "S-mid", "S-end"), ...) {
    df_wrap <- df %>% wrap_ord(n_wrap = 2)

    cov_model <- loess(early_late_cov ~ ord2, data = df_wrap, span = 0.2)
    trend_df_cov <- df_wrap %>%
        mutate(cov_trend = predict(cov_model)) %>%
        filter(ord2 >= 1, ord2 <= 2) %>%
        mutate(ord2 = ord2 - 1)


    cov_lm_model <- lm(cov_trend ~ ord2, data = trend_df_cov)
    seg_model <- segmented::segmented(cov_lm_model, seg.Z = ~ord2, npsi = n_breaks, psi = psi, ...)

    # plot(seg_model)
    # points(trend_df_cov$ord2, trend_df_cov$cov_trend, col="red", cex=0.1)
    # points(df$ord2, df$early_late_cov, col="blue", cex=0.1)
    # title(df$type[1])

    psi_breaks <- seg_model$psi[, 2]

    df <- df %>%
        mutate(ord_grp = cut(ord2, breaks = c(0, psi_breaks, 1), labels = labels))

    for (i in 2:length(labels)) {
        df[[labels[i]]] <- psi_breaks[i - 1]
    }

    p_cov <- df %>%
        ggplot(aes(x = ord2, y = early_late_cov, color = ord_grp, group = ord_grp)) +
        geom_point() +
        geom_line(aes(x = ord2, y = cov_trend), data = trend_df_cov, inherit.aes = FALSE, color = "red") +
        geom_vline(xintercept = psi_breaks, color = "gray", linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        guides(color = FALSE) +
        ggtitle(df$type[1]) +
        theme_bw()

    print(p_cov)



    return(df)
}

get_cc_stats <- function(df, span = 0.2) {
    df_wrap <- df %>% wrap_ord(n_wrap = 2)

    df_long <- df_wrap %>%
        pivot_longer(cols = c(avg_early, avg_late), names_to = "el_type", values_to = "meth") %>%
        mutate(el_type = gsub("^avg_", "", el_type)) %>%
        arrange(el_type, ord2)

    early_df <- df_long %>% filter(el_type == "early")
    late_df <- df_long %>% filter(el_type == "late")

    early_trend <- loess(meth ~ ord2, data = early_df, span = 0.2) %>% predict()
    late_trend <- loess(meth ~ ord2, data = late_df, span = 0.2) %>% predict()

    trend_df <-
        bind_rows(
            early_df %>% select(ord2) %>% mutate(trend_meth = early_trend) %>% mutate(el_type = "early") %>% filter(ord2 >= 1, ord2 <= 2),
            late_df %>% select(ord2) %>% mutate(trend_meth = late_trend) %>% mutate(el_type = "late") %>% filter(ord2 >= 1, ord2 <= 2)
        )

    stats <- trend_df %>%
        group_by(el_type) %>%
        summarise(
            peak = max(trend_meth),
            peak_pos = ord2[trend_meth == peak] - 1,
            trough = min(trend_meth),
            trough_pos = ord2[trend_meth == trough] - 1,
            amplitude = peak - trough
        )

    # p_early <- early_df %>% filter(ord2 <= 1) %>% ggplot(aes(x=ord2, y=meth)) + geom_point() + geom_vline(xintercept = c(stats$peak_pos[1], stats$trough_pos[1])) + geom_line(aes(x=ord2 - 1, y=trend_meth), color="red", data = trend_df %>% filter(el_type == "early")) + ggtitle(df$type[1])
    # p_late <- late_df %>% filter(ord2 <= 1) %>% ggplot(aes(x=ord2, y=meth)) + geom_point() + geom_vline(xintercept = c(stats$peak_pos[2], stats$trough_pos[2])) + geom_line(aes(x=ord2 - 1, y=trend_meth), color="red", data = trend_df %>% filter(el_type == "late"))
    # print(p_early / p_late)


    return(stats)
}

get_cc_early_late_meth_trend <- function(df, span = 0.2) {
    df <- df %>% wrap_ord(n_wrap = 2)

    model <- loess(avg ~ ord2, data = df, span = span)

    broom::augment(model, df) %>% filter(ord2 >= 1, ord2 <= 2)
}

plot_cc_early_late_meth <- function(df, point_size = 0.5, add_trend_lines = FALSE, phases = NULL, plot_phase_lines = FALSE, y_lim = c(0.65, 1), rm_legend = TRUE, trend_ylim = NULL) {
    alpha <- 1
    if (add_trend_lines) {
        df <- df %>% wrap_ord(n_wrap = 2)
        alpha <- 0.3
    }

    p <- df %>%
        pivot_longer(cols = c(avg_early, avg_late), names_to = "el_type", values_to = "meth") %>%
        mutate(el_type = gsub("^avg_", "", el_type)) %>%
        ggplot(aes(x = ord2, y = meth, color = el_type))


    if (plot_phase_lines) {
        if (is.null(phases)) {
            phases <- c(df[["S-mid"]][1], df[["S-end"]][1]) + 1
        }

        p <- p +
            # geom_vline(xintercept = phases[c(-1, -length(phases))], linetype = "dashed", color="darkgray")
            geom_vline(xintercept = phases, linetype = "dashed", color = "darkgray")
    }


    p <- p +
        geom_point(size = point_size, alpha = alpha) +
        scale_color_manual(name = "", values = c("early" = "cyan", "late" = "darkblue")) +
        xlab("Cell cycle phase") +
        ylab("Meth.") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        ylim(y_lim)



    if (add_trend_lines) {
        p <- p + geom_smooth(method = "loess", span = 0.2, se = FALSE, size = 0.5) + coord_cartesian(xlim = c(1, 2), expand = 0)
        p_early_late_cov <- df %>%
            filter(!is.na(ord2)) %>%
            ggplot(aes(x = ord2, y = early_late_cov)) +
            geom_point(size = point_size, alpha = 1) +
            geom_smooth(method = "loess", span = 0.2, se = FALSE, size = 0.4, color = "red") +
            coord_cartesian(xlim = c(1, 2), expand = 0) +
            geom_vline(xintercept = phases, linetype = "dashed", color = "darkgray") #+
        # theme(axis.text.y = element_text(size = 5))
        if (!is.null(trend_ylim)) {
            p_early_late_cov <- p_early_late_cov + scale_y_continuous(limits = trend_ylim, breaks = trend_ylim)
        }
    } else {
        p_early_late_cov <- df %>%
            filter(!is.na(ord2)) %>%
            ggplot(aes(x = ord2, y = early_late_cov)) +
            geom_point(size = point_size, alpha = alpha)
    }

    if (rm_legend) {
        p <- p + guides(color = FALSE)
    }

    p_early_late_cov <- p_early_late_cov +
        # ylab("Early/Late cov.") +
        theme(
            axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(0, 0, 0, 0)
        )


    p_all <- wrap_plots(p_early_late_cov, p, ncol = 1, heights = c(0.3, 0.7))

    return(p_all)
}

calc_cc_early_late_ab_diff <- function(df, low = "(-1.24,-0.709]", high = "(0.267,1]") {
    df <- df %>%
        mutate(ab_score = forcats::fct_recode(ab_score, "l_ab" = low, "h_ab" = high)) %>%
        filter(ab_score %in% c("h_ab", "l_ab")) %>%
        select(-cov_late, -cov_early, -meth_late, -meth_early) %>%
        pivot_wider(names_from = c(ab_score), values_from = c(avg_late, avg_early)) %>%
        mutate(d_ab_early = avg_early_h_ab - avg_early_l_ab, d_ab_late = avg_late_h_ab - avg_late_l_ab) %>%
        # select(-starts_with("avg_")) %>%
        pivot_longer(cols = starts_with("d_ab"), names_to = "el_type", values_to = "d_ab") %>%
        mutate(el_type = gsub("^d_ab_", "", el_type))



    df <- df %>% wrap_ord(n_wrap = 2)

    return(df)
}


plot_cc_early_late_ab_diff <- function(df, ret_df = FALSE) {
    df <- calc_cc_early_late_ab_diff(df)

    p <- df %>%
        ggplot(aes(x = ord2, y = d_ab, color = el_type)) +
        scale_color_manual(name = "", values = c("early" = "cyan", "late" = "darkblue")) +
        geom_smooth(method = "loess", span = 0.2, se = FALSE, size = 0.5) +
        coord_cartesian(xlim = c(1, 2), expand = 0) +
        ylab("B-phil - A-phil") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )

    if (ret_df) {
        return(df)
    }

    return(p)


    # p <- df %>%
    #     ggplot(aes(x=ord2, y=d_ab, color=el_type)) +
    #         geom_point(size = point_size) +
    #         scale_color_manual(name="", values=c("early" = "cyan", "late" = "darkblue")) +
    #         xlab("Cell cycle phase") +
    #         ylab("A-phil - B-phil")

    # p1 <- cowplot::ggdraw(p + guides(color=FALSE)) +
    #     cowplot::draw_plot(
    #         ggplot( df %>% mutate(ord_grp = cut(ord2, seq(0,1,l=10), include.lowest=TRUE))) +
    #             geom_boxplot(aes(x=ord_grp, y=d_ab, fill=el_type), outlier.shape=NA, width=0.3, position=position_dodge(width=0.5)) +
    #             scale_fill_manual(name="", values=c("early" = "cyan", "late" = "darkblue")) +
    #             guides(fill=FALSE) +
    #             xlab("") +
    #             ylab("") +
    #             theme(plot.background = element_rect(fill = "transparent",colour = NA)))

    # df <- df %>%
    #     mutate(ord_grp = cut(ord2, breaks=phases, labels=phase_names, include.lowest=TRUE))

    # p2 <- df %>%
    #     # mutate(ord_grp = cut(ord2, seq(0,1,l=n_bins), include.lowest=TRUE)) %>%
    #     ggplot(aes(x=ord_grp, y=d_ab, fill=el_type, color=el_type)) +
    #         ggforce::geom_sina(size = point_size, alpha=0.2, position=position_dodge(width=1)) +
    #         geom_boxplot(width=0.3, position=position_dodge(width=1), color="black", outlier.shape=NA, fatten = 0.5, lwd = 0.5) +
    #         # geom_boxplot(color="black", outlier.shape=NA) +
    #         scale_fill_manual(name="", values=c("early" = "cyan", "late" = "darkblue")) +
    #         scale_color_manual(name="", values=c("early" = "cyan", "late" = "darkblue")) +
    #         # theme(axis.text.x=element_blank(),
    #         #     axis.ticks.x=element_blank()) +
    #         xlab("Cell cycle phase") +
    #         ylab("B-phil - A-phil")


    # if (ret_df){
    #     return(df)
    # }

    # return(p2)
}


plot_cell_cycle_ord_diagnostics <- function(l, title = NULL) {
    cell_ord <- l$ord
    df <- l$df
    p_curve <- plot_principal_curve(l$pc, l$mat_norm)
    if (!is.null(title)) {
        p_curve <- p_curve + ggtitle(title)
    }
    print(p_curve)
    p_early_late_cov <- df %>%
        filter(!is.na(ord)) %>%
        ggplot(aes(x = ord, y = early_late_cov)) +
        geom_point()
    p_early_late_avg <- df %>%
        select(cell_id, early = avg_early, late = avg_late, ord) %>%
        gather("type", "meth", -cell_id, -ord) %>%
        ggplot(aes(x = ord, y = meth, color = type)) +
        geom_point(size = 0.5) +
        ggsci::scale_color_lancet(name = "")
    print(p_early_late_cov / p_early_late_avg)
}


get_S_cpgs <- function(segmented = fread(here("output/cell_cycle/segmented.tsv")), early = TRUE, min_cov = 10, min_cov_bulk = 20, min_meth_bulk = 0.7) {
    if (early) {
        segmented_wt <- segmented %>%
            group_by(type) %>%
            mutate(phase = cut(avg_early, breaks = quantile(avg_early, c(0, 0.3, 0.7, 1)), labels = c("S", "other", "G1"), include.lowest = TRUE, na.rm = TRUE)) %>%
            mutate(phase = ifelse(ord2 <= 0.75 & ord2 >= 0.25, as.character(phase), "other")) %>%
            ungroup()
    } else {
        segmented_wt <- segmented %>%
            group_by(type) %>%
            mutate(phase = cut(avg_late, breaks = quantile(avg_late, c(0, 0.3, 0.7, 1)), labels = c("S", "other", "G1"), include.lowest = TRUE, na.rm = TRUE)) %>%
            ungroup()
    }


    message("number of cells:")
    print(segmented_wt %>% count(phase))

    message("extracting sc methylation")
    s_cpgs_raw <- db_f %>%
        left_join_cells(segmented_wt %>% select(cell_id, phase) %>% filter(phase != "other")) %>%
        filter_cells(!is.na(phase)) %>%
        group_by_cells(phase) %>%
        summarise()

    if (early) {
        s_cpgs <- s_cpgs_raw %>%
            filter(tor > 0, cov >= min_cov) %>%
            mutate(avg = meth / cov) %>%
            pivot_wider(names_from = "phase", values_from = c("meth", "cov", "avg"))
    } else {
        s_cpgs <- s_cpgs_raw %>%
            filter(tor < 0, cov >= min_cov) %>%
            mutate(avg = meth / cov) %>%
            pivot_wider(names_from = "phase", values_from = c("meth", "cov", "avg"))
    }


    print(nrow(s_cpgs))

    message("extracting bulk methylation")
    bulk_meth <- gextract_meth(c("Zhang_Nature_Genetics_2017.Ect_mCG", "Zhang_Nature_Genetics_2017.Mes_mCG", "Zhang_Nature_Genetics_2017.End_mCG"), names = c("ecto", "meso", "endo"), intervals = s_cpgs %>% select(chrom, start, end), iterator = s_cpgs %>% select(chrom, start, end), min_cov = min_cov_bulk)


    s_cpgs_f <- s_cpgs %>%
        left_join(bulk_meth %>% select(chrom:end, ecto, meso, endo)) %>%
        filter(ecto >= min_meth_bulk)

    message("downsampling")
    dsn <- min_cov
    s_cpgs_dsn <- s_cpgs_f %>%
        filter(!is.na(meth_S)) %>%
        filter(cov_S >= dsn) %>%
        select(chrom, start, end, meth = meth_S) %>%
        mutate(m = 1) %>%
        uncount(meth) %>%
        bind_rows(s_cpgs_f %>%
            filter(!is.na(meth_S)) %>%
            filter(cov_S >= dsn) %>%
            mutate(unmeth = cov_S - meth_S) %>%
            select(chrom, start, end, unmeth) %>%
            mutate(m = 0) %>% uncount(unmeth)) %>%
        group_by(chrom, start, end) %>%
        sample_n(dsn) %>%
        summarise(meth = sum(m), unmeth = dsn - meth, cov = dsn, .groups = "drop") %>%
        left_join(s_cpgs_f %>% select(chrom:b_score, avg_S, avg_G1, ecto))

    s_cpgs_dsn <- s_cpgs_dsn %>% mutate(avg = meth / cov)

    message(glue("# of CpGs: {nrow(s_cpgs_dsn)}"))
    return(list(s_cpgs = s_cpgs_dsn, segmented = segmented_wt))
}


get_G1_cpgs <- function(segmented = fread(here("output/cell_cycle/segmented.tsv")), early = TRUE, min_cov = 10, min_cov_bulk = 20, min_meth_bulk = 0.7) {
    if (early) {
        segmented_wt <- segmented %>%
            group_by(type) %>%
            mutate(phase = cut(avg_early, breaks = quantile(avg_early, c(0, 0.3, 0.7, 1)), labels = c("S", "other", "G1"), include.lowest = TRUE, na.rm = TRUE)) %>%
            mutate(phase = ifelse(phase == "S", ifelse(ord2 <= 0.75 & ord2 >= 0.25, as.character(phase), "other"), as.character(phase))) %>%
            ungroup()
    } else {
        segmented_wt <- segmented %>%
            group_by(type) %>%
            mutate(phase = cut(avg_late, breaks = quantile(avg_late, c(0, 0.3, 0.7, 1)), labels = c("S", "other", "G1"), include.lowest = TRUE, na.rm = TRUE)) %>%
            ungroup()
    }


    message("number of cells:")
    print(segmented_wt %>% count(phase))

    message("extracting sc methylation")
    s_cpgs_raw <- db_f %>%
        left_join_cells(segmented_wt %>% select(cell_id, phase) %>% filter(phase != "other")) %>%
        filter_cells(!is.na(phase)) %>%
        group_by_cells(phase) %>%
        summarise()


    if (early) {
        s_cpgs <- s_cpgs_raw %>%
            filter(tor > 0, cov >= min_cov) %>%
            mutate(avg = meth / cov) %>%
            pivot_wider(names_from = "phase", values_from = c("meth", "cov", "avg"))
    } else {
        s_cpgs <- s_cpgs_raw %>%
            filter(tor < 0, cov >= min_cov) %>%
            mutate(avg = meth / cov) %>%
            pivot_wider(names_from = "phase", values_from = c("meth", "cov", "avg"))
    }

    print(nrow(s_cpgs))

    message("extracting bulk methylation")
    bulk_meth <- gextract_meth(c("Zhang_Nature_Genetics_2017.Ect_mCG", "Zhang_Nature_Genetics_2017.Mes_mCG", "Zhang_Nature_Genetics_2017.End_mCG"), names = c("ecto", "meso", "endo"), intervals = s_cpgs %>% select(chrom, start, end), iterator = s_cpgs %>% select(chrom, start, end), min_cov = min_cov_bulk)


    s_cpgs_f <- s_cpgs %>%
        left_join(bulk_meth %>% select(chrom:end, ecto, meso, endo)) %>%
        filter(ecto >= min_meth_bulk)

    message("downsampling")
    dsn <- min_cov
    s_cpgs_dsn <- s_cpgs_f %>%
        filter(!is.na(meth_G1)) %>%
        filter(cov_G1 >= dsn) %>%
        select(chrom, start, end, meth = meth_G1) %>%
        mutate(m = 1) %>%
        uncount(meth) %>%
        bind_rows(s_cpgs_f %>%
            filter(!is.na(meth_G1)) %>%
            filter(cov_G1 >= dsn) %>%
            mutate(unmeth = cov_G1 - meth_G1) %>%
            select(chrom, start, end, unmeth) %>%
            mutate(m = 0) %>% uncount(unmeth)) %>%
        group_by(chrom, start, end) %>%
        sample_n(dsn) %>%
        summarise(meth = sum(m), unmeth = dsn - meth, cov = dsn, .groups = "drop") %>%
        left_join(s_cpgs_f %>% select(chrom:b_score, avg_S, avg_G1))

    s_cpgs_dsn <- s_cpgs_dsn %>% mutate(avg = meth / cov)

    message(glue("# of CpGs: {nrow(s_cpgs_dsn)}"))
    return(list(g1_cpgs = s_cpgs_dsn, segmented = segmented_wt))
}
