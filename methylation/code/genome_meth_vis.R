plot_ab_diff_cis_example <- function(intervals, tracks, names, colors, genome_res, plot_res, trend_track = NULL, trend_res = 500, d_expand=NULL, min_cov=10, add_legend = TRUE, med_ab_score = gquantiles("DNMT.ab_score")) {

    scope <- intervals %>% gintervals.centers() %>% mutate(start = start - genome_res, end = end + genome_res)
    plot_scope <- intervals %>% gintervals.centers() %>% mutate(start = start - plot_res, end = end + plot_res)

    df_trend <- gextract_meth(tracks, names = names, intervals = plot_scope, d_expand = d_expand, min_cov = min_cov, annot_tracks = "DNMT.ab_score", annot_tracks_names = "ab_score") %>% select(-intervalID, -ends_with(".cov"))

    df_cpgs <- gextract_meth(tracks, names = names, intervals = plot_scope, d_expand = NULL, min_cov = min_cov, annot_tracks = "DNMT.ab_score", annot_tracks_names = "ab_score") %>% select(-intervalID, -ends_with(".cov"))
    
    df_long_trend <- df_trend %>% gather("samp", "meth", -(chrom:end), -ab_score) %>% mutate(samp = factor(samp, levels=names))
    df_long_cpgs <- df_cpgs %>% gather("samp", "meth", -(chrom:end), -ab_score) %>% mutate(samp = factor(samp, levels=names))
    

    if (!is.null(trend_track)){
        df_trend <- gextract_meth(trend_track, names = "trend", intervals = scope, iterator = trend_res, min_cov = min_cov) %>% select(-intervalID)
    }


    p_scatter <- df_long_cpgs %>% 
        filter(!is.na(meth)) %>% 
        ggplot(aes(x=start, y=meth, color=samp)) + 
            geom_rect(data = NULL, inherit.aes=FALSE, xmin = intervals$start[1], xmax = intervals$end[1], fill="yellow",  alpha=0.005, ymin=-Inf, ymax=Inf, color=NA) + 
            scale_x_continuous(limits = c(plot_scope$start, plot_scope$end), expand = c(0.01,0.01)) +
            geom_line(data = df_long_trend, size=0.5) + 
            geom_point(data = df_long_cpgs, size=0.5, alpha=0.5) + 
            scale_color_manual(name = "", values=colors) + 
            ylim(0,1) + 
            xlab(glue("{plot_scope$chrom} ({plot_scope$start}-{plot_scope$end})")) + 
            ylab("Methylation") + 
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )

    max_y_genome_high_res <- df_long_trend %>% mutate(y = ab_score - med_ab_score) %>% pull(y) %>% max(na.rm=TRUE)
    
    p_genome_high_res <- df_long_trend %>%        
    # p_genome <- gene_cors1 %>%
        distinct(start, .keep_all=TRUE) %>%
        # ggplot(aes(x = start, xend = start, y = 1, yend = 0)) +
        # geom_segment() +
        ggplot(aes(x=start, y =ab_score - med_ab_score)) + 
        geom_hline(yintercept=med_ab_score) + 
        geom_rect(data = NULL, inherit.aes=FALSE, xmin = intervals$start[1], xmax = intervals$end[1], fill="yellow",  alpha=0.005, ymin=-Inf, ymax=max_y_genome_high_res * 0.7, color=NA) +          
        # geom_point(size = 2, shape=21, fill="darkgray", color="black") +
        geom_point(size = 0.5, shape=21, aes(fill=ab_score), color="black") +
        scale_fill_gradient2(low="darkblue", high="darkred", limits=c(-1,1)) + 
        # geom_point(size = 0.5) +
        # ylim(0, 2) +
        # scale_y_continuous(limits = c(0,1), expand = c(0, 0)) + 
        # ylim(0, 1) +
        scale_x_continuous(limits = c(plot_scope$start, plot_scope$end), expand = c(0.01,0.01)) +
        xlab("") +
        ylab("3b/3a\ns-score") +
        ggpubr::theme_pubr(base_size = 4) +
        theme(
            text = element_text(family = "ArialMT", size = 6),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA), 
            plot.background = element_rect(fill = "transparent",colour = NA)
        )

    cpg_intervs <- giterator.intervals(intervals=scope, iterator="intervs.global.seq_CG")

    cpg_intervs <- cpg_intervs %>% 
        distinct(start, .keep_all=TRUE)   
    
    p_genome <- cpg_intervs %>%         
        ggplot(aes(x = start, xend = start, y = 1, yend = 0)) +
        geom_rect(data = NULL, inherit.aes=FALSE, xmin = intervals$start[1], xmax = intervals$end[1], fill="green", alpha=0.005, ymin=0, ymax=1, color=NA) + 
        geom_hline(yintercept=0) + 
        # geom_rect(xmin = plot_scope$start, xmax = plot_scope$end, fill="gray",  alpha=1, ymin=-0.5, ymax=0.5, color=NA) + 
        # scale_y_continuous(limits = c(-4,4), expand = c(0, 0)) +         
        scale_x_continuous(limits = c(scope$start, scope$end), expand = c(0,0)) +
        # geom_segment(inherit.aes = FALSE,  x = scope$start, xend = plot_scope$start, y=-4, yend = 0, linetype = "dashed") +         
        # geom_spoke(inherit.aes = FALSE,  angle = 270, linetype = "dashed", color="red", radius=0.5, x=scope$start, y=0) +         
        # geom_segment(inherit.aes = FALSE,  x = scope$end, xend = plot_scope$end, y=-4, yend = 0, linetype = "dashed") + 
        geom_segment(inherit.aes = FALSE,  x = scope$start, xend = scope$start, y=0, yend = 1) + 
        xlab("") +
        ylab("%5mc") +
        ggpubr::theme_pubr(base_size = 4) +
        theme(
            text = element_text(family = "ArialMT", size = 6),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )

   

    if (!is.null(trend_track)){
        p_genome <- p_genome + geom_line(data = df_trend, inherit.aes=FALSE, aes(x=start, y=trend), color="gray")
    }


    genes <- gintervals.neighbors1("intervs.global.tss", scope %>% select(chrom, start, end)) %>% 
        filter(dist == 0) %>%        
        distinct(geneSymbol, .keep_all = TRUE) %>%
        select(chrom, start, end, strand, geneSymbol)

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
                xend = genes$start[i] + genes$strand[i] * (scope$end - scope$start) * 0.05,
                y = 1,
                yend = 1,
                color = "blue",
                arrow = arrow(length = unit(0.05, "inches"))
            ) +
            annotate("text", label = genes$geneSymbol[i], x = genes$start[i], y = 1.7, size = 2, family = "ArialMT")
        # if (plot_all_tss_lines) {
        #     p_scatter <- p_scatter +
        #         geom_vline(color = "darkblue", xintercept = genes$start[i], linetype = "dashed")
        # }
    }

    exons <- gintervals.neighbors1("intervs.global.exons", scope %>% select(chrom, start, end)) %>% 
        filter(dist == 0) 

    p_genome <- p_genome + geom_rect(data = exons, aes(xmin=start, xmax=end, ymin=-0.2, ymax=0.2), fill="blue", color="black")

    # p_genome <- p_genome + 
    #     annotate("segment", x = scope$start, y=-0.4, xend = scope$start + 1e3, yend = -0.4, arrow = arrow(length = unit(0.05, "inches"))) + 
    #     annotate("segment", xend = scope$start, y=-0.4, x = scope$start + 1e3, yend = -0.4, arrow = arrow(length = unit(0.05, "inches"))) + 
    #     annotate("text", label = "1kbp", x = scope$start + 500, y=-0.7, size=0.5, family = "ArialMT") 
        

    # annotate('text', label = df$gene[1], x = df$start_expr[1], y = 3.5, size = 2) +
            # annotate('text', label = scales::unit_format(scale = 1e-3, unit = 'Kbp', big.mark = ',')(scope$end - scope$start), x = scope$start + 1e3, y = -0.5, size = 2) +
            # annotate('text', label = '⟷', x = scope$start + 1e3, y = -1, size = 3)

    

    # margin <- margin(t = 0, b = 0, r = 3, l = 3, unit = "pt")
   
    p_scat <- (p_genome_high_res + guides(fill="none")) / (p_scatter + guides(color="none")) + patchwork::plot_layout(heights=c(0.2, 0.8) )
    p <- p_genome / p_scat  +  patchwork::plot_layout(heights=c(0.2, 0.8) )
     # browser()
    # p_genome_comb <- cowplot::ggdraw(p_genome) + cowplot::draw_plot(p_genome_high_res + guides(fill=FALSE), height = 0.3, vjust=-0.3)
    # p_genome_comb <- p_genome_comb + 
    #     cowplot::draw_label('⟷', x=0.07, y=0.4) + 
    #     cowplot::draw_label(scales::unit_format(scale = 1e-3, unit = 'Kbp', big.mark = ',')(scope$end - scope$start), x=0.07, y=0.48)

    # p <- cowplot::plot_grid(p_genome_comb, p_scatter + guides(color = FALSE), ncol = 1, align = "v", rel_heights = c(0.6, 0.4))

    # if (add_legend){
    #     legends <- cowplot::plot_grid(cowplot::get_legend(p_scatter), cowplot::get_legend(p_genome_high_res), ncol=1)
    #     p <- cowplot::plot_grid(p, legends, nrow = 1, rel_widths = c(0.8, 0.2))    
    # }      
    

    return(p)
}
