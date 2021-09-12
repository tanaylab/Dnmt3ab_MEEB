get_mc_proj2d <- function(mc){
    graph <- mc$tads_clust$coclust$co_cluster %>% rename(cell1 = node1, cell2 = node2)    
    proj2d <- sc5mc.mc2d_force_knn(mc$mc_map, graph)
    return(proj2d)
}


plot_metacell_avg_meth <- function(mc, low_cgc){
    d <- get_mc_map(mc) %>% left_join(low_cgc)
    d <-  d %>% abbreviate_mc() %>% mutate(meth = cut(low_cgc_score, 9) ) %>% count(mc, meth) %>% group_by(mc) %>% mutate(p = n / sum(n))
    p <- d %>% 
        mutate(meth = forcats::fct_rev(meth)) %>% 
        ggplot(aes(x=mc, y=p, fill=meth)) + 
            geom_col() + 
            scale_y_continuous(labels=scales::percent) + 
            xlab('') + 
            ylab('') + 
            scale_fill_brewer(palette='RdBu') + 
            theme(strip.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())    

    return(p)
}

plot_mc_map_summary <- function(mc_map, field){    
    field <- rlang::sym(field)
    colors <- config$color_conf[[field]] %||% c(categorical_colors, categorical_colors)   
    
    mc_map <- abbreviate_mc(mc_map)
    mc_map <- mc_map %>% summarise_mc_map(field)
    
    n_data <- mc_map %>% distinct(mc, n)
    
    p <- mc_map %>% 
        ggplot(aes(x=mc, y=p, fill=!!field)) + 
            geom_col() + 
            scale_fill_manual(values=colors) + 
            scale_y_continuous(labels=scales::percent) + 
            ylab('') + 
            xlab('') + 
            geom_text(data=n_data, inherit.aes=FALSE, aes(x=mc, label=n, y=-0.05), size=3) + 
            guides(fill=guide_legend(ncol=2))

    return(p)
}

summarise_mc_map <- function(mc_map, field){    
    field <- rlang::sym(field)
    
    field_levels <- names(config$color_conf[[field]])
    if (is.null(field_levels)){
        field_levels <- unique(mc_map[[field]])
    }    

    # cells_n <- db@cells %>% count(!!field) %>% mutate(p_t = n / sum(n)) %>% select(!!field, p_t)    
    
    mc_map <- mc_map %>% 
        count(mc, !!field) %>% 
        # left_join(cells_n) %>% 
        group_by(mc) %>% 
        # mutate(p = (n / p_t) / sum(n / p_t), n = sum(n)) %>%
        mutate(p = n / sum(n), n = sum(n)) %>%
        ungroup() %>% 
        mutate(!!field := factor(!! field, levels=field_levels)) %>% 
        factorize_mc()

    return(mc_map)
}

plot_proj2d_sources_annot <- function(proj2d, mc_map, annot){
    plots <- map(unique(mc_map$cell_source), ~ 
        sc5mc.plot_proj2d(proj2d, alpha=0.5, point_color='gray', point_size=0.5, plot_mc_points=FALSE, plot_graph=FALSE, plot_mc_names=FALSE) + 
        geom_point(data=proj2d$cells %>% filter(cell_source == .x), inherit.aes=FALSE, aes_string(x='x', y='y', color=annot), size=0.5) + scale_color_manual(values=color_conf[[annot]], guide=FALSE) + ggtitle(.x))
    p <- cowplot::plot_grid(plotlist=plots, align='hv')
    return(p)
}

plot_metacell_matrices <- function(mc, outdir, tads_cor=NULL, hc_tads=NULL, res_scale=1, raster_quality=1){
    if (is.null(mc$tads_vm)){
        mc$tads_vm <- get_mc_tads_vm(mc$tads_ds)
    }

    list2env(mc$tads_clust, envir=environment())
    
    mc_map <- mc$mc_map
    
    annot_df <- data.frame(cell_id = colnames(tads_m)) %>% 
        left_join(db@cells %>% select(cell_id, cell_source, cell_type, treatment)) %>% 
        left_join(mc_map %>% select(cell_id)) %>%
        select(-cell_id)

    
    ha <- ComplexHeatmap::rowAnnotation(df = annot_df)#, col=color_conf)

    colors <- circlize::colorRamp2(seq(quantile(tads_m, 0.05, na.rm=TRUE), quantile(tads_m, 0.95, na.rm=TRUE),length.out=5), c("#0571b0","#92c5de","white","#f4a582","#ca0020"))

    split_mc <- tibble(cell_id = colnames(tads_m)) %>% left_join(mc_map) %>% pull(mc)
    

    tads_cor <- tads_cor %||% tgs_cor(t(tads_m), spearman=TRUE, pairwise.complete.obs=TRUE) 
    hc_tads <-  hc_tads %||% (tads_cor %>% tgs_dist() %>% hclust(method='ward.D2'))

    hm <- ComplexHeatmap::Heatmap(t(tads_m), name = 'Tad methylation', show_column_names=FALSE, show_row_names=FALSE, cluster_columns=hc_tads, cluster_rows=TRUE, split = split_mc, column_title = glue('{ncol(tads_m)} cells X {nrow(tads_m)} TADS'), row_dend_reorder = FALSE, col=colors, column_dend_reorder=FALSE, use_raster=TRUE, raster_quality=raster_quality) 

    shades <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(1000)

    cell_ord <- ComplexHeatmap::row_order(hm) %>% flatten() %>% as_vector()

    hm_coclust <- ComplexHeatmap::Heatmap(m_coc[cell_ord, cell_ord], cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE, col=c("white", "pink", "red", "black", "brown", "orange"))
    
    hm_cor <-  ComplexHeatmap::Heatmap(cm[cell_ord, cell_ord], cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE, use_raster=TRUE, raster_quality=raster_quality)
    
    if (!is.null(cm_norm)){
        hm_cor_norm <-  ComplexHeatmap::Heatmap(cm_norm[cell_ord, cell_ord], cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE, use_raster=TRUE, raster_quality=raster_quality)

        # norm_knn_mat <- norm_knn_mat[rownames(cm_norm), ]
        # hm_cell_cycle <- ComplexHeatmap::Heatmap(norm_knn_mat[cell_ord, ], show_row_names=FALSE, show_column_names=FALSE, cluster_columns = FALSE, cluster_rows = FALSE, width = unit(5, "cm"))
    } 

    m_knn <- cm_knn %>% select(col1, col2, val) %>% tidyr::complete(col1, col2) %>% as.data.frame() %>% spread(col2, val) %>% column_to_rownames('col1') %>% as.matrix()
    hm_knn <- ComplexHeatmap::Heatmap(m_knn[cell_ord, cell_ord], cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE, use_raster=TRUE, raster_quality=raster_quality, , na_col='white', col=shades)

    m_g <- g %>% select(col1, col2, weight) %>% tidyr::complete(col1, col2) %>% as.data.frame() %>% spread(col2, weight) %>% column_to_rownames('col1') %>% as.matrix()
    hm_g <- ComplexHeatmap::Heatmap(m_g[cell_ord, cell_ord], cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE, use_raster=TRUE, raster_quality=raster_quality, na_col='white', col=shades)
    
    message('printing plots')
    png(glue('{outdir}/tads_mat.png'), width = 900*res_scale, height = 1100*res_scale)
    ComplexHeatmap::draw(hm + ha)
    dev.off()

    png(glue('{outdir}/matrices.png'), width = 4000*res_scale, height = 1100*res_scale)
    if (!is.null(cm_norm)){
        # ComplexHeatmap::draw(hm_cor + hm_cell_cycle + hm_cor_norm + hm_knn + hm_g + hm_coclust)        
        ComplexHeatmap::draw(hm_cor + hm_cor_norm + hm_knn + hm_g + hm_coclust)        
    } else {
        ComplexHeatmap::draw(hm_cor + hm_knn + hm_g + hm_coclust)        
    }

    dev.off()


}

