create_vis_cfg <- function(min_cov = 10, d_expand = 100, meth_track_height = 0.3, meth_track_cex = 0.3) {
  meth_md <- tracks_key %>%
    filter(sort %in% c("meso", "endo", "epi", "epi_endo"), group == 6) %>%
    mutate(name = glue("eb_{day}_{line}_{sort}")) %>%
    select(line, day, sort, track_name, name)
  meth_md <- meth_md %>%
    bind_rows(
      ext_tracks_key %>% filter(line == "wt", sort %in% c("meso", "endo", "epi", "epi_endo"), group == "Zhang") %>% mutate(day = gsub("d", "E", day), name = glue("mouse_{day}_{sort}")) %>% select(line, day, sort, track_name, name) %>% mutate(line = "mouse")
    ) %>%
    as.data.frame()

  meth_md <- meth_md %>%
    mutate(meth_track = glue("{track_name}.meth"), cov_track = glue("{track_name}.cov"), meth_vtrack = glue("{name}.meth"), cov_vtrack = glue("{name}.cov")) %>%
    mutate(expr = glue("{meth_vtrack} / {cov_vtrack}")) %>%
    mutate(screen_expr = glue("{cov_vtrack} >= {min_cov}"))

  sort_colors <- tribble(
    ~sort, ~color,
    "meso", "darkblue",
    "endo", "darkgreen",
    "epi", "darkred",
    "ecto", "purple"
  )

  meth_md <- meth_md %>% left_join(sort_colors)

  chip_tracks <- c(
    "Xiang_Nature_Genetics_2019.Mes_H3K27ac",
    "Xiang_Nature_Genetics_2019.End_H3K27ac",
    "Xiang_Nature_Genetics_2019.Ect_H3K27ac",
    "Xiang_Nature_Genetics_2019.E65Epi_H3K27ac",
    "Xiang_Nature_Genetics_2019.E65VE_H3K27ac",
    "Xiang_Nature_Genetics_2019.PS_H3K27ac"
  )
  chip_names <- c(
    "e7_meso_k27ac",
    "e7_endo_k27ac",
    "e7_ecto_k27ac",
    "e6_epi_k27ac",
    "e6_VE_k27ac",
    "e7_PS_k27ac"
  )
  chip_exprs <- glue::glue("-log2(1 - {chip_names})")
  chip_md <- tibble(track = chip_tracks, name = chip_names, expr = chip_exprs) %>% mutate(color = "black")
  # chip_md <- chip_md %>%
  #     mutate(
  #         color = case_when(
  #             grepl("meso", name) ~ "darkblue",
  #             grepl("endo", name) ~ "darkgreen",
  #             grepl("epi", name) ~ "darkred",
  #             grepl("ecto", name) ~ "purple",
  #             TRUE ~ "black")
  #     )

  meth_md <- meth_md %>% filter(sort %in% c("meso", "endo"))
  trend_md <- meth_md %>%
    group_by(line, day) %>%
    summarise(expr = glue("pmean({expr[1]}, {expr[2]})"), screen_expr = glue("{cov_vtrack[1]} >= {min_cov} & {cov_vtrack[2]} >= {min_cov}")) %>%
    mutate(name = glue("{day}_{line}")) %>%
    ungroup() %>%
    mutate(name = factor(name, levels = c("d6_wt", "d6_ko3a", "d6_ko3b", "E7_mouse"))) %>%
    arrange(name) %>%
    mutate(name = as.character(name))


  d4_tracks <- tracks_key %>%
    filter(line %in% c("wt", "ko3a", "ko3b"), day == "d4") %>%
    mutate(meth_track = glue("{track_name}.meth"), cov_track = glue("{track_name}.cov"), meth_vtrack = glue("{name}.meth"), cov_vtrack = glue("{name}.cov"))

  d4_joined <- d4_tracks %>%
    group_by(line, day) %>%
    summarise(meth_expr = glue("psum({tracks}, na.rm=TRUE)", tracks = paste(meth_vtrack, collapse = ", ")), cov_expr = glue("psum({tracks}, na.rm=TRUE)", tracks = paste(cov_vtrack, collapse = ", ")), expr = glue("{meth_expr} / {cov_expr}")) %>%
    unite("name", day, line) %>%
    mutate(name = factor(name, levels = c("d4_wt", "d4_ko3a", "d4_ko3b"))) %>%
    arrange(name) %>%
    mutate(name = as.character(name)) %>%
    mutate(screen_expr = glue("{cov_expr} >= {min_cov}"))

  trend_line <- trend_md[4, ]

  v <- vis_create() %>%
    vis_add_vtrack(vtrack = c(meth_md$meth_vtrack, d4_tracks$meth_vtrack), src = c(meth_md$meth_track, d4_tracks$meth_track), sshift = -d_expand, eshift = d_expand, func = "sum") %>%
    vis_add_vtrack(vtrack = c(meth_md$cov_vtrack, d4_tracks$cov_vtrack), src = c(meth_md$cov_track, d4_tracks$cov_track), sshift = -d_expand, eshift = d_expand, func = "sum") %>%
    vis_add_ideogram(height = 0.5) %>%
    vis_add_genome_axis(height = 0.5, fontsize = 5) %>%
    vis_add_genes(height = 1.5, fonsize.group = 2) %>%
    vis_add_track(track = c(trend_line$expr, meth_md$expr[1:2]), screen_expr = c(trend_line$screen_expr, meth_md$screen_expr[1:2]), iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", meth_md$color[1:2]), cex = meth_track_cex, group_tracks = TRUE, legend = FALSE, plot_type = list(c("l"), c("p"), c("p")), group_name = trend_md$name[1], ylim = list(c(0, 1)), height = meth_track_height) %>%
    vis_add_track(track = c(trend_line$expr, meth_md$expr[3:4]), screen_expr = c(trend_line$screen_expr, meth_md$screen_expr[3:4]), iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", meth_md$color[3:4]), cex = meth_track_cex, group_tracks = TRUE, legend = FALSE, plot_type = list(c("l"), c("p"), c("p")), group_name = trend_md$name[2], ylim = list(c(0, 1)), height = meth_track_height) %>%
    vis_add_track(track = c(trend_line$expr, meth_md$expr[5:6]), screen_expr = c(trend_line$screen_expr, meth_md$screen_expr[5:6]), iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", meth_md$color[5:6]), cex = meth_track_cex, group_tracks = TRUE, legend = FALSE, plot_type = list(c("l"), c("p"), c("p")), group_name = trend_md$name[3], ylim = list(c(0, 1)), height = meth_track_height) %>%
    vis_add_track(track = c(trend_line$expr, meth_md$expr[7:8]), screen_expr = c(trend_line$screen_expr, meth_md$screen_expr[7:8]), iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", meth_md$color[7:8]), cex = meth_track_cex, group_tracks = TRUE, legend = FALSE, plot_type = list(c("l"), c("p"), c("p")), group_name = trend_md$name[4], ylim = list(c(0, 1)), height = meth_track_height) %>%
    vis_add_track(track = d4_joined$expr[1:2], screen_expr = d4_joined$screen_expr[1:2], iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", "darkred"), cex = meth_track_cex, legend = FALSE, plot_type = list(c("l", "p")), name = d4_joined$name[1:2], group_name = "d4_ko3a", ylim = list(c(0, 1)), height = meth_track_height, group_tracks = TRUE) %>%
    vis_add_track(track = d4_joined$expr[c(1, 3)], screen_expr = d4_joined$screen_expr[c(1, 3)], iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = c("black", "darkred"), cex = meth_track_cex, legend = FALSE, plot_type = list(c("l", "p")), group_name = "d4_ko3b", name = d4_joined$name[c(1, 3)], ylim = list(c(0, 1)), height = meth_track_height, group_tracks = TRUE) %>%
    vis_add_vtrack(vtrack = chip_md$name, src = chip_md$track, func = "global.percentile") %>%
    vis_add_track(track = chip_md$expr[1:4], name = chip_md$name[1:4], plot_type = "histogram", iterator = "dense", color = chip_md$color[1:4], ylim = list(c(3, 12)), grid = FALSE, height = 0.3) %>%
    vis_save(here("output/endo_meso_vis.yaml"))
  return(v)
}

plot_vis_meso_endo_intervals <- function(...){
    plot_vis_intervals(..., vis_cfg_func = create_vis_cfg)
}


plot_vis_intervals <- function(intervals, dir, filename, vis_cfg_func, zooms = c(5, 8), meth_track_cex = c(0.1, 0.3, 0.3), min_cov = 10, d_expand = 100, plot_height = 4, plot_width = 7, hl_d_expand = 100) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(intervals, filename)

    plot_intervals <- function(intervs) {
        hl_intervs <- intervs %>% mutate(start = start - hl_d_expand, end = end + hl_d_expand)
        try({            
            v <- vis_cfg_func(min_cov = min_cov, d_expand = d_expand, meth_track_cex = meth_track_cex[1]) %>%
                vis_add_title(paste(intervs$geneSymbol[1], glue("diff={round(intervs$diff, digits=2)}"))) %>%
                vis_modify(intervals = intervs %>% mutate(start = start - 1e5, end = end + 1e5), highlight = hl_intervs) %>%
                vis_save_plot(glue("{dir}/{intervs$geneSymbol[1]}_{intervs$chrom}_{intervs$start}_{intervs$end}.png"), width = plot_width, height = plot_height)
            for (i in 1:length(zooms)) {
                v <- vis_cfg_func(min_cov = min_cov, d_expand = d_expand, meth_track_cex = meth_track_cex[i + 1]) %>%
                    vis_add_title(paste(intervs$geneSymbol[1], glue("diff={round(intervs$diff, digits=2)}"))) %>%
                    vis_modify(intervals = intervs %>% mutate(start = start - 1e5, end = end + 1e5), highlight = hl_intervs) %>%
                    vis_zoom_in(zooms[i]) %>%
                    vis_save_plot(glue("{dir}/{intervs$geneSymbol[1]}_{intervs$chrom}_{intervs$start}_{intervs$end}_zoom{zooms[i]}.png"), width = plot_width, height = plot_height)
            }
        })
    }

    # plot the first interval in order to cache virtual tracks
    plot_intervals(intervals[1, ])

    plyr::a_ply(intervals[-1, ], 1, plot_intervals, .parallel = TRUE)
}

create_vis_cfg_germ_layers_diff <- function(min_cov = 10, d_expand = 100, meth_track_height = 0.3, meth_track_cex = 0.3) {    
    # meth_md_eb <- tracks_key %>%
    #     filter(sort %in% c("meso", "endo", "epi", "epi_endo"), group == 6) %>%
    #     mutate(name = glue("eb_{day}_{line}_{sort}")) %>%
    #     select(line, day, sort, track_name, name)
  
    meth_md_gl <- tribble(
        ~track_name, ~name, ~sort,        
        "Zhang_Nature_Genetics_2017.Ect_mCG", "ecto_d7.5", "ecto",        
        "mEBDNMT.PBAT000215_wt_emb", "zm_ecto_d7.5", "ecto",
        "Zhang_Nature_Genetics_2017.Mes_mCG", "meso_d7.5", "meso",
        "mEBDNMT.PBAT000213_wt_emb", "zm_meso_d7.5", "meso",
        "mEBDNMT.PBAT000214_wt_emb", "zm_endo_d7.5", "endo",
        "Zhang_Nature_Genetics_2017.E65Epi_mCG", "epi_d6.5", "epi",
        "ENCODE.Ecker.e10_5.forebrain.wgbs_sum", "forebrain_d10.5", "brain",
        "ENCODE.Ecker.e10_5.heart.wgbs_sum", "heart_d10.5", "heart",
        "ENCODE.Ecker.e10_5.limb.wgbs_sum", "limb_d10.5", "limb",
        "ENCODE.Ecker.e11_5.liver.wgbs_sum", "liver_d11.5", "liver",
        "mEBDNMT.PBAT000192_wt_d6", "eb_meso_d6", "meso",
        "mEBDNMT.PBAT000193_wt_d6", "eb_endo_d6", "endo",
        "mEBDNMT.PBAT000181_wt_d6", "eb_epi_endo_d6", "epi_endo",
        "mEBDNMT.PBAT000194_wt_d6", "eb_epi_d6", "epi",
        "mEBDNMT.PBAT000207_wt_d5", "eb_endo_d5", "endo",
        "mEBDNMT.PBAT000208_wt_d5", "eb_epi_d5", "epi"
    )

    # meth_md <- bind_rows(meth_md_gl, meth_md_eb)
    meth_md <- meth_md_gl %>%
        mutate(meth_track = glue("{track_name}.meth"), cov_track = glue("{track_name}.cov"), meth_vtrack = glue("{name}.meth"), cov_vtrack = glue("{name}.cov")) %>%
        mutate(expr = glue("{meth_vtrack} / {cov_vtrack}")) %>%
        mutate(screen_expr = glue("{cov_vtrack} >= {min_cov}"))
    

    sort_colors <- tribble(
        ~sort, ~color,
        "meso", "darkblue",
        "endo", "darkgreen",
        "epi", "darkred",
        "epi_endo", "darkgreen",
        "ecto", "purple",
        "brain", "purple",
        "heart", "darkblue",
        "liver", "darkgreen",
        "limb", "darkblue"
    )

    meth_md <- meth_md %>% left_join(sort_colors)

    chip_tracks <- c(
        "Xiang_Nature_Genetics_2019.Mes_H3K27ac",
        "Xiang_Nature_Genetics_2019.End_H3K27ac",
        "Xiang_Nature_Genetics_2019.Ect_H3K27ac",
        "Xiang_Nature_Genetics_2019.E65Epi_H3K27ac",
        "Xiang_Nature_Genetics_2019.E65VE_H3K27ac",
        "Xiang_Nature_Genetics_2019.PS_H3K27ac"
    )
    chip_names <- c(
        "e7_ecto_k27ac",
        "e7_meso_k27ac",
        "e7_endo_k27ac",        
        "e6_epi_k27ac",
        "e6_VE_k27ac",
        "e7_PS_k27ac"
    )
    chip_exprs <- glue::glue("-log2(1 - {chip_names})")
    chip_md <- tibble(track = chip_tracks, name = chip_names, expr = chip_exprs) %>% mutate(color = "black")

    v <- vis_create() %>%
        vis_add_vtrack(vtrack = c(meth_md$meth_vtrack), src = c(meth_md$meth_track), sshift = -d_expand, eshift = d_expand, func = "sum") %>%
        vis_add_vtrack(vtrack = c(meth_md$cov_vtrack), src = c(meth_md$cov_track), sshift = -d_expand, eshift = d_expand, func = "sum") %>%
        vis_add_ideogram(height = 0.5) %>%
        vis_add_genome_axis(height = 0.5, fontsize = 5) %>%
        vis_add_genes(height = 1.5, fonsize.group = 2) %>%
        vis_add_track(track  = meth_md$expr, screen_expr = meth_md$screen_expr, color = meth_md$color, cex = meth_track_cex, height = meth_track_height, ylim = list(c(0, 1)), plot_type="p", name = meth_md$name) %>%    
        vis_add_vtrack(vtrack = chip_md$name, src = chip_md$track, func = "global.percentile") %>%
        vis_add_track(track = chip_md$expr[1:4], name = chip_md$name[1:4], plot_type = "histogram", iterator = "dense", color = chip_md$color[1:4], ylim = list(c(3, 12)), grid = FALSE, height = 0.3) %>%
        vis_save(here("output/germ_layer_vis.yaml"))
    return(v)
}

vis_add_meth_track <- function(v, tracks, name="meth", min_cov=NULL, d_expand = NULL, iterator = "intervs.global.seq_CG", screen_iterator = "intervs.global.seq_CG", color = "black", height = 0.3, cex = 0.3, ylim=list(c(0, 1)), plot_type = 'p', group_tracks = TRUE, group_name = NULL, plot_cov = FALSE, cov_name = "cov", ylim_cov = c(0, 250), ...){
    
    if (!is.list(tracks)){
        tracks_list <- list(tracks)
    } else {
        tracks_list <- tracks
    }

    exprs <- c()
    cov_exprs <- c()
    screen_exprs <- c()

    for (i in 1:length(tracks_list)){
        tracks <- tracks_list[[i]]

        cov_vtracks <- glue("{tracks}_smoo.cov")
        walk2(cov_vtracks, glue("{tracks}.cov"), ~ gvtrack.create(.x, .y, "sum"))
        meth_vtracks <- glue("{tracks}_smoo.meth")
        walk2(meth_vtracks, glue("{tracks}.meth"), ~ gvtrack.create(.x, .y, "sum"))

        if (!is.null(d_expand)) {
            walk(cov_vtracks, ~ gvtrack.iterator(.x, sshift = -d_expand, eshift = d_expand))
            walk(meth_vtracks, ~ gvtrack.iterator(.x, sshift = -d_expand, eshift = d_expand))
        }


        meth_expr <- glue("psum({tracks}, na.rm=TRUE)", tracks = paste(meth_vtracks, collapse = ", "))
        cov_expr <- glue("psum({tracks}, na.rm=TRUE)", tracks = paste(cov_vtracks, collapse = ", "))
        exprs <- c(exprs, glue("{meth_expr} / {cov_expr}"))
        cov_exprs <- c(cov_exprs, cov_expr)
        

        if (!is.null(min_cov)){
            screen_exprs <- c(screen_exprs, glue("{cov_expr} >= {min_cov}"))
        }

        if (!is.null(d_expand)){
            v <- v %>% vis_add_vtrack(vtrack = glue("{tracks}_smoo.meth"), src = glue("{tracks}.meth"), sshift = -d_expand, eshift = d_expand, func = "sum")
            v <- v %>% vis_add_vtrack(vtrack = glue("{tracks}_smoo.cov"), src = glue("{tracks}.cov"), sshift = -d_expand, eshift = d_expand, func = "sum")
        } else {
            v <- v %>% vis_add_vtrack(vtrack = glue("{tracks}_smoo.meth"), src = glue("{tracks}.meth"), func = "sum")
            v <- v %>% vis_add_vtrack(vtrack = glue("{tracks}_smoo.cov"), src = glue("{tracks}.cov"), func = "sum")
        }    
    }
    

    v <- v %>% vis_add_track(track  = exprs, screen_expr = screen_exprs, screen_iterator = screen_iterator, color = color, cex = cex, height = height, ylim = ylim, plot_type="p", name = name, iterator=iterator, group_tracks = group_tracks, group_name = group_name, ...)

    if (plot_cov){        
        v <- v %>% vis_add_track(track  = cov_exprs,  color = color, cex = cex, height = height, ylim = list(ylim_cov), plot_type="p", name = name, iterator=iterator, group_tracks = group_tracks, group_name = cov_name, ...)
    }

    return(v)
}

create_vis_cfg_3a_3b_diff <- function(min_cov = 10, d_expand = 250){
    df_wt <- tracks_key  %>% filter(day == "d6") %>% filter(line == "wt")
    df_ko3a <- tracks_key  %>% filter(day == "d6") %>% filter(line == "ko3a")
    df_ko3b <- tracks_key  %>% filter(day == "d6") %>% filter(line == "ko3b")

    df_wt_d5 <- tracks_key  %>% filter(day == "d5") %>% filter(line == "wt")
    df_ko3a_d5 <- tracks_key  %>% filter(day == "d5") %>% filter(line == "ko3a")
    df_ko3b_d5 <- tracks_key  %>% filter(day == "d5") %>% filter(line == "ko3b")

    df_wt_d1_to_4 <- tracks_key  %>% filter(day %in% paste0("d", 1:4)) %>% filter(line == "wt")
    df_ko3a_d1_to_4 <- tracks_key  %>% filter(day %in% paste0("d", 1:4)) %>% filter(line == "ko3a")
    df_ko3b_d1_to_4 <- tracks_key  %>% filter(day %in% paste0("d", 1:4)) %>% filter(line == "ko3b")

    invivo_tracks <- c(
        "Zhang_Nature_Genetics_2017.Ect_mCG",
        "Zhang_Nature_Genetics_2017.Mes_mCG",
        "Zhang_Nature_Genetics_2017.End_mCG")

    invivo_names <- c("ecto", "meso", "endo")

    chip_tracks <- c(
        "Xiang_Nature_Genetics_2019.Mes_H3K27ac",
        "Xiang_Nature_Genetics_2019.End_H3K27ac",
        "Xiang_Nature_Genetics_2019.Ect_H3K27ac",
        "Xiang_Nature_Genetics_2019.E65Epi_H3K27ac",
        "Xiang_Nature_Genetics_2019.E65VE_H3K27ac",
        "Xiang_Nature_Genetics_2019.PS_H3K27ac"
    )
    chip_names <- c(
        "e7_ecto_k27ac",
        "e7_meso_k27ac",
        "e7_endo_k27ac",        
        "e6_epi_k27ac",
        "e6_VE_k27ac",
        "e7_PS_k27ac"
    )
    chip_exprs <- glue::glue("-log2(1 - {chip_names})")
    chip_md <- tibble(track = chip_tracks, name = chip_names, expr = chip_exprs) %>% mutate(color = "black")

    v <- vis_create() %>%
        vis_add_ideogram(height = 0.2) %>%
        vis_add_genome_axis(height = 0.5, fontsize = 5) %>%
        vis_add_genes(height = 0.5, fonsize.group = 2) %>% 
        vis_add_meth_track(list(df_wt$track_name, df_ko3a$track_name, df_ko3b$track_name), min_cov = min_cov, d_expand = 250, name = c("wt", "ko3a", "ko3b"), color = c("blue", "darkgreen", "darkred"), group_name = "d6", cex = 0.1, plot_cov = FALSE, legend = FALSE) %>%                 
        vis_add_meth_track(as.list(invivo_tracks), min_cov = min_cov, d_expand = 250, name = invivo_names, color = c("blue", "#1B1919FF", "#ED0000FF"), group_name = "invivo", cex = 0.1, legend = FALSE) %>%         
        vis_add_meth_track(list(df_wt_d5$track_name, df_ko3a_d5$track_name, df_ko3b_d5$track_name), min_cov = min_cov, d_expand = 250, name = c("wt", "ko3a", "ko3b"), color = c("blue", "darkgreen", "darkred"), group_name = "d5", cex = 0.1, legend = FALSE) %>%         
        # vis_add_meth_track(list(df_wt_d1_to_4$track_name, df_ko3a_d1_to_4$track_name, df_ko3b_d1_to_4$track_name), min_cov = min_cov, d_expand = 250, name = c("wt", "ko3a", "ko3b"), color = c("blue", "darkgreen", "darkred"), group_name = "d1-4", cex = 0.1, legend = FALSE) %>%    
        vis_add_track(track = "DNMT.ab_score", name = "AB score", plot_type = list("p"), iterator = "intervs.global.seq_CG", color = "black", height = 0.5, cex = 0.1) %>% 
        vis_add_vtrack(vtrack = chip_md$name, src = chip_md$track, func = "global.percentile") %>%
        vis_add_track(track = chip_md$expr[1:4], name = chip_md$name[1:4], plot_type = "histogram", iterator = "dense", color = chip_md$color[1:4], ylim = list(c(3, 12)), grid = FALSE, height = 0.3)
        # vis_add_meth_track(df_wt$track_name, min_cov = min_cov, d_expand = 250, name = "wt") %>% 
        # vis_add_meth_track(df_ko3a$track_name, min_cov = min_cov, d_expand = 250, name = "ko3a") %>%  
        # vis_add_meth_track(df_ko3b$track_name, min_cov = min_cov, d_expand = 250, name = "ko3b")

    return(v)
}




