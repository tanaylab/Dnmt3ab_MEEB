

#' Screen for "dips" in methylation
#' 
#' @param track methylation track prefix (without the ".cov" or ".meth" suffix)
#' @param intervals intervals set 
#' @param iterator iterator. Default - CpGs
#' @param dip_size size of the dip region 
#' @param flank_size size of flank regions 
#' @param margin margin between the dip region and the flanking regions
#' @param min_diff minimal average methylation difference between dip regions and flanking regions. A "dip" would be called only if a difference exists between the center and both right and left flanks.
#' @param min_cov minimal coverage for center
#' @param min_nhood_cov minimal coverage for flanks
#' @param canonize merge overlapping regions in the resulting intervals
#' @param ret_values if FALSE - returns only the "dip" intervals. if TRUE - returns also the "cov" and "meth" values for the center and flanking regions with intervals shifted by dip_size/2.
#' 
#' 
meth_dip_screen <- function(track, intervals = gintervals.all(), iterator="intervs.global.seq_CG", dip_size=300, flank_size=500, margin=50, min_diff=0.5, min_cov=10, min_nhood_cov=min_cov, canonize=FALSE, ret_values = FALSE){
    cov_track <- glue("{track}.cov")
    meth_track <- glue("{track}.meth")
    
    center_shift <- dip_size / 2
    
    gvtrack.create("cov_vt_center", cov_track, func="sum")
    gvtrack.create("cov_vt_right", cov_track, func="sum")
    gvtrack.create("cov_vt_left", cov_track, func="sum", )
    gvtrack.iterator("cov_vt_center", sshift = -center_shift, eshift = center_shift)
    gvtrack.iterator("cov_vt_right", sshift = -flank_size -center_shift - margin, eshift = -center_shift - 1 - margin)
    gvtrack.iterator("cov_vt_left", sshift = center_shift + 1 + margin, eshift = center_shift + flank_size + margin)

    gvtrack.create("meth_vt_center", meth_track, func="sum")
    gvtrack.create("meth_vt_left", meth_track, func="sum")
    gvtrack.create("meth_vt_right", meth_track, func="sum")
    gvtrack.iterator("meth_vt_center", sshift = -center_shift, eshift = center_shift)
    gvtrack.iterator("meth_vt_right", sshift = -flank_size -center_shift - 1 - margin, eshift = -center_shift - 1 - margin)
    gvtrack.iterator("meth_vt_left", sshift = center_shift + 1 + margin, eshift = center_shift + flank_size + 1 + margin)
    
    center_expr <- "(meth_vt_center / cov_vt_center)"
    left_expr <- "(meth_vt_left / cov_vt_left)"
    right_expr <- "(meth_vt_right / cov_vt_right)"        
    
    left_screen_expr <- glue("{left_expr} - {center_expr} >= {min_diff}")
    right_screen_expr <- glue("{right_expr} - {center_expr} >= {min_diff}")
    cov_screen_expr <- glue("cov_vt_center >= {min_cov} & cov_vt_left >= {min_nhood_cov} & cov_vt_right >= {min_nhood_cov}")
    screen_expr <- glue("{left_screen_expr} & {right_screen_expr} & {cov_screen_expr}")

    df <- gscreen(screen_expr, iterator = iterator, intervals = intervals)
    
    if (canonize){
        df_canonic <- df %>% mutate(start = start - center_shift, end = end + center_shift) %>% gintervals.force_range() %>% gintervals.canonic()
        df <- df %>% 
            gintervals.neighbors1(df_canonic) %>% 
            gintervals.neighbors1(df_canonic %>% gintervals.centers()) %>% 
            arrange(chrom1, start1, end1, dist1) %>% 
            group_by(chrom1, start1, end1) %>% 
            slice(1) %>% 
            ungroup() %>% 
            select(chrom, start, end)
    }
    
    if (ret_values){
        df <- gextract(c(center_expr, left_expr, right_expr, "cov_vt_center", "cov_vt_left", "cov_vt_right", "meth_vt_center", "meth_vt_left", "meth_vt_right"), iterator=df %>% select(chrom, start, end), intervals = df %>% select(chrom, start, end), colnames = c("center", "left", "right", "center.cov", "left.cov", "right.cov", "center.meth", "left.meth", "right.meth")) %>% arrange(intervalID) %>% select(-intervalID) %>% as_tibble()
        df <- df  %>% 
            mutate(start = start - center_shift, end = end + center_shift) 
    }
    
    return(df)
}

