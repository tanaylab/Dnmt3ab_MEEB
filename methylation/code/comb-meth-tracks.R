comb_line_day_tracks <- function() {
    gdir.create("mEBDNMT/comb", showWarnings = FALSE)

    plyr::ddply(tracks_key %>% filter(line != "mouse"), c("day", "line"), function(x) {
        new_track <- glue("mEBDNMT.comb.{x$day[1]}_{x$line[1]}")
        print(new_track)
        comb_meth_tracks(x$track_name, new_track)
    }, .parallel = TRUE)
    gdb.reload()
}

comb_meth_tracks <- function(tracks, new_track, d_expand = NULL) {
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
    expr <- glue("{meth_expr} / {cov_expr}")

    iter <- gscreen(glue("{cov_expr} > 0"), intervals = gintervals.all(), iterator = "intervs.global.seq_CG") %>% arrange(chrom, start, end)

    gdir.create(gsub("\\.", "/", new_track), showWarnings = FALSE)

    gtrack.create(expr = cov_expr, track = paste0(new_track, ".cov"), iterator = iter, description = glue("combind tracks: {trs}", trs = paste(tracks, collapse = ", ")))
    gtrack.create(expr = meth_expr, track = paste0(new_track, ".meth"), iterator = iter, description = glue("combind tracks: {trs}", trs = paste(tracks, collapse = ", ")))
    gtrack.create(expr = expr, track = paste0(new_track, ".avg"), iterator = iter, description = glue("combind tracks: {trs}", trs = paste(tracks, collapse = ", ")))
}
