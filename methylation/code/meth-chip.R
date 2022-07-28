extract_meth_annot <- function(tracks, names, chip_tracks = NULL, ...) {
    define_chip_vtracks()
    # k4me3_peaks = define_k4me3_peaks() %cache_df% here("output/k4me3_peaks.tsv")
    # k27me3_peaks = define_k27me3_peaks() %cache_df% here("output/k27me3_peaks.tsv")
    # gvtrack.create("d_k4me3", k4me3_peaks, "distance")
    # gvtrack.create("d_k27me3", k27me3_peaks, "distance")
    gvtrack.create("d_exon", "intervs.global.exons", "distance")
    gvtrack.create("d_tss", "intervs.global.tss", "distance")
    gvtrack.create("tor", "Encode.esd3.replichip.rep2")
    gvtrack.create("ab_score", "DNMT.ab_score")
    gvtrack.create("a_score", "DNMT.a_score")
    gvtrack.create("b_score", "DNMT.b_score")
    gvtrack.iterator("tor", sshift = -15000, eshift = 15000)
    gvtrack.create("cg_cont", "seq.CG_500_mean")
    gvtrack.create("gc_cont", "seq.GC_500_mean")

    annot_tracks <- c("d_exon", "d_tss", "tor", "ab_score", "a_score", "b_score", "cg_cont", "gc_cont")

    annot_names <- annot_tracks
    if (!is.null(annot_tracks)) {
        annot_tracks <- c(annot_tracks, glue("-log2(1-{chip_tracks})"))
        annot_names <- c(annot_names, chip_tracks)
    }

    d <- gextract_meth(tracks = tracks, names = names, annot_tracks = annot_tracks, annot_tracks_names = annot_names, ...)
    return(d)
}
