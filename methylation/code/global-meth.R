# calc_global_meth <- function() {
#     tracks <- gpatterns::gpatterns.ls("EarlyDevZME000003")
#     enh_intervs <- define_enhancer_intervs() %cache_df% here("output/enh_intervs.csv") %>% as_tibble()
#     prom_intervs <- gpatterns::gpatterns.get_promoters(upstream = 1e3, downstream=50) %>% as_tibble()
#     off_tar_intervs <- gintervals.diff(gintervals.all(), "intervs.captPBAT_probes.ES_EB_V1") %>% gintervals.diff(enh_intervs) %>% gintervals.diff(prom_intervs)

#     gmeth_low <- gpatterns::gpatterns.global_meth_trend(tracks = tracks, intervals = off_tar_intervs, strat_breaks = c(0.01, 0.02, 0.03))$trend
#     gmeth_low_high <- gpatterns::gpatterns.global_meth_trend(tracks = tracks, intervals = off_tar_intervs, strat_breaks = c(0, 0.03, 0.08, 0.15))$trend
#     gmeth <- bind_rows(gmeth_low, gmeth_low_high) %>%
#         select(track, breaks, cg_num, meth) %>%
#         left_join(get_sample_annotation() %>% select(track, samp_id, eb_day, mut))
#     return(gmeth)
# }

calc_global_meth <- function() {
    tracks <- tracks_key$track_name
    enh_intervs <- define_enhancer_intervs() %cache_df% here("output/enh_intervs.csv") %>% as_tibble()
    prom_intervs <- gpatterns::gpatterns.get_promoters(upstream = 1e3, downstream = 50) %>% as_tibble()
    off_tar_intervs <- gintervals.diff(gintervals.all(), "intervs.captPBAT_probes.ES_EB_V1") %>%
        gintervals.diff(enh_intervs) %>%
        gintervals.diff(prom_intervs) %>%
        gintervals.diff("intervs.captPBAT_probes.ES_EB_V2")

    gmeth_low <- gpatterns::gpatterns.global_meth_trend(tracks = tracks, intervals = off_tar_intervs, strat_breaks = c(0.01, 0.02, 0.03))$trend
    gmeth_low_high <- gpatterns::gpatterns.global_meth_trend(tracks = tracks, intervals = off_tar_intervs, strat_breaks = c(0, 0.03, 0.08, 0.15))$trend
    gmeth <- bind_rows(gmeth_low, gmeth_low_high) %>%
        select(track, breaks, cg_num, meth) %>%
        organize_samp_id()
    return(gmeth)
}
