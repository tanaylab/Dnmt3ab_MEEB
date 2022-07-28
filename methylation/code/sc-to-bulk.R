get_sc_germ_layer_metadata <- function() {
    metadata <- db_f@cells %>% select(cell_id, day, line, gate, sort, experiment)
    invivo_sort <- fread(here("data/cells_germ_layer_invivo.tsv")) %>% as_tibble()

    metadata <- bind_rows(
        metadata %>%
            filter(!grepl("^e", day)) %>%
            mutate(
                sort =
                    ifelse(
                        !is.na(gate),
                        case_when(
                            gate == "P10" ~ "CXCR4+EPCAM-",
                            gate == "P11" ~ "CXCR4+EPCAM+",
                            gate == "P12" ~ "CXCR4-EPCAM+"
                        ),
                        sort
                    )
            ) %>%
            mutate(sort = case_when(
                sort == "CXCR4+EPCAM+" ~ "endo",
                sort == "CXCR4+EPCAM-" ~ "meso",                
                sort == "CXCR4-EPCAM+" ~ "epi"
            )),
        metadata %>%
            filter(grepl("^e", day)) %>%
            select(-sort) %>%
            left_join(invivo_sort %>% rename(sort = germ_layer))
    ) %>%
        mutate(sort = ifelse(day == "d4", "epi", sort))

    metadata <- metadata %>%
        mutate(sort = ifelse(day != "e8.5", gsub("ecto", "epi/ecto", sort), sort))

    return(metadata)
}

sc_to_bulk_tracks <- function() {
    load_cgdb()
    db1 <- db_f %>%
        select_cells(plate, cell_id, cell_num) %>%
        left_join_cells(get_sc_germ_layer_metadata()) %>%
        filter_cells(!is.na(day), !is.na(line), !is.na(sort)) %>%
        filter_cells(sort != "other") %>%
        mutate_cells(
            track = glue("mEBDNMT.sc_bulk.{day}_{line}_{sort}"),
            track = gsub("/", "_", track),
            track = gsub("e7.5", "e7_5", track),
            track = gsub("\\?", "", track),
            track = gsub("e6.5", "e6_5", track),
            track = gsub("e8.5", "e8_5", track)
        )

    db1@cells <- db1@cells %>% add_count(track, name = "n_group")
    db1 <- db1 %>% filter_cells(n_group >= 20)

    gdir.create("mEBDNMT/sc_bulk", showWarnings = FALSE)

    df <- db1 %>%
        group_by_cells(track) %>%
        summarise()

    tracks <- unique(db1@cells$track)
    for (track in tracks) {
        print(track)
        if (!gtrack.exists(paste0(track, ".cov"))) {
            dff <- df %>%
                filter(track == !!track) %>%
                filter(cov > 0)
            gdir.create(gsub("\\.", "/", track), showWarnings = FALSE)
            gtrack.create_sparse(track = paste0(track, ".cov"), intervals = dff %>% select(chrom, start, end), values = dff$cov, description = "xxx")
            gtrack.create_sparse(track = paste0(track, ".meth"), intervals = dff %>% select(chrom, start, end), values = dff$meth, description = "xxx")
            gtrack.create_sparse(track = paste0(track, ".avg"), intervals = dff %>% select(chrom, start, end), values = dff$meth / dff$cov, description = "xxx")
        }
    }
}
