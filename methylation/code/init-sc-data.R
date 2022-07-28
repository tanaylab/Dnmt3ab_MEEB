sc_md <- fread(here("data/sc_meth_plates.tsv")) %>% as_tibble()

load_cgdb <- function() {
    library(sc5mc)
    message("Loading cgdb...")
    db <<- cgdb_load(here("data/cgdb"))

    db_f <<- db %>%
        filter_cells(CHH <= 0.02) %>%
        left_join_cells(sc_md)

    invisible(db)
}

fill_sort_column <- function(db) {
    db %>%
        mutate_cells(
            sort = ifelse(
                sort == "index",
                case_when(
                    gate == "P11" ~ "CXCR4+EPCAM+",
                    gate == "P10" ~ "CXCR4+EPCAM-",
                    gate == "P12" ~ "CXCR4-EPCAM+"
                ),
                sort
            )
        )
}

add_gating_data <- function(db, gates) {
    plates <- unique(db@cells$plate)
    gate_data <- map_dfr(plates, function(plate) {
        if (tgcg::has_index_sort(plate)) {
            d <- tgcg:::get_plate_index_sort_csv(plate, gating = TRUE)
            d <- d %>%
                gather("gate", "val", -plate_pos) %>%
                filter(gate %in% gates, val) %>%
                select(plate_pos, gate) %>%
                mutate(plate = plate)
        } else {
            return(NULL)
        }
    })

    db <- db %>% left_join_cells(gate_data %>% distinct(), by = c("plate", "plate_pos"))
    return(db)
}

load_plpdb <- function() {
    library(sc5mc)
    message("Loading plpdb...")
    db_plp <<- cgdb_load(here("data/plpdb"))

    db_plp_f <<- db_plp %>%
        filter_cells(total_cov >= 5e4)

    invisible(db_plp)
}
