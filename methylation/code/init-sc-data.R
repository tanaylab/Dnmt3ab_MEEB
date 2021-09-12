sc_md <- fread(here("data/sc_meth_plates.tsv")) %>% as_tibble()
sc_md_dev <- fread(here("data/sc_meth_plates_dev.tsv")) %>% as_tibble()


init_ebdnmt_cgdb <- function(path = here("data/cgdb"), overwrite = FALSE){
    library(sc5mc)
    db <- tgcg::get_cgdb("mm9", reload_from_disk = TRUE)
    
    gvtrack.create("d_exon", "intervs.global.exons", "distance")
    gvtrack.create("d_tss", "intervs.global.tss", "distance")
    gvtrack.create("tor", "Encode.esd3.replichip.rep2", "avg")
    gvtrack.create("ab_score", "DNMT.ab_score_xgb_plus", "avg")
    gvtrack.create("a_score", "DNMT.a_score_xgb_plus", "avg")
    gvtrack.create("b_score", "DNMT.b_score_xgb_plus", "avg")    
    gvtrack.iterator("tor", sshift=-15000, eshift=15000)
    annot_tracks <- c("d_exon", "d_tss", "tor", "ab_score", "a_score", "b_score")
    db <- db %>% g_mutate_cpgs(annot_tracks)

    db <- db %>% filter_cells(plate %in% sc_md$plate)
    
    db <- db %>% filter_cells(plate != 'PZM00015')
    db <- db %>% mutate_cells(empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)

    db <- db %>% mutate_cells(passage = as.numeric(gsub('ES_line1_passage', '', cell_source)),
                     treatment = ifelse(treatment == '2i', '0d_2i', treatment), 
                     type = ifelse(cell_type == 'SemiMEF', '2d_semiMEF', treatment),
                     type = ifelse(type == '2d', paste0(type, '_passage', passage), type),
                     empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)
    bad_cpgs <- db %>% summarise_cells() %>% filter(cov > nrow(db@cells)) %>% pull(id)
    db <- db %>% filter_cpgs(!(id %in% bad_cpgs))
    
    db <- add_cgdb_sumamry_stats(db)

    db <- db %>% add_gating_data(gates = c("P10", "P11", "P12"))

    cgdb_save(db, path, force = overwrite)
    # filter_db(db)

    freemem(db)    
}

update_ebdnmt_cgdb <- function(path = here("data/cgdb")){
    library(sc5mc)
    genome_db <- tgcg::get_cgdb("mm9", reload_from_disk = TRUE)
    cur_db <- cgdb_load(path)

    cur_cells <- cur_db@cells
    all_cells <- genome_db@cells

    diff_plates <- setdiff(unique(all_cells$plate), unique(cur_cells$plate))
    ignore_plates <- c("eb5mc", paste0("PZM", str_pad(string = 1:80, pad = "0", width = 5) ) )    

    diff_db <- genome_db %>% filter_cells(plate %in% diff_plates, !(plate %in% ignore_plates))    
    diff_db <- diff_db %>% filter_cells(plate != 'PZM00015', grepl("^PZM", plate))
    diff_db <- diff_db %>% mutate_cells(empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)

    diff_db <- diff_db %>% 
        mutate_cells(passage = as.numeric(gsub('ES_line1_passage', '', cell_source)),
                     treatment = ifelse(treatment == '2i', '0d_2i', treatment), 
                     type = ifelse(cell_type == 'SemiMEF', '2d_semiMEF', treatment),
                     type = ifelse(type == '2d', paste0(type, '_passage', passage), type),
                     empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)    

    diff_db@cpgs <- cur_db@cpgs

    diff_plates <- unique(diff_db@cells$plate)    
    if (length(diff_plates) == 0){
        message("Nothing to update. In order to create the cgdb from scratch please run init_ebdnmt_cgdb")
        return()
    }    
    print(diff_plates)

    diff_db <- add_cgdb_sumamry_stats(diff_db)

    diff_db <- diff_db %>% add_gating_data(gates = c("P10", "P11", "P12"))

    diff_plates <- unique(diff_db@cells$plate)

    cur_db@cells <- rbind(cur_db@cells, diff_db@cells) %>% arrange(plate, cell_id)

    fwrite(cur_db@cells, paste0(path, "/cells.csv"))

    freemem(cur_db)
}

add_cgdb_sumamry_stats <- function(db, min_cov = 300){
    for (col in c("total_cov", "total_meth", "glob_avg", "low_cgc_m", "mid_cgc_m", "high_cgc_m")){
        if (col %in% colnames(db@cells)){
            db@cells[[col]] <- NULL
        }
    }    
    
    stats <- db %>% summarise_cpgs() %>% rename(total_cov = cov, total_meth = meth) %>% mutate(glob_avg = total_meth / total_cov)
    db <- db %>% left_join_cells(stats)
    
    cgc_trend <- get_cgc_trend(db, min_cov = min_cov, breaks=c(0,0.03,0.08,0.15)) %>% select(cell_id, low_cgc_m=`(0,0.03]`, mid_cgc_m = `(0.03,0.08]`, high_cgc_m = `(0.08,0.15]`)
    db <- db %>% left_join_cells(cgc_trend)
    
    return(db)    
}

load_cgdb <- function(){    
    library(sc5mc)   
    message('Loading cgdb...')
    db <<- cgdb_load(here("data/cgdb"))
        
    db_f <<- db %>% 
        filter_cells(CHH <= 0.02) %>% 
        left_join_cells(sc_md)
    
    invisible(db)
}

load_cgdb_dev <- function(){    
    library(sc5mc)   
    message('Loading cgdb...')
    db <<- cgdb_load(here("data/cgdb"))
        
    db_f <<- db %>% 
        filter_cells(CHH <= 0.02) %>% 
        left_join_cells(sc_md_dev)
    
    invisible(db)
}


fill_sort_column <- function(db){
    db %>% 
        mutate_cells(
            sort = ifelse(
                sort == "index", 
                case_when(
                    gate == "P11" ~ "CXCR4+EPCAM+", 
                    gate == "P10" ~ "CXCR4+EPCAM-", 
                    gate == "P12" ~ "CXCR4-EPCAM+"), 
            sort)
        )
}

add_gating_data <- function(db, gates){
    plates <- unique(db@cells$plate)
    gate_data <- map_dfr(plates, function(plate){
        if (tgcg::has_index_sort(plate)){
            d <- tgcg:::get_plate_index_sort_csv(plate, gating=TRUE)
            d <- d %>% gather("gate", "val", -plate_pos) %>% filter(gate %in% gates, val) %>% select(plate_pos, gate) %>% mutate(plate = plate)    
        } else {
            return(NULL)
        }        
    })

    db <- db %>% left_join_cells(gate_data %>% distinct(), by=c("plate", "plate_pos"))
    return(db)    
}

init_ebdnmt_plpdb <- function(path = here("data/plpdb"), overwrite = FALSE){
    library(sc5mc)
    db <- tgcg::get_plpdb("mm9", reload_from_disk = TRUE)
    
    gvtrack.create("tor", "Encode.esd3.replichip.rep2", "avg")
    gvtrack.iterator("tor", sshift=-15000, eshift=15000)
    annot_tracks <- c("tor")
    db <- db %>% g_mutate_cpgs(annot_tracks)

    db <- db %>% filter_cells(plate %in% sc_md$plate)
    
    db <- db %>% filter_cells(plate != 'PZM00015')
    db <- db %>% mutate_cells(empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)

    db <- db %>% mutate_cells(passage = as.numeric(gsub('ES_line1_passage', '', cell_source)),
                     treatment = ifelse(treatment == '2i', '0d_2i', treatment), 
                     type = ifelse(cell_type == 'SemiMEF', '2d_semiMEF', treatment),
                     type = ifelse(type == '2d', paste0(type, '_passage', passage), type),
                     empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)
    bad_cpgs <- db %>% summarise_cells() %>% filter(cov > nrow(db@cells)) %>% pull(id)
    db <- db %>% filter_cpgs(!(id %in% bad_cpgs))

    stats <- db %>% summarise_cpgs() %>% select(cell_id, total_cov = cov)
    db <- db %>% left_join_cells(stats)
    
    cgdb_save(db, path, force = overwrite)    

    freemem(db)    
}

update_ebdnmt_plpdb <- function(path = here("data/plpdb")){
    library(sc5mc)
    genome_db <- tgcg::get_plpdb("mm9", reload_from_disk = TRUE)
    cur_db <- plpdb_load(path)

    cur_cells <- cur_db@cells
    all_cells <- genome_db@cells

    diff_plates <- setdiff(unique(all_cells$plate), unique(cur_cells$plate))
    ignore_plates <- c("eb5mc", paste0("PZM", str_pad(string = 1:80, pad = "0", width = 5) ) )    

    diff_db <- genome_db %>% filter_cells(plate %in% diff_plates, !(plate %in% ignore_plates))  
    diff_db <- diff_db %>% filter_cells(plate != 'PZM00015', grepl("^PZM", plate))
    diff_db <- diff_db %>% mutate_cells(empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)

    diff_db <- diff_db %>% 
        mutate_cells(passage = as.numeric(gsub('ES_line1_passage', '', cell_source)),
                     treatment = ifelse(treatment == '2i', '0d_2i', treatment), 
                     type = ifelse(cell_type == 'SemiMEF', '2d_semiMEF', treatment),
                     type = ifelse(type == '2d', paste0(type, '_passage', passage), type),
                     empty = ifelse(is.na(empty), FALSE, empty)) %>%
                 filter_cells(!empty)    

    diff_db@cpgs <- cur_db@cpgs

    
    diff_plates <- unique(diff_db@cells$plate)    
    if (length(diff_plates) == 0){
        message("Nothing to update. In order to create the cgdb from scratch please run init_ebdnmt_plpdb")
        return()
    }        
    print(diff_plates)

    stats <- diff_db %>% summarise_cpgs() %>% select(cell_id, total_cov = cov)

    diff_db <- diff_db %>% left_join_cells(stats)
   
    cur_db@cells <- rbind(cur_db@cells, diff_db@cells) %>% arrange(plate, cell_id)

    fwrite(cur_db@cells, paste0(path, "/cells.csv"))

    freemem(cur_db)
}

load_plpdb <- function(){     
    library(sc5mc) 
    message('Loading plpdb...')
    db_plp <<- cgdb_load(here("data/plpdb"))
        
    db_plp_f <<- db_plp %>% 
        filter_cells(total_cov >= 5e4)        
    
    invisible(db_plp)
}



# update_db <- function(path=here("data/cgdb"), update_plates = NULL){
#     db <- cgdb_load(path)
#     core_db_cells <- fread(glue(conf$cgdb_data_root, '/cells.csv')) %>% as_tibble()
#     gsetroot(conf$groot)

#     if (!is.null(update_plates)){
#         db@cells <- db@cells %>% filter(!(plate %in% update_plates))
#     }

#     cells_md <- core_db_cells %>% filter(!(cell_id %in% db@cells$cell_id)) %>% mutate(empty = ifelse(is.na(empty), FALSE, empty)) %>% filter(!empty)

#     stats <- get_plate_stats()
#     cells_md <- cells_md %>% left_join(stats %>% select(cell_id, cg_num, total_reads, mapped_reads, mapped_frac, uniq_reads, uniq_frac, CHH, CHG, CpG))

#     new_cells_md <- db@cells %>% bind_rows(cells_md) %>% distinct(cell_id, .keep_all = TRUE)
    
#     db <- add_cgdb_sumamry_stats(db)
    
#     fwrite(new_cells_md, glue(path, '/cells.csv'))  
# }

