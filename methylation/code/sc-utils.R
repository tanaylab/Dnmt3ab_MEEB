
calc_cgc_ab_score_sc <- function(db_f, score_var){
    k <- 20

    score_var <- enquo(score_var)

    message("Classifying CpGs...")
    db1 <- db_f %>%
        mutate_cells(bucket = ntile(n = k)) %>%
        mutate_cpgs(
            score = cut(!! score_var, quantile(!! score_var, c(0,0.1,0.9,1))),
            cg_cont = cut(cg500, c(0,0.02,0.08,0.2)),
        tor_grp = case_when(tor <= 0 ~ "late", tor >= 0 ~ "early")) %>%
        filter_cpgs(chrom != "chrX", chrom != "chrY", !is.na(tor_grp), !is.na(cg_cont), !is.na(score))


    message("Extracting methylation...")
    res <- plyr::ldply(1:k, function(i) {
        print(i)
        freemem(db1)
        db1 %>%
            filter_cells(bucket == i) %>%
            group_by_cpgs(score, cg_cont, tor_grp) %>%
            summarise() %>%
            mutate(avg = meth / cov) %>%
            pivot_wider(c(cell_id, cg_cont, score), names_from=tor_grp, values_from=cov:avg)
        }, .parallel = TRUE)    

    return(res)
}


