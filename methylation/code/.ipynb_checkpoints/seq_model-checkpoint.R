get_seq_df <- function(intervals, flank_bp){    
    intervals <- as_tibble(intervals)
    seq_df <- intervals %>% 
        distinct(chrom, start, end) %>%        
        # Remove positions outside of chromosome boundries
        left_join(gintervals.all() %>% select(chrom, start_chrom = start, end_chrom = end), by="chrom") %>% 
        filter(start - flank_bp - 1 >= start_chrom, end + flank_bp + 2 <= end_chrom) %>%
        select(-start_chrom, -end_chrom) %>%
        # Extract sequences
        mutate(s = gseq.extract(mutate(., start = start - flank_bp - 1, end = end + flank_bp + 2))) %>%
        mutate(nuc = strsplit(toupper(s), split="")) %>% 
#         unnest(c(nuc)) %>% 
        unnest_legacy(nuc) %>% 
        select(-s) %>% 
        group_by(chrom, start, end) %>% 
        mutate(pos = 1:n())    

    # Remove CpG positions
    seq_df <- seq_df %>% filter(!(pos %in% c(flank_bp + 2, flank_bp + 3) ))
    
    # Get dinucleotides (note that we join together position before and after the CpG)
    seq_df <- seq_df %>% mutate(dnuc = ifelse(is.na(lead(nuc)), NA, paste0(nuc, lead(nuc)))) 
    
    # Remove last position (we extract it in order to get the last dinuc)
    seq_df <- seq_df %>% filter(pos != max(pos))
    
    # Remove first position (we extract it in order to get the first dinuc)
    seq_df <- seq_df %>% filter(pos != min(pos))
    
    seq_df <- ungroup(seq_df)
    
    return(seq_df)
}

seq_df_to_wide <- function(seq_df, flank_bp, dinuc=TRUE){
    if (dinuc){
        seq_df_wide <- seq_df %>% select(-nuc) %>% mutate(n = 1) %>% pivot_wider(names_from = c(pos, dnuc), values_from=n, values_fill=list(n = 0)) %>% select(-contains("N", ignore.case=FALSE))        
        stopifnot(ncol(seq_df_wide) - 3 == 16 * flank_bp * 2)    
    } else {
        seq_df_wide <- seq_df %>% select(-dnuc) %>% mutate(n = 1) %>% pivot_wider(names_from = c(pos, nuc), values_from=n, values_fill=list(n = 0)) %>% select(-contains("N", ignore.case=FALSE))
        stopifnot(ncol(seq_df_wide) - 3 == 4 * flank_bp * 2)    
    }
    
    return(seq_df_wide)
}

gen_seq_model <- function(seq_df_wide, raw_mat, column){
    column <- enquo(column)
    mat <- seq_df_wide %>% inner_join(raw_mat %>% filter(!is.na(!!column)) %>% distinct(chrom, start, end, !!column), by=c("chrom", "start", "end"))
    feats <- mat %>% select(-(chrom:end), -!!column) %>%  as.matrix()
    feats <- feats[, sort(colnames(feats))] 
    y <- mat %>% pull(!! column)

    fit_cv <- glmnet::cv.glmnet(feats, y)
    pred <- predict(fit_cv, feats, s="lambda.min")
    return(list(mat = mat, feats = feats, y = y, fit_cv = fit_cv, pred = pred))
}

compute_intervals_ab_scores <- function(intervals, flank_bp, model_a, model_b, model_ab){    
    intervals <- intervals %>% select(chrom, start, end)

    message("Extracting sequence...")
    seq_df <- get_seq_df(intervals, flank_bp=flank_bp)

    message("Preparing features...")
    seq_df_wide <- seq_df_to_wide(seq_df, flank_bp=flank_bp, dinuc=TRUE)
    feats <- seq_df_wide %>% select(-(chrom:end)) %>%  as.matrix()
    feats <- feats[, colnames(model_ab$feats)]
    
    message("Predicting...")
    ab_score <- predict(model_ab$fit_cv, feats, s="lambda.min")
    a_score <- predict(model_a$fit_cv, feats, s="lambda.min")
    b_score <- predict(model_b$fit_cv, feats, s="lambda.min")
    
    res <- seq_df_wide %>% select(chrom, start, end) %>% mutate(ab_score = ab_score[, 1], a_score = a_score[, 1], b_score = b_score[, 1])
    return(res)    
}

compute_intervals_ab_scores_parallel <- function(intervals, flank_bp, model_a, model_b, model_ab){
    n_buckets <- round(nrow(intervals) / 1e5)
    intervals <- intervals %>% mutate(bucket = ntile(n = n_buckets))
    message(glue("# of buckets: {n_buckets}"))
    if (n_buckets > 1){        
        res <- plyr::ddply(intervals, "bucket", compute_intervals_ab_scores, flank_bp, model_a, model_b, model_ab, .parallel=TRUE) %>% 
            select(-bucket) %>% 
            as_tibble()
        return(res)
    }  

    return(compute_intervals_ab_scores(intervals, flank_bp, model_a, model_b, model_ab))
}
