get_seq_df <- function(intervals, flank_bp, strand = 1) {
    intervals <- as_tibble(intervals)
    seq_df <- intervals %>%
        distinct(chrom, start, end) %>%
        # Remove positions outside of chromosome boundries
        left_join(gintervals.all() %>% select(chrom, start_chrom = start, end_chrom = end), by = "chrom") %>%
        filter(start - flank_bp - 1 >= start_chrom, end + flank_bp + 2 <= end_chrom) %>%
        select(-start_chrom, -end_chrom)

    seq_df <- seq_df %>%
        mutate(strand = strand) %>%
        select(chrom, start, end, everything()) %>%
        # Extract sequences
        mutate(s = gseq.extract(mutate(., start = start - flank_bp - 1, end = end + flank_bp + 2))) %>%
        select(-strand) %>%
        mutate(nuc = strsplit(toupper(s), split = "")) %>%
        #         unnest(c(nuc)) %>%
        unnest_legacy(nuc) %>%
        select(-s) %>%
        group_by(chrom, start, end) %>%
        mutate(pos = 1:n())

    # Remove CpG positions
    seq_df <- seq_df %>% filter(!(pos %in% c(flank_bp + 2, flank_bp + 3)))

    # Get dinucleotides (note that we join together position before and after the CpG)
    seq_df <- seq_df %>% mutate(dnuc = ifelse(is.na(lead(nuc)), NA, paste0(nuc, lead(nuc))))

    # Remove last position (we extract it in order to get the last dinuc)
    seq_df <- seq_df %>% filter(pos != max(pos))

    # Remove first position (we extract it in order to get the first dinuc)
    seq_df <- seq_df %>% filter(pos != min(pos))

    seq_df <- ungroup(seq_df)

    return(seq_df)
}

get_kmers_mat <- function(intervals, flank_bp, kmer_len) {
    intervals <- as_tibble(intervals)
    cpgs <- intervals %>%
        distinct(chrom, start, end) %>%
        # Remove positions outside of chromosome boundries
        left_join(gintervals.all() %>% select(chrom, start_chrom = start, end_chrom = end), by = "chrom") %>%
        filter(start - flank_bp - 1 >= start_chrom, end + flank_bp + 2 <= end_chrom) %>%
        select(-start_chrom, -end_chrom)

    cpgs <- cpgs %>% mutate(s = toupper(gseq.extract(mutate(., start = start - flank_bp - 1, end = end + flank_bp + 2))))
    alphabet <- c("A", "C", "G", "T")
    kmers <- expand.grid(map(1:kmer_len, ~alphabet)) %>%
        unite("kmer", starts_with("Var"), sep = "") %>%
        pull(kmer)

    res <- t(sapply(cpgs$s, stringi::stri_count_fixed, pattern = kmers, overlap = TRUE, USE.NAMES = FALSE))
    colnames(res) <- kmers
    res <- cbind(cpgs %>% select(chrom, start, end), res) %>% intervs_to_mat()

    return(res)
}

seq_df_to_wide <- function(seq_df, flank_bp, dinuc = TRUE) {
    if (dinuc) {
        seq_df_wide <- seq_df %>%
            select(-nuc) %>%
            mutate(n = 1) %>%
            pivot_wider(names_from = c(pos, dnuc), values_from = n, values_fill = list(n = 0)) %>%
            select(-contains("N", ignore.case = FALSE))
        stopifnot(ncol(seq_df_wide) - 3 == 16 * flank_bp * 2)
    } else {
        seq_df_wide <- seq_df %>%
            select(-dnuc) %>%
            mutate(n = 1) %>%
            pivot_wider(names_from = c(pos, nuc), values_from = n, values_fill = list(n = 0)) %>%
            select(-contains("N", ignore.case = FALSE))
        stopifnot(ncol(seq_df_wide) - 3 == 4 * flank_bp * 2)
    }

    return(seq_df_wide)
}

gen_seq_model <- function(seq_df_wide, raw_mat, column, ...) {
    column <- enquo(column)
    mat <- seq_df_wide %>% inner_join(raw_mat %>% filter(!is.na(!!column)) %>% distinct(chrom, start, end, !!column), by = c("chrom", "start", "end"))
    feats <- mat %>%
        select(-(chrom:end), -!!column) %>%
        as.matrix()
    feats <- feats[, sort(colnames(feats))]
    y <- mat %>% pull(!!column)

    fit_cv <- glmnet::cv.glmnet(feats, y, ...)
    pred <- predict(fit_cv, feats, s = "lambda.min")
    return(list(mat = mat, feats = feats, y = y, fit_cv = fit_cv, pred = pred))
}

hypertune_xgb <- function(seq_df_wide, raw_mat, column) {
    column <- enquo(column)
    mat <- seq_df_wide %>% inner_join(raw_mat %>% filter(!is.na(!!column)) %>% distinct(chrom, start, end, !!column), by = c("chrom", "start", "end"))
    feats <- mat %>%
        select(-(chrom:end), -!!column) %>%
        as.matrix()
    feats <- feats[, sort(colnames(feats))]
    y <- mat %>% pull(!!column)

    library(xgboost)
    slice <- dplyr::slice
    dtrain <- xgb.DMatrix(feats, label = y)

    nrounds <- 500
    # learning rate + nrounds

    tune_grid <- expand.grid(
        nrounds = seq(50, nrounds, by = 50),
        eta = c(0.01, 0.015, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
        gamma = 0,
        max_depth = c(2, 4, 6),
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
    )

    tune_control <- caret::trainControl(
        method = "cv",
        number = 3,        
        verboseIter = FALSE,
        allowParallel = FALSE
    )

    xgb_tune <- caret::train(
        x = dtrain,
        y = y,
        trControl = tune_control,
        tuneGrid = tune_grid,
        method = "xgbTree",
        verbose = TRUE,
        metric = "RMSE"
    )

    message(paste("hp1 RMSE:", (xgb_tune$results %>% inner_join(xgb_tune$bestTune))$RMSE))

    xgb_tune$bestTune



    # max depth min child weight
    max_depth <- (xgb_tune$bestTune$max_depth - 1):(xgb_tune$bestTune$max_depth + 1)

    if (xgb_tune$bestTune$max_depth == 2) {
        max_depth <- 2:4
    }

    tune_grid2 <- expand.grid(
        nrounds = seq(50, nrounds, by = 50),
        eta = xgb_tune$bestTune$eta,
        max_depth = max_depth,
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = c(1, 2, 3),
        subsample = 1
    )

    xgb_tune2 <- caret::train(
        x = dtrain,
        y = y,
        trControl = tune_control,
        tuneGrid = tune_grid2,
        method = "xgbTree",
        verbose = TRUE,
        metric = "RMSE"
    )


    message(paste("hp2 RMSE:", (xgb_tune2$results %>% inner_join(xgb_tune2$bestTune))$RMSE))
    xgb_tune2$bestTune

    # column and row sampling

    tune_grid3 <- expand.grid(
        nrounds = seq(50, nrounds, by = 50),
        eta = xgb_tune$bestTune$eta,
        max_depth = xgb_tune2$bestTune$max_depth,
        gamma = 0,
        colsample_bytree = c(0.4, 0.6, 0.8, 1),
        min_child_weight = xgb_tune2$bestTune$min_child_weight,
        subsample = c(0.5, 0.75, 1)
    )

    xgb_tune3 <- caret::train(
        x = dtrain,
        y = y,
        trControl = tune_control,
        tuneGrid = tune_grid3,
        method = "xgbTree",
        verbose = TRUE,
        metric = "RMSE"
    )

    message(paste("hp3 RMSE:", (xgb_tune3$results %>% inner_join(xgb_tune3$bestTune))$RMSE))
    xgb_tune3$bestTune

    # choosing gamma

    tune_grid4 <- expand.grid(
        nrounds = seq(50, nrounds, by = 50),
        eta = xgb_tune$bestTune$eta,
        max_depth = xgb_tune2$bestTune$max_depth,
        gamma = c(0, 0.1, 0.3, 0.7, 0.9, 1),
        colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
        min_child_weight = xgb_tune2$bestTune$min_child_weight,
        subsample = xgb_tune3$bestTune$subsample
    )


    xgb_tune4 <- caret::train(
        x = dtrain,
        y = y,
        trControl = tune_control,
        tuneGrid = tune_grid4,
        method = "xgbTree",
        verbose = TRUE,
        metric = "RMSE"
    )

    message(paste("hp4 RMSE:", (xgb_tune4$results %>% inner_join(xgb_tune4$bestTune))$RMSE))
    xgb_tune4$bestTune

    # reducing the learning rate

    max_rounds <- 2000    

    tune_grid5 <- expand.grid(
        nrounds = seq(50, max_rounds, by = 50),
        eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
        max_depth = xgb_tune2$bestTune$max_depth,
        gamma = xgb_tune4$bestTune$gamma,
        colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
        min_child_weight = xgb_tune2$bestTune$min_child_weight,
        subsample = xgb_tune3$bestTune$subsample
    )



    xgb_tune5 <- caret::train(
        x = dtrain,
        y = y,
        trControl = tune_control,
        tuneGrid = tune_grid5,
        method = "xgbTree",
        verbose = TRUE,
        metric = "RMSE"
    )



    message(paste("hp5 RMSE:", (xgb_tune5$results %>% inner_join(xgb_tune5$bestTune))$RMSE))
    xgb_tune5$bestTune



    params <- list(
        booster = "gbtree",
        objective = "reg:squarederror",
        subsample = xgb_tune3$bestTune$subsample,
        max_depth = xgb_tune2$bestTune$max_depth,
        colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
        gamma = xgb_tune4$bestTune$gamma,
        eta = xgb_tune5$bestTune$eta,
        eval_metric = "rmse",
        min_child_weight = xgb_tune2$bestTune$min_child_weight
    )

    return(list(params = params, nrounds = xgb_tune5$bestTune$nrounds))
}

gen_seq_model_xgboost <- function(seq_df_wide, raw_mat, column, params=NULL, ...) {
    column <- enquo(column)
    mat <- seq_df_wide %>% inner_join(raw_mat %>% filter(!is.na(!!column)) %>% distinct(chrom, start, end, !!column), by = c("chrom", "start", "end"))
    feats <- mat %>%
        select(-(chrom:end), -!!column) %>%
        as.matrix()
    feats <- feats[, sort(colnames(feats))]
    y <- mat %>% pull(!!column)

    library(xgboost)
    slice <- dplyr::slice
    dtrain <- xgb.DMatrix(feats, label = y)

    
    if (is.null(params)){
        params <- list(params = list(booster = "gbtree", objective = "reg:squarederror", eta = 0.3, gamma = 0, max_depth = 6, min_child_weight = 1, subsample = 1, colsample_bytree = 1), nrounds = 200)
    }    
    
    xgbcv <- xgb.cv(params = params$params, data = dtrain, nrounds = params$nrounds, nfold = 10, prediction = TRUE)

    xgbtrain <- xgb.train(params = params$params, data = dtrain, nrounds = params$nrounds)

    # pred <- predict(xgbtrain, dtrain)

    return(list(mat = mat, feats = feats, y = y, fit_cv = xgbtrain, xgbcv = xgbcv, pred = xgbcv$pred))
}

compute_intervals_ab_scores <- function(intervals, flank_bp, model_a, model_b, model_ab, model_ab_xgb, model_a_xgb, model_b_xgb) {
    gc()
    intervals <- intervals %>% select(chrom, start, end)

    message("Extracting sequence...")
    seq_df_plus <- get_seq_df(intervals, flank_bp = flank_bp, strand = 1)
    seq_df_minus <- get_seq_df(intervals, flank_bp = flank_bp, strand = -1)

    message("Preparing features...")
    seq_df_wide_plus <- seq_df_to_wide(seq_df_plus, flank_bp = flank_bp, dinuc = TRUE)
    seq_df_wide_minus <- seq_df_to_wide(seq_df_minus, flank_bp = flank_bp, dinuc = TRUE)
    feats_plus <- seq_df_wide_plus %>%
        select(-(chrom:end)) %>%
        as.matrix()
    feats_plus <- feats_plus[, colnames(model_ab$feats)]

    feats_minus <- seq_df_wide_minus %>%
        select(-(chrom:end)) %>%
        as.matrix()
    feats_minus <- feats_minus[, colnames(model_ab$feats)]

    message("Predicting...")
    ab_score_plus <- predict(model_ab$fit_cv, feats_plus, s="lambda.min")
    ab_score_xgb_plus <- predict(model_ab_xgb$fit_cv, feats_plus, s = "lambda.min")
    a_score_plus <- predict(model_a$fit_cv, feats_plus, s="lambda.min")
    a_score_xgb_plus <- predict(model_a_xgb$fit_cv, feats_plus, s = "lambda.min")
    b_score_plus <- predict(model_b$fit_cv, feats_plus, s="lambda.min")
    b_score_xgb_plus <- predict(model_b_xgb$fit_cv, feats_plus, s = "lambda.min")

    ab_score_minus <- predict(model_ab$fit_cv, feats_minus, s="lambda.min")
    ab_score_xgb_minus <- predict(model_ab_xgb$fit_cv, feats_minus, s = "lambda.min")
    a_score_minus <- predict(model_a$fit_cv, feats_minus, s="lambda.min")
    a_score_xgb_minus <- predict(model_a_xgb$fit_cv, feats_minus, s = "lambda.min")
    b_score_minus <- predict(model_b$fit_cv, feats_minus, s="lambda.min")
    b_score_xgb_minus <- predict(model_b_xgb$fit_cv, feats_minus, s = "lambda.min")

    message("Merging predictions...")
    # Note that a model for methylation of A-/- predicts the activity of DNMT3B, and B-/- predict the activity of DNMT3A. 
    res <- seq_df_wide_plus %>%
        select(chrom, start, end) %>%
        mutate(
            ab_score_glm_plus = ab_score_plus[, 1],
            ab_score_xgb_plus = ab_score_xgb_plus,
            b_score_glm_plus = a_score_plus[, 1],
            b_score_xgb_plus = a_score_xgb_plus,
            a_score_glm_plus = b_score_plus[, 1],
            a_score_xgb_plus = b_score_xgb_plus,
            ab_score_glm_minus = ab_score_minus[, 1],
            ab_score_xgb_minus = ab_score_xgb_minus,
            b_score_glm_minus = a_score_minus[, 1],
            b_score_xgb_minus = a_score_xgb_minus,
            a_score_glm_minus = b_score_minus[, 1],
            a_score_xgb_minus = b_score_xgb_minus,
        )

    return(res)
}

compute_intervals_ab_scores_parallel <- function(intervals, flank_bp, model_a, model_b, model_ab, model_ab_xgb, model_a_xgb, model_b_xgb) {
    library(glmnet)
    doMC::registerDoMC(10)
    n_buckets <- round(nrow(intervals) / 1e5)
    intervals <- intervals %>% mutate(bucket = ntile(n = n_buckets))
    message(glue("# of buckets: {n_buckets}"))
    if (n_buckets > 1) {
        res <- plyr::ddply(intervals, "bucket", compute_intervals_ab_scores, flank_bp, model_a, model_b, model_ab, model_ab_xgb, model_a_xgb, model_b_xgb, .parallel = TRUE) %>%
            select(-bucket) %>%
            as_tibble()
        return(res)
    }

    return(compute_intervals_ab_scores(intervals, flank_bp, model_a, model_b, model_ab))
}


get_coef_df <- function(model, relevel = TRUE) {
    library(glmnet)
    tmp_coeffs <- coef(model$fit_cv, s = "lambda.min")
    coef_df <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x) %>%
        dplyr::slice(-1) %>%
        as_tibble() %>%
        separate(name, c("pos", "dinuc"), sep = "_")

    # We remove last position from plotting
    coef_df <- coef_df %>% filter(pos != 13)

    if (relevel) {
        new_levels <- as.character(c(2:6, 9:12))
        names(new_levels) <- c(-4:4)
        coef_df <- coef_df %>%
            mutate(pos = factor(pos, levels = new_levels), pos = fct_recode(pos, !!!new_levels))
    } else {
        coef_df <- coef_df %>%
            mutate(pos = factor(as.numeric(pos)))
    }

    # names(new_levels) <- c(-4:-1, 0, 1:4)

    dinuc_levels <- c(
        "AT",
        "TA",
        "CG",
        "GC",
        "AA",
        "TT",
        "AC",
        "GT",
        "AG",
        "CT",
        "CA",
        "TG",
        "CC",
        "GG",
        "GA",
        "TC"
    )

    stopifnot(length(dinuc_levels) == 16)

    coef_df <- coef_df %>%
        mutate(dinuc = factor(dinuc, levels = rev(dinuc_levels))) %>%
        tidyr::complete(pos, dinuc, fill = list(coefficient = NA))

    return(coef_df)
}

sym_coefs <- function(model) {
    library(glmnet)
    tmp_coeffs <- coef(model$fit_cv, s = "lambda.min")
    coef_df <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x) %>%
        slice(-1) %>%
        as_tibble() %>%
        separate(name, c("pos", "dinuc"), sep = "_") %>%
        mutate(pos = as.numeric(pos))

    rc <- function(s) {
        return(stringi::stri_reverse(chartr("acgtACGT", "tgcaTGCA", s)))
    }

    dinuc_levels <- c(
        "AT",
        "TA",
        "CG",
        "GC",
        "AA",
        "TT",
        "AC",
        "GT",
        "AG",
        "CT",
        "CA",
        "TG",
        "CC",
        "GG",
        "GA",
        "TC"
    )

    coef_df <- coef_df %>%
        mutate(dinuc = factor(dinuc, levels = rev(dinuc_levels))) %>%
        tidyr::complete(pos, dinuc, fill = list(coefficient = NA))

    # coef_df <- coef_df %>% filter(pos != 13)

    new_levels <- as.character(c(2:6, 9:12))
    names(new_levels) <- c(-4:4)

    coef_df <- coef_df %>%
        mutate(pos = factor(pos, levels = new_levels), pos = fct_recode(pos, !!!new_levels))

    sym_coef_df <- coef_df %>%
        mutate(sym_pos = -as.numeric(as.character(pos)), sym_dinuc = rc(dinuc)) %>%
        mutate(sym_pos = factor(sym_pos, levels = levels(pos))) %>%
        left_join(coef_df %>% select(sym_pos = pos, sym_dinuc = dinuc, sym_coefficient = coefficient)) %>%
        replace_na(replace = list(coefficient = 0, sym_coefficient = 0)) %>%
        mutate(sum_coef = coefficient + sym_coefficient)
    return(sym_coef_df)
}

coef_df_to_matrix <- function(coef_df, model, intercept = NULL) {
    new_levels <- as.character((-4:4))
    names(new_levels) <- as.character(c(2:6, 9:12))

    coef_df <- coef_df %>%
        mutate(pos = factor(pos, levels = new_levels), pos = fct_recode(pos, !!!new_levels))


    orig_coef_mat <- coef(model$fit_cv, s = "lambda.min")

    coef_mat <- coef_df %>%
        filter(!is.na(pos), !is.na(dinuc)) %>%
        unite("rowname", pos, dinuc) %>%
        select(rowname, coefficient) %>%
        as.data.frame() %>%
        column_to_rownames("rowname") %>%
        as.matrix()

    coef_mat <- rbind(orig_coef_mat[1, ], coef_mat)
    rownames(coef_mat)[1] <- rownames(orig_coef_mat)[1]
    colnames(coef_mat) <- colnames(orig_coef_mat)

    if (!is.null(intercept)) {
        coef_mat[1, ] <- intercept
    }

    mat_rownames <- rownames(orig_coef_mat)[!grepl("13_", rownames(orig_coef_mat))]

    coef_mat <- as.matrix(coef_mat[mat_rownames, ])
    coef_mat[is.na(coef_mat)] <- 0

    return(coef_mat)
}

compute_interval_model_mat_score <- function(intervals, model, mat_ab, mat_a, mat_b, flank_bp = 5) {
    library(glmnet)
    intervals <- intervals %>% select(chrom, start, end)

    message("Extracting sequence...")
    seq_df_plus <- get_seq_df(intervals, flank_bp = flank_bp, strand = 1)
    seq_df_minus <- get_seq_df(intervals, flank_bp = flank_bp, strand = -1)

    message("Preparing features...")
    seq_df_wide_plus <- seq_df_to_wide(seq_df_plus, flank_bp = flank_bp, dinuc = TRUE)
    seq_df_wide_minus <- seq_df_to_wide(seq_df_minus, flank_bp = flank_bp, dinuc = TRUE)
    feats_plus <- seq_df_wide_plus %>%
        select(-(chrom:end)) %>%
        as.matrix()
    feats_plus <- feats_plus[, colnames(model$feats)]

    feats_minus <- seq_df_wide_minus %>%
        select(-(chrom:end)) %>%
        as.matrix()
    feats_minus <- feats_minus[, colnames(model$feats)]

    orig_mat <- coef(model$fit_cv, s = "lambda.min")
    orig_mat_trans <- sym_coefs(model) %>%
        select(pos, dinuc, coefficient) %>%
        coef_df_to_matrix(model_ab)

    message("Predicting...")
    res <- seq_df_wide_plus %>%
        select(chrom, start, end) %>%
        mutate(
            score_plus = predict_mat(mat_ab, feats_plus),
            score_minus = predict_mat(mat_ab, feats_minus),
            score_a_plus = predict_mat(mat_a, feats_plus),
            score_a_minus = predict_mat(mat_a, feats_minus),
            score_b_plus = predict_mat(mat_b, feats_plus),
            score_b_minus = predict_mat(mat_b, feats_minus),
            score_model_plus = predict_mat(orig_mat_trans, feats_plus),
            score_model_minus = predict_mat(orig_mat_trans, feats_minus),
            score_orig_plus = predict(model$fit_cv, feats_plus, s = "lambda.min"),
            score_orig_minus = predict(model$fit_cv, feats_minus, s = "lambda.min")
        )

    # # verification
    # score <- predict(model$fit_cv, feats, s="lambda.min")
    # res <- seq_df_wide %>%
    #     select(chrom, start, end) %>%
    #     mutate(
    #         score_pred = score[, 1],
    #         score_mat = predict_mat(mat, feats)
    #     )


    return(res)
}

predict_mat <- function(w, A) {
    A <- A[, rownames(w)[-1]]
    y <- t(w[-1, ]) %*% t(A) + w[1]
    y <- as.vector(y)
    return(y)
}


compute_interval_model_mat_score_parallel <- function(intervals, model_ab, mat_ab, mat_a, mat_b) {
    library(glmnet)
    n_buckets <- round(nrow(intervals) / 1e5)
    intervals <- intervals %>% mutate(bucket = ntile(n = n_buckets))
    message(glue("# of buckets: {n_buckets}"))
    if (n_buckets > 1) {
        res <- plyr::ddply(intervals, "bucket", compute_interval_model_mat_score, model_ab, mat_ab, mat_a, mat_b, .parallel = TRUE) %>%
            select(-bucket) %>%
            as_tibble()
        return(res)
    } else {
        res <- compute_interval_model_mat_score(intervals, model_ab, mat_ab, mat_a, mat_b)
    }


    # gtrack.create_sparse(track = "DNMT.baubec_ab_score", intervals = res %>% select(chrom:end), values=res$score, description="")
    # gtrack.create_sparse(track = "DNMT.baubec_a_score_strand", intervals = res %>% select(chrom:end), values=res$score_a, description="")
    # gtrack.create_sparse(track = "DNMT.baubec_b_score_strand", intervals = res %>% select(chrom:end), values=res$score_b, description="")
    # browser()

    return(res)
}


plot_model_scatter <- function(model, x_lab="", y_lab="", xlim=NULL, ylim=NULL, bandwidth=0.08, point_size=0.001){
    tibble(pred = model$pred, y = model$y) %>% 
        mutate(col = densCols(., bandwidth=0.06,colramp=colorRampPalette(c("white","lightblue", "blue", "darkblue", "yellow", "gold","orange","red", "darkred" )))) %>% 
        ggplot(aes(x=pred, y=y, col=col)) + 
            geom_point(shape=19, size=point_size) + 
            scale_color_identity() + 
            coord_cartesian(xlim = xlim, ylim = ylim) +                 
            xlab(x_lab) + 
            ylab(y_lab) +         
            theme(aspect.ratio=1, panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
            labs(subtitle = glue("R^2 = {cor}", cor = round(cor(model$pred, model$y)^2, digits=2))) + 
            theme(plot.subtitle = ggtext::element_markdown())
}

plot_model_scatter_legend <- function(model, colors=c("white","lightblue", "blue", "darkblue", "yellow", "gold","orange","red", "darkred" )){
    dc <- densCols(tibble(pred = model$pred, y = model$y), bandwidth=0.06,colramp=colorRampPalette(colors))
    dd <- grDevices:::.smoothScatterCalcDensity(data.frame(pred = model$pred, y = model$y), bandwidth = 0.06, nbin=128)
    dens <- as.numeric(dd$fhat)
    dens <- dens[dens>0]

    n_colors <- 1000
    legend <- data.frame(density=seq(min(dens), max(dens), len=n_colors), color=colorRampPalette(colors)(n_colors))
    plot(NA, xlim=c(0,n_colors), ylim=c(0,n_colors+1), type="n", ann=FALSE, axes=FALSE)
    rect(0, 1:n_colors, 50, 2:(n_colors+1), border=NA, col=legend$col)
    # text(2, (1:n_colors)+0.5, signif(legend$density, 2))  
    
    
}