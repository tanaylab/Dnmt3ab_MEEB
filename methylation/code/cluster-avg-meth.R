filter_meth_mat <- function(mat, sd_thresh, k_trend) {
    d <- tibble(m = rowMeans(mat, na.rm = TRUE), sd = matrixStats::rowSds(mat, na.rm = TRUE)) %>%
        mutate(i = 1:n()) %>%
        arrange(m) %>%
        mutate(
            trend =
                zoo::rollmedian(sd, k_trend, fill = c(
                    median(sd[1:k_trend], na.rm = TRUE),
                    NA,
                    median(sd[(n() - k_trend):n()], na.rm = TRUE)
                ))
        ) %>%
        mutate(chosen = sd >= trend + sd_thresh)
}
