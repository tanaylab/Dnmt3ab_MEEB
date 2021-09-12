gintervals.neighbors1 <- gpatterns::gintervals.neighbors1
fread <- function(...) as_tibble(tgutil::fread(...))

gintervals.centers <- function(inv) {
    inv %>%
        mutate(
            center = floor((start + end) / 2), start = center,
            end = center + 1
        ) %>%
        select(-center)
}


psum <- function(..., na.rm=FALSE) { 
    dat <- do.call(cbind, list(...))
    res <- rowSums(dat, na.rm=na.rm) 
    idx_na <- !rowSums(!is.na(dat))
    res[idx_na] <- NA
    return(res)
}

pmean <- function(..., na.rm=FALSE) { 
    dat <- do.call(cbind, list(...))
    res <- rowMeans(dat, na.rm=na.rm) 
    idx_na <- !rowMeans(!is.na(dat))
    res[idx_na] <- NA
    return(res)
}
