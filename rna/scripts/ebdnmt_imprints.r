
test_imprints <- function(d) {
    imprint <- fread("config/imprinted.txt")
    rownames(imprint) <- imprint$gname

    # build skeleton - best correlation of each imprinted gene

    imp_cor <- tgs_cor(t(d$legc_iv[imprint$gname, ]))
    diag(imp_cor) <- 0
    best_mate <- rownames(imp_cor)[apply(imp_cor, 1, which.max)]
    names(best_mate) <- imprint$gname

    lines <- fread("config/lines.txt")

    eps <- 2e-5
    imp_type <- list()
    for (mut in c("wt", "3a", "3b", "dko")) {
        mc <- d[[sprintf("mc_%s", mut)]]
        ls <- lines[lines$type == mut, "line"]
        mut_l <- comp_egc_on_lines(d,
            type = mut,
            line1 = ls[1], line2 = ls[2]
        )
        egc_type1 <- t(tgs_matrix_tapply(mut_l$egc1, mc@colors[colnames(mut_l$egc1)], sum))
        rownames(egc_type1) <- rownames(mut_l$egc1)
        egc_type2 <- t(tgs_matrix_tapply(mut_l$egc2, mc@colors[colnames(mut_l$egc2)], sum))
        rownames(egc_type2) <- rownames(mut_l$egc2)
        tot_type1 <- colSums(egc_type1)
        tot_type2 <- colSums(egc_type2)
        tot_type1[tot_type1 < 100000] <- NA
        tot_type2[tot_type2 < 100000] <- NA
        egc_type1 <- t(t(egc_type1) / tot_type1)
        egc_type2 <- t(t(egc_type2) / tot_type2)
        imp_type[[mut]] <- list(egc_type1[imprint$gname, ], egc_type2[imprint$gname, ])
        f1 <- colSums(mut_l$egc1) > 80000
        f2 <- colSums(mut_l$egc2) > 80000
        for (g1 in names(best_mate)) {
            g2 <- best_mate[g1]
            max_x <- imprint[g1, "ymax"]
            max_y <- imprint[g2, "ymax"]
            png(sprintf(
                "paper_figs/imprints/%s.%s.%s.%s.png",
                g1, g2, mut, ls[1]
            ))
            plot(log2(mut_l$egc1_n[g1, f1] + eps), log2(mut_l$egc1_n[g2, f1] + eps), xlab = g1, ylab = g2, pch = 19, col = mc@colors[f1], cex = 1.6, xlim = c(-16.7, max_x), ylim = c(-16.7, max_y))
            dev.off()
            png(sprintf(
                "paper_figs/imprints/%s.%s.%s.%s.png",
                g1, g2, mut, ls[2]
            ))
            plot(log2(mut_l$egc2_n[g1, f2] + eps), log2(mut_l$egc2_n[g2, f2] + eps), xlab = g1, ylab = g2, pch = 19, col = mc@colors[f2], cex = 1.6, xlim = c(-16.7, max_x), ylim = c(-16.7, max_y))
            dev.off()
        }
    }
    return(imp_type)
}

barplot_imprints <- function(imp_type) {
    imprint <- fread("config/imprinted.txt")
    rownames(imprint) <- imprint$gname

    for (gname in imprint$gname) {
        png(sprintf("paper_figs/imprints/all_%s.png", gname), w = 1200, h = 300)
        stat <- data.frame(
            type = colnames(imp_type[["wt"]][[1]]),
            wt1 = imp_type[["wt"]][[1]][gname, ]
        )
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["wt"]][[2]]),
            wt2 = imp_type[["wt"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["3a"]][[2]]),
            d3a1 = imp_type[["3a"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["3a"]][[2]]),
            d3a2 = imp_type[["3a"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["3b"]][[2]]),
            d3b1 = imp_type[["3b"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["3b"]][[2]]),
            d3b2 = imp_type[["3b"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["dko"]][[2]]),
            dko1 = imp_type[["dko"]][[2]][gname, ]
        ))
        stat <- stat %>% left_join(data.frame(
            type = colnames(imp_type[["dko"]][[2]]),
            dko2 = imp_type[["dko"]][[2]][gname, ]
        ))
        barplot(as.matrix(log2(stat[, -1] + 1e-5)) - log2(1e-5), beside = T, col = stat$type, main = gname)
        grid()
        dev.off()
    }
}
