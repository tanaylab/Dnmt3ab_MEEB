
gen_cc_gmods = function()
{
	mat = scdb_mat("eseb_07_wt")
	mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")
	md = mat@cell_metadata
	mat_ds = scm_downsamp(mat@mat, 3000)

	cols = unique(mc@colors)
	all_pcna = data.frame(name = rownames(mat_ds))
	all_mki67= data.frame(name = rownames(mat_ds))
	allnames = c("gene")
	for(col in cols) {
	  col_cells = names(mc@mc)[mc@colors[mc@mc]==col]
		for(t in 1:7) {
			cells = intersect(col_cells, rownames(md)[md$EB_day==t])
			cells = intersect(colnames(mat_ds), cells)
			if(length(cells) > 100) {
				#cor pcna, Ki6
				lumis = log2(1+7*mat_ds[,cells])
				c_mki67 = apply(lumis, 1, cor, lumis["Mki67",])
				c_pcna = apply(lumis, 1, cor, lumis["Pcna",])
				c_mki67[c_mki67>0.9]=0
				c_pcna[c_pcna>0.9]=0
				coltag = sub("#","",col)
				message("doing ", coltag, " t ", t, " n = ", length(cells))
		
				rnk1 = rank(-abs(c_mki67))	
				rnk2 = rank(-abs(c_pcna))
				df = data.frame(mki67=c_mki67, pcna=c_pcna, lab=names(c_mki67),
							col = ifelse(rnk1<40 | rnk2<40, col, "gray"))
				p = ggplot(df, aes(mki67, pcna)) + 
					geom_point(col=df$col, size=5) +
					geom_text_repel(data = filter(df, df$col!="gray"), 
							size=8,force=10, max.iter=1e+4,aes(label=lab))
				png(sprintf("paper_figs/cc/%s.%s.png", coltag, t), w=800, h=800)
				print(p)
				dev.off()
				all_pcna = cbind(all_pcna, c_pcna)
				all_mki67= cbind(all_mki67, c_mki67)
				allnames = c(allnames, sprintf("c%s.%s", coltag, t))
			}
		}
	}
	colnames(all_pcna) = allnames
	colnames(all_mki67) = allnames
	gmod_m = c("Mki67", names(tail(sort(rowMeans(all_mki67[,-1])),30)))
	gmod_s = c("Pcna", names(tail(sort(rowMeans(all_pcna[,-1])),30)))
	tot_m = colSums(mat_ds[gmod_m,])
	tot_s = colSums(mat_ds[gmod_s,])

	for(col in cols) {
	  col_cells = names(mc@mc)[mc@colors[mc@mc]==col]
		for(t in 1:7) {
			cells = intersect(col_cells, rownames(md)[md$EB_day==t])
			cells = intersect(colnames(mat_ds), cells)
			if(length(cells) > 100) {
				coltag = sub("#","",col)
				foc_mat = mat_ds[,cells]
				foc_s = tot_s[cells]
				foc_m = tot_m[cells]
				png(sprintf("paper_figs/cc/phase.%s.%s.png", coltag, t), h=400, w=1200)
				layout(matrix(1:3,nrow=1))
				plot(foc_m, foc_s, pch=19, col=col, cex=1)
				fmean = function(x) { return(ifelse(length(x)>10, mean(x), NA)) }
				b_prf1 = tapply(foc_mat["Dnmt3b",],cut(foc_s,seq(10,120,10)), fmean)
				b_prf2 = tapply(foc_mat["Dnmt3b",],cut(foc_m,seq(0,60,6)),fmean)
				max_y = max(c(5,b_prf1, b_prf2),na.rm=T)
				plot(b_prf1, type="b", ylim=c(0,max_y), col="gold", lwd=2, xlab="s_score")
				lines(tapply(foc_mat["Dnmt3a",],cut(foc_s,seq(10,120,10)),fmean), col="magenta", lwd=2)
				lines(tapply(foc_mat["Dnmt3l",],cut(foc_s,seq(10,120,10)),fmean), col="black", lwd=2)
				lines(tapply(foc_mat["Uhrf1",],cut(foc_s,seq(10,120,10)),fmean), col="darkgreen", lwd=2)
				plot(b_prf2, type="b", ylim=c(0,max_y), col="gold", lwd=2, xlab="m_score")
				lines(tapply(foc_mat["Dnmt3a",],cut(foc_m,seq(0,60,6)), fmean), col="magenta",lwd=2)
				lines(tapply(foc_mat["Dnmt3l",],cut(foc_m,seq(0,60,6)), fmean), col="black",lwd=2)
				lines(tapply(foc_mat["Uhrf1",],cut(foc_m,seq(0,60,6)), fmean), col="darkgreen",lwd=2)
				dev.off()
			}
		}
	}
}

plot_all_cc = function()
{
	ivcc = fread("/net/mraid14/export/tgdata/users/aviezerl/proj/ebdnmt/output/cell_cycle/invivo_cgc.tsv")
	ivcc_ab = fread("/net/mraid14/export/tgdata/users/aviezerl/proj/ebdnmt/output/cell_cycle/invivo_ab_score.tsv")
	ebcc = fread("/net/mraid14/export/tgdata/users/aviezerl/proj/ebdnmt/output/cell_cycle/eb_cgc.tsv")
	ebcc_ab = fread("/net/mraid14/export/tgdata/users/aviezerl/proj/ebdnmt/output/cell_cycle/eb_ab_score.tsv")

	for(layer in c("meso", "ecto", "endo")) {
		png(sprintf("paper_figs/fig5/e75_cc_el_%s.png", layer), w=1000,h=300)
		plot_mm_cc(ivcc, layer=layer, d="e7.5", cg="(0,0.02]", ylim=c(0.8,0.95))
		dev.off()

		png(sprintf("paper_figs/fig5/e75_cc_ab_%s.png", layer), w=1000,h=300)
		plot_mm_ab_dlt(ivcc_ab, layer=layer, d="e7.5", cg="(0,0.02]", ylim=c(-0.03,0.06))
		dev.off()
	}
	f = (0.99-ebcc$early_late_cov*0.11-ebcc$avg_early)<0
	cells = unique(ebcc$cell_id[!f])
	ebcc = ebcc[ebcc$cell_id %in% cells,]
	ebcc_ab = ebcc_ab[ebcc_ab$cell_id %in% cells,]
	for(line in c("wt", "ko3a", "ko3b")) {
		png(sprintf("paper_figs/fig5/eb_cc_el_%s.png", line), w=1000,h=300)
		plot_mm_cc(ebcc, foc_line=line, d="d6", cg="(0,0.02]", ylim=c(0.8,0.95))
		dev.off()
		png(sprintf("paper_figs/fig5/eb5_cc_el_%s.png", line), w=1000,h=300)
		plot_mm_cc(ebcc, foc_line=line, d="d5", cg="(0,0.02]", ylim=c(0.8,0.95))
		dev.off()

		png(sprintf("paper_figs/fig5/eb_cc_ab_%s.png", line), w=1000,h=300)
		plot_mm_ab_dlt(ebcc_ab, foc_line=line, d="d6", cg="(0,0.02]", ylim=c(-0.03,0.06))
		dev.off()
		png(sprintf("paper_figs/fig5/eb5_cc_ab_%s.png", layer), w=1000,h=300)
		plot_mm_ab_dlt(ebcc_ab[cells,], foc_line=line, d="d5", cg="(0,0.02]", ylim=c(-0.03,0.06))
		dev.off()
	}

}

plot_mm_cc = function(tab, layer = NA, d, cg=NA, foc_line = NA, ylim=NULL, foc_cells = NULL)
{
	tab = tab %>% filter(day==d)
	if(!is.na(layer)) {
		tab = tab %>% filter(germ_layer==layer)
	}
	if(!is.na(cg)) {
		tab = tab %>% filter(cg_cont==cg)
	}
	if(!is.na(foc_line)) {
		tab = tab %>% filter(line==foc_line)
	}

	shades = colorRampPalette(c("white", "blue", "red", "yellow"))(1000)

	layout(matrix(c(1:3),nrow=1,byrow=T), w=c(1,1,3))

	trend = rollmean(tab$early_late_cov[order(tab$ord1)],20)
	i_mid = rollmean(tab$ord1[order(tab$ord1)],20)[which.min(trend)]

	par(mar=c(2,2,2,1))	
	message("i_mid ", i_mid)	
	tab$ord2 = (i_mid - tab$ord1)
	tab$ord2 = tab$ord2-floor(tab$ord2)
	plot(tab$avg_early, tab$avg_late, pch=21,bg=shades[1+999*tab$ord2],xlab=NA)
	par(mar=c(2,1,2,1))	
	plot(tab$avg_late - tab$avg_early, tab$early_late_cov, pch=21,bg=shades[1+999*tab$ord2],xlab=NA)

	par(mar=c(2,1,0,2))	
	plot(tab$ord2, tab$avg_early, pch=19, col="darkblue", 
				ylim=c(min(c(tab$avg_early,tab$avg_late)), 
							max(c(tab$avg_early, tab$avg_late))))
	points(tab$ord2, tab$avg_late, pch=19, col="cyan")
	
}

plot_mm_ab_dlt = function(tab, layer = NA, d, cg=NA, foc_line = NA, ylim=NULL, foc_cells= NULL)
{
	tab = tab %>% filter(day==d)
	if(!is.na(layer)) {
		tab = tab %>% filter(germ_layer==layer)
	}
	if(!is.na(cg)) {
		tab = tab %>% filter(cg_cont==cg)
	}
	if(!is.na(foc_line)) {
		tab = tab %>% filter(line==foc_line)
	}
	if(!is.null(foc_cells)) {
		tab = tab[tab$cell_id %in% foc_cells,]
	}
	tab_high_ab = tab %>% filter(ab_score=="(0.267,1]")
	tab_low_ab = tab %>% filter(ab_score=="(-1.24,-0.709]")
	tab_high_ab = tab_high_ab %>% mutate(hab_late = avg_late, hab_early=avg_early)
	tab_low_ab = tab_low_ab %>% mutate(lab_late = avg_late, lab_early=avg_early)
	tab_merg = tab_low_ab %>% select(cell_id, cg_cont, day, ord1, early_late_cov, lab_late, lab_early) %>% left_join(tab_high_ab %>% select(cell_id, cg_cont, day, hab_late, hab_early))

	tab = tab_merg
	if(nrow(tab) > 50) {

	shades = colorRampPalette(c("white", "blue", "red", "yellow"))(1000)

	trend = rollmean(tab$early_late_cov[order(tab$ord1)],20)
	i_mid = rollmean(tab$ord1[order(tab$ord1)],20)[which.min(trend)]

	message("i_mid ", i_mid)	
	tab$ord2 = (i_mid - tab$ord1)
	tab$ord2 = tab$ord2-floor(tab$ord2)

	tab$d_early = (tab$hab_early-tab$lab_early)
	tab$d_late= (tab$hab_late-tab$lab_late)
	grp_early = split(tab$d_early, cut(tab$ord2, c(0,0.33,0.66,1)))	
	grp_late= split(tab$d_late, cut(tab$ord2, c(0,0.33,0.66,1)))	
	ks_e = ks.test(grp_early[[1]], grp_early[[3]])
	ks_l = ks.test(grp_late[[1]], grp_late[[3]])

	if(is.null(ylim)) {
		ylim = c(min(tab$d_early, tab$d_late), max(tab$d_early, tab$d_late))
	}	
	plot(tab$ord2, tab$hab_early-tab$lab_early, pch=19, col="darkblue", cex=0.6,
						main = sprintf("e pv %e l pv %e", ks_e$p.value, ks_l$p.value),
						ylim=ylim)
	points(tab$ord2, tab$hab_late-tab$lab_late, pch=19, col="cyan", cex=0.6)
	boxplot(split(tab$d_early, cut(tab$ord2, seq(0,1,l=10))), add=T, 
								at=seq(0.02,0.92,l=9), boxwex=0.02, col="darkblue")
	boxplot(split(tab$d_late, cut(tab$ord2, seq(0,1,l=10))), add=T, 
								at=seq(0.06,0.96,l=9), boxwex=0.02, col="cyan")
	grid()

	}
	return(tab)

}
