gen_sfig2_atlas = function()
{

#gotg atlas
	atlas_key = read.table("config/atlas_type_order.txt", h=T, stringsAsFactors=F)
	color_ord = 1:nrow(atlas_key)
	names(color_ord) = atlas_key$colour
	mat = scdb_mat("emb_gotg_75")
	mc = scdb_mc("emb_gotg_75_bs500f")
	mc_day = table(mc@mc, as.numeric(as.factor(mat@cell_metadata[names(mc@mc), "stage"])))
	mc_day_max = apply(mc_day,1, function(x) sum(x*0:5)/sum(x))
	mc_ord = order(color_ord[mc@colors] + 0.01*mc_day_max)

	mcell_mc_plot_marks("emb_gotg_75_bs500f", "eseb_key_tfs", mat="emb_gotg_75", 
					plot_cells=F, reorder_marks=T, fold_burn=1.2, mc_ord=mc_ord)

	mc2d = scdb_mc2d("emb_gotg_75_bs500f")
	mc = scdb_mc("emb_gotg_75_bs500f")
	png("paper_figs/sfig2/ed1_gotg_2d.png", w=800, h=800)
	plot(mc2d@mc_x, mc2d@mc_y, pch=21, bg=mc@colors, xlab=NA, ylab=NA)
	segments(mc2d@mc_x[mc2d@graph$mc1], mc2d@mc_y[mc2d@graph$mc1], 
					mc2d@mc_x[mc2d@graph$mc2], mc2d@mc_y[mc2d@graph$mc2])
	points(mc2d@mc_x, mc2d@mc_y, pch=21, bg=mc@colors, cex=1.5)
	dev.off()

}

gen_sfig2_cmp_proj_and_epi = function()
{
	browser()
	mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")
	legc = log2(mc@e_gc)
  	mat = scdb_mat("eseb_07_wt")
	md= mat@cell_metadata
	mc_d = table(mc@mc, as.character(md[names(mc@mc),"EB_day"]))
head(mc_d)
	mc_d_e = apply(mc_d, 1, function(x) return(sum(x*(0:7))/sum(x)))
	f = mc_d_e > 3.5 & mc@colors != "#635547" & mc@colors != "darkkhaki"
	
	s2m_proj_egc = mcatlas_annotate_mc_by_sc2mc_projection("emb_gotg_75", 
									"eseb_07_wt", "eseb_07_wt_bs500f_got_s2m", 
									qmat_naming_type="mars")
	got_legc = log2(s2m_proj_egc+1e-5)
	s2m_proj_egc = mcatlas_annotate_mc_by_sc2mc_projection("sing_emb_wt10", 
									"eseb_07_wt", "eseb_07_wt_bs500f_tem_s2m", 
									qmat_naming_type="mars")
	tem_legc = log2(s2m_proj_egc+1e-5)

	gs = c("Utf1", "Tdgf1", "T", "Eomes", "Snai1", "Foxa2", "Epcam", "Dnmt3a", "Sox17", "Msx1", "Hand2", "Pitx1", "Prrx2", "Twist1", "Fn1", "Dnmt3b", "Emb", "Bmp4", "Igf2", "Lefty2", "Cer1", "Fgf5")
	
	png("paper_figs/sfig2/got_cmp_gs.png",w=1600,h=600)
	layout(matrix(1:24, byrow=T,nrow=3))
	for(g in gs) {
		plot(got_legc[g,f], legc[g,f], pch=19, col=mc@colors[f], cex.main=2, main=g, xlab=NA, ylab=NA, cex.axis=1.2)
	}
	dev.off()
	png("paper_figs/sfig2/tem_cmp_gs.png", w=1600, h=600)
	layout(matrix(1:24, byrow=T,nrow=3))
	for(g in gs) {
		plot(tem_legc[g,f], legc[g,f], pch=19, col=mc@colors[f], cex.main=2, main=g, xlab=NA, ylab=NA, cex.axis=1.2)
	}
	dev.off()

	gs = c("Klf4", "Esrrb", "Ifitm1", "Ifitm3", "Dppa3", "Sox2", "Dnmt3a", "Dnmt3b","Eomes", "Fst")
	png("paper_figs/sfig2/epi_dynamics.png",w=1200,h=400)
	layout(matrix(1:12, byrow=T,nrow=2))
	par(mar=c(3,3,2,1))
	for(g in gs) {
		plot(mc_d_e[!f], legc[g,!f], pch=19, col=mc@colors[!f], cex=1.5, xlab=NA, ylab=NA, main=g, cex.main=2,cex.axis=1.5)
	}
	dev.off()
}

plot_compare_eb_to_iv_states = function()
{
	
	a = mcatlas_annotate_mc_by_sc2mc_projection("sing_emb_wt10", 
									"eseb_07_wt", "eseb_07_wt_bs500f_tem_s2m", 
									qmat_naming_type="mars")

	mc = scdb_mc("eseb_07_wt_bs500f_tem_s2m")
	egc = mc@e_gc[rownames(a),]
	
	atl_cor = c()
	for(i in 1:ncol(a)) { 
		atl_cor = c(atl_cor, cor(wt_egc[,i],a[,i])) 
	}
	names(atl_cor) = 1:length(atl_cor)

	best_mc_hit = tapply(atl_cor, d$mc_wt@colors, 
										function(x) names(x)[which.max(x)])

#epi, ps, nm, aps
	foc_colors = c("#635547", "#DABE99", "#C594BF", "#c19f70")

	best_mc_hit = best_mc_hit[foc_colors]
	png("paper_figs/sfig1/cmp_eb_iv_mcs.png",w=1000,h=1000)
	layout(matrix(1:16,nrow=4))
	par(mar=c(2,2,2,2))
	for(i in 1:4) {
		for(j in 1:4) {
			pair_c = cor(log2(1e-5+egc[,best_mc_hit[i]]), log2(1e-5+a[,best_mc_hit[j]]))
			plot(log2(1e-5+egc[,best_mc_hit[i]]), log2(1e-5+a[,best_mc_hit[j]]), pch=19, col=ifelse(i==j, foc_colors[i], "gray"), cex=1, main=sprintf("r=%.3f",pair_c), xlab=NA,ylab=NA)
		}
	}
	dev.off()
}
	
