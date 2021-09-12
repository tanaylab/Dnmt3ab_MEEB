
#we use ebdnmt_module
source("paper_scripts/ebdnmt_util_plots.r")
source("paper_scripts/ebdnmt_module_analysis.r")

tgconfig::override_params("config/eseb.yaml", package="metacell")

do_fig2 = F
do_fig3 = F

d = load_data()
gmods = gen_gmods(d)

plot_gmods_cors(d, gmods)

plot_cmp_top_gmods(d, gmods, n_genes=40, exclude_genes=cmp_exclude_genes, text_size=14, force=30, w=1000, h=1000)

plot_gmods_cors(d, gmods, lead_g=c("Foxa1", "Mesp1", "Eomes", "Utf1", "Tal1", "Dppa3", "Hand1", "Foxc1"))

if(do_fig2) {
ebdnmt_plot_kit(mc_id = "eseb_07_3ab_bs500f_got_s2m",
            mat_id = "eseb_07_3ab",
            mc2d_id="eseb_07_3ab_bs500f_got_s2m",
            base_dir = "paper_figs/fig2/dko",
				do_gene_maps=T)

ebdnmt_plot_kit(mc_id = "eseb_07_TKO_bs500f_got_s2m",
            mat_id = "eseb_07_TKO",
            mc2d_id="eseb_07_TKO_bs500f_got_s2m",
            base_dir = "paper_figs/fig2/tko",
				do_gene_maps=T)

#plotting epiblast/gene comparisons
if(!dir.exists("paper_figs/fig2/cmp_mod_g")) {
	dir.create("paper_figs/fig2/cmp_mod_g")
}

wt_l = comp_egc_on_lines(d, "wt", "J1", "N15")
f_cov = colSums(wt_l$egc1)>80000
wt_j1_legc = log2(1e-5+wt_l$egc1_n)
for(g in c("Dnmt3l", "Tet1", "Lefty1", "Lefty2", "Nodal", "Ifitm1", "Ifitm3", "Dppa3")) {
#restrict to j1 wt
	plot_g_gmod_compare(d, gmods=gmods, gmod="epi", gene=g, 
				legc1=wt_j1_legc[,f_cov], legc2=d$legc_dko,
				col1=d$mc_wt@colors[f_cov], col2=d$mc_dko@colors, 
				base_dir="paper_figs/fig2/cmp_mod_g/", suf="_dko")
}

#testing genes on specific lines

dko_l = strat_gmod_on_lines(d, "dko", "DKO", "DKO16",
									gmods[["epi"]], bin_vals=seq(0,60,10))
wt_l = strat_gmod_on_lines(d, "wt", "J1", "N15", 
									gmods[["epi"]],bin_vals=seq(0,60,10))
for(g in c("Lefty1", "Lefty2", "Nodal", "Cer1", "Cfc1", "Tdgf1","Eomes", "Tal1", "Foxa2", "Mesp1")) {
	png(sprintf("paper_figs/fig2/epi_stratif_lines/%s.png",g), w=300,h=300)
	plot(log2(dko_l$strati1_n[g,]+1e-5), type="l", col="darkblue", 
								ylim=c(-16.7,-10), lwd=3, xlab="epi module", ylab=g)
	lines(log2(dko_l$strati2_n[g,]+1e-5), type="l",col="darkred", lwd=3)
	lines(log2(wt_l$strati1_n[g,]+1e-5), type="l", lty=1, col="gray",lwd=3)
	dev.off()
}
#distribution of program per time and line : move from cmp dir to fig dir

#generating comparison of epi and dppa3

#generating comparison of gene program for wt vs dko


} #if(do_fig2)


if(do_fig3) {
	ebdnmt_plot_kit(mc_id = "eseb_07_3a_bs500f_got_s2m",
            mat_id = "eseb_07_3a",
            mc2d_id="eseb_07_3a_bs500f_got_s2m",
            base_dir = "paper_figs/fig3/3a",
				do_gene_maps=T)

	ebdnmt_plot_kit(mc_id = "eseb_07_3b_bs500f_got_s2m",
            mat_id = "eseb_07_3b",
            mc2d_id="eseb_07_3b_bs500f_got_s2m",
            base_dir = "paper_figs/fig3/3b",
				do_gene_maps=T)

#for figure 3ab
	for(nm in names(gmods)) { 
		plt_cmp_3a3b_lines(d, gmods, prog=nm, k_cells=160, text_size=10, force=20)
	}

}

plot_ab_rna_compare = function()
{
#plot 3a 3b RNA comparisons
	iv_mc = scdb_mc("sing_emb_wt10_recolored")
	eb_mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")

	layer_colors = c("#65A83E","#65A83E", "#635547", "#DABE99", "#EF5A9D","#C594BF","#cc7818","#FF0000", "pink")

	names(layer_colors) = c("Ecto","EctoMeso", "Epi","MesoEndo", "Endo", "Meso", "ExeMeso", "Heme","EXE")
	type_key = read.table("config/atlas_type_order_temporal.txt", 
							sep="\t", h=T, stringsAsFactors=F)
	color_layer = type_key$layer
	names(color_layer) = type_key$color

	mc_layer_col = layer_colors[color_layer[iv_mc@colors]]
	eb_mc_layer_col = layer_colors[color_layer[eb_mc@colors]]

	iv_legc = log2(iv_mc@e_gc)

	eb_legc = log2(eb_mc@e_gc)

	min_stage = 4
	md = scdb_mat("sing_emb_wt10")@cell_metadata
	max_t = max(md$age_group,na.rm=T)
	mc_t = table(iv_mc@mc, md[names(iv_mc@mc),"age_group"])
	mc_t_e = apply(mc_t,1,function(x) return(sum(x*(1:max_t))/sum(x)))
	f = mc_t_e > 3

	eb_md = scdb_mat("eseb_07_wt")@cell_metadata
	max_t = 7
	eb_md$EB_day = as.character(eb_md$EB_day)
	mc_t = table(eb_mc@mc, eb_md[names(eb_mc@mc),"EB_day"])
	mc_t_e = apply(mc_t,1,function(x) return(sum(x*(0:max_t))/sum(x)))
	f_eb = mc_t_e > 3

	png("paper_figs/fig3/rna_3a3b.png", w=1000,h=500)
	layout(matrix(1:2,nrow=1))
	plot(iv_legc["Dnmt3a",f], iv_legc["Dnmt3b",f], pch=21, bg=mc_layer_col[f], 
					xlim=c(-17,-10), ylim=c(-15,-7.5), cex=1.5)
	plot(eb_legc["Dnmt3a",f_eb], eb_legc["Dnmt3b",f_eb], pch=21, 
				bg=eb_mc_layer_col[f_eb], xlim=c(-17,-10), ylim=c(-15,-7.5),cex=1.5)
	dev.off()
}
