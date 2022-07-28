
library("devtools"); 
library(metacell)
scdb_init("scrna_db/", force_reinit=T); 
scfigs_init("figs/")

source("paper_scripts/ebdnmt_util_plots.r")

tgconfig::override_params("config/eseb.yaml", package="metacell")

ebdnmt_plot_kit(mc_id = "eseb_07_wt_bs500f_got_s2m",
            mat_id = "eseb_07_wt",
            mc2d_id="eseb_07_wt_bs500f_got_s2m",
            base_dir = "paper_figs/fig1/")

if(0) {
mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")

#plot TF heat map

atlas_key = read.table("config/atlas_type_order.txt", h=T, stringsAsFactors=F)
color_ord = 1:nrow(atlas_key)
names(color_ord) = atlas_key$colour
mat = scdb_mat("eseb_07_wt")
mc_day = table(mc@mc, mat@cell_metadata[names(mc@mc), "EB_day"])
mc_day_max = apply(mc_day,1, function(x) sum(x*0:7)/sum(x))
mc_ord = order(color_ord[mc@colors] + 0.01*mc_day_max)

mcell_mc_plot_marks("eseb_07_wt_bs500f_got_s2m", "eseb_key_tfs", 
			mat="eseb_07_wt", plot_cells=F, 
			mc_ord = mc_ord,
			fig_fn = "paper_figs/fig1/wt_tfs_heat.png",
			reorder_marks=F, fold_burn=1.5)

#plot 2D without ESCs and singletons
f_mc = mc@colors!="antiquewhite3" & mc@colors!="antiquewhite1"

mgr = scdb_mc2d("eseb_07_wt_bs500f_got_s2m")
f_gr = f_mc[mgr@graph$mc1] & f_mc[mgr@graph$mc2]
deg = tabulate(c(mgr@graph$mc1[f_gr], mgr@graph$mc2[f_gr]))

f_mc = f_mc & deg!=0

tgconfig::set_param("mcell_mc2d_cex",0.7,"metacell")

mcell_mc2d_plot("eseb_07_wt_bs500f_got_s2m",show_mcid=F,cell_outline=T, sc_cex=0.8,edge_w=0.4, filt_mc=f_mc)

# plot type time trace

plot_type_trace("eseb_07_wt", "eseb_07_wt_bs500f_got_s2m")

# plot gene 2D projections

plot_gene_maps(mat_id = "eseb_07_wt", 
					mc_id = "eseb_07_wt_bs500f_got_s2m", 
					mc2d_id = "eseb_07_wt_bs500f_got_s2m",
					base_dir
					filt_mc = f_mc)
}
