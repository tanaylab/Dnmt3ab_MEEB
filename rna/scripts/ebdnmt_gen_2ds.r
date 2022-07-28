library("devtools")
library(metacell)
scdb_init("scrna_db", force_reinit=T)
scfigs_init("figs/")

source("paper_scripts/ebdnmt_proj_atlas_util.r")

tgconfig::override_params("config/eseb.yaml",package="metacell")

#source("scripts/mars_10x_translate.r")
do_wt = 1
do_3a = 0
do_3b = 0
do_3ab = 0
do_D1 = 0
do_TKO = 0
do_wt_sing = 0

if(do_wt) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.015, "metacell")
tgconfig::set_param("mcell_mc2d_K", 30, "metacell")
tgconfig::set_param("mcell_mc2d_max_confu_deg", 4, "metacell")

mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")
lfp = log2(mc@mc_fp)
legc = log2(mc@e_gc +1e-5)

m = scdb_mat("eseb_07_wt")

md = m@cell_metadata

mc_wt_d = table(mc@mc, md[names(mc@mc), "EB_day"])

mc_wt_dscore = apply(mc_wt_d, 1, function(x) mean(rep(0:7, times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & legc["Dppa3",]> -14 & legc["Klf4",]< -13, "darkkhaki", mc@colors)
mc@colors[legc["Psme1",] > -11.5 & legc["Sox2",]< -14] = "#FACB12" 
mc@colors = ifelse(mc_wt_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[legc["Zscan4c",]> -15] = "antiquewhite1"

scdb_add_mc("eseb_07_wt_bs500f_got_s2m", mc)
cut_edges = data.frame(mc1=c(496,491,491),
                        mc2=c(29,120,79))
mcell_mc2d_force_knn("eseb_07_wt_bs500f_got_s2m", "eseb_07_wt_bs500f_got_s2m", 
				graph_id="eseb_07_wt",ignore_edges=cut_edges)

#mcell_mc2d_rotate("eseb_07_wt_bs500f_got_s2m", alpha=-180)

mcell_mc2d_plot("eseb_07_wt_bs500f_got_s2m")

cut_mcs = unique(c(cut_edges$mc1,cut_edges$mc2))
mcell_mc_plot_submc_marks("eseb_07_wt_bs500f_got_s2m", "eseb_07_wt", 
					foc_mcs=cut_mcs,
					fig_fn="paper_figs/wt_merge_edge.png",n_max_marks=15, w=900)

mcell_mc_plot_marks("eseb_07_wt_bs500f_got_s2m", "eseb_07_wt_bs500f", "eseb_07_wt", add_metadata="EB_day")

mcell_mc2d_plot_by_factor("eseb_07_wt_bs500f_got_s2m","eseb_07_wt", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_wt_bs500f_got_s2m","eseb_07_wt", "EB_day",
						ncols=4, single_plot=F, neto_points=T)

day_gel = sub("WT_","",md[names(mc@mc),"exp_state"])
md_vals = paste(md[names(mc@mc),"rep"], day_gel, sep="_")
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
names(md_vals) = names(mc@mc)
mcell_mc2d_plot_by_factor("eseb_07_wt_bs500f_got_s2m","eseb_07_wt", 
					meta_field="rep_and_day",
					meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)

atlas_key = read.table("config/atlas_type_order.txt", h=T, stringsAsFactors=F)
color_ord = 1:nrow(atlas_key)
names(color_ord) = atlas_key$colour
mat = scdb_mat("eseb_07_wt")
mc_day = table(mc@mc, mat@cell_metadata[names(mc@mc), "EB_day"])
mc_day_max = apply(mc_day,1, function(x) sum(x*0:7)/sum(x))
mc_ord = order(color_ord[mc@colors] + 0.01*mc_day_max)
mcell_mc_plot_marks("eseb_07_wt_bs500f_got_s2m", "eseb_key_tfs", "eseb_07_wt",mc_ord=mc_ord, plot_cells=F)

}
if(do_3b) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.016, "metacell")
tgconfig::set_param("mcell_mc2d_K", 50, "metacell")

mc = scdb_mc("eseb_07_3b_bs500f_got_s2m")
lfp = log2(mc@mc_fp)

m = scdb_mat("eseb_07_3b")

md = m@cell_metadata

mc_3b_d = table(mc@mc, md[names(mc@mc), "EB_day"])

mc_3b_dscore = apply(mc_3b_d, 1, function(x) mean(rep(0:7, times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & lfp["Dppa3",]>0.2 & lfp["Klf4",]<1, "darkkhaki", mc@colors)
mc@colors = ifelse(mc_3b_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[lfp["Zscan4c",]>0.5] = "antiquewhite1"

scdb_add_mc("eseb_07_3b_bs500f_got_s2m", mc)

cut_edges = data.frame(mc1=c(62,47,47),
									mc2=c(331,18,300))

mcell_mc2d_force_knn("eseb_07_3b_bs500f_got_s2m", "eseb_07_3b_bs500f_got_s2m", 
				graph_id="eseb_07_3b",
				ignore_edges=cut_edges)

mcell_mc2d_rotate("eseb_07_3b_bs500f_got_s2m", flipx=T)
mcell_mc2d_plot("eseb_07_3b_bs500f_got_s2m")

cut_mcs = unique(c(cut_edges$mc1,cut_edges$mc2))
mcell_mc_plot_submc_marks("eseb_07_3b_bs500f_got_s2m", "eseb_07_3b", 
					foc_mcs=cut_mcs,
					fig_fn="paper_figs/3b_merge_edge.png",n_max_marks=15, w=900)

mcell_mc_plot_marks("eseb_07_3b_bs500f_got_s2m", "eseb_07_3b_bs500f", "eseb_07_3b", add_metadata="EB_day")

mcell_mc2d_plot_by_factor("eseb_07_3b_bs500f_got_s2m","eseb_07_3b", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_3b_bs500f_got_s2m","eseb_07_3b", "EB_day",
						ncols=4, single_plot=F,neto_points=T)

md_vals = paste(md[names(mc@mc),"rep"], md[names(mc@mc),"EB_day"], sep="_")
names(md_vals) = names(mc@mc)
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
mcell_mc2d_plot_by_factor("eseb_07_3b_bs500f_got_s2m","eseb_07_3b", meta_field="rep_and_day",
								meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)

}

if(do_3a) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.016, "metacell")

mc = scdb_mc("eseb_07_3a_bs500f_got_s2m")
lfp = log2(mc@mc_fp)

m = scdb_mat("eseb_07_3a")

md = m@cell_metadata

day = as.integer(as.character(md[names(mc@mc), "EB_day"]))

mc_3a_d = table(mc@mc, day)

mc_3a_dscore = apply(mc_3a_d, 1, function(x) mean(rep(c(0:7), times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & lfp["Dppa3",]>0.2 & lfp["Klf4",]<1, "darkkhaki", mc@colors)
mc@colors = ifelse(mc_3a_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[lfp["Zscan4c",]>0.5] = "antiquewhite1"

scdb_add_mc("eseb_07_3a_bs500f_got_s2m", mc)

cut_edges=data.frame(mc1=c(240,81,186,186,179,315),
							mc2=c(185,180,180,315,183,373))
mcell_mc2d_force_knn("eseb_07_3a_bs500f_got_s2m", "eseb_07_3a_bs500f_got_s2m", 
			graph_id="eseb_07_3a",
			ignore_edges=cut_edges)

mcell_mc2d_rotate("eseb_07_3a_bs500f_got_s2m", alpha=-90)
mcell_mc2d_plot("eseb_07_3a_bs500f_got_s2m")

cut_mcs = unique(c(cut_edges$mc1,cut_edges$mc2))
mcell_mc_plot_submc_marks("eseb_07_3a_bs500f_got_s2m", "eseb_07_3a", 
					foc_mcs=cut_mcs,
					fig_fn="paper_figs/3a_merge_edge.png",n_max_marks=15, w=900)

mcell_mc_plot_marks("eseb_07_3a_bs500f_got_s2m", "eseb_07_3a_bs500f", "eseb_07_3a", add_metadata="EB_day")

mcell_mc2d_plot_by_factor("eseb_07_3a_bs500f_got_s2m","eseb_07_3a", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_3a_bs500f_got_s2m","eseb_07_3a", "EB_day",
						ncols=4, single_plot=F, neto_points=T)
md_vals = paste(md[names(mc@mc),"rep"], md[names(mc@mc),"EB_day"], sep="_")
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
names(md_vals) = names(mc@mc)
mcell_mc2d_plot_by_factor("eseb_07_3a_bs500f_got_s2m","eseb_07_3a", meta_field="rep_and_day",
								meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)

}

if(do_3ab) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.006, "metacell")
tgconfig::set_param("mcell_mc2d_K", 50, "metacell")
tgconfig::set_param("mcell_mc2d_max_confu_deg", 5, "metacell")

mc = scdb_mc("eseb_07_3ab_bs500f_got_s2m")
lfp = log2(mc@mc_fp)
legc = log2(mc@e_gc+1e-5)

m = scdb_mat("eseb_07_3ab")

md = m@cell_metadata

day = as.integer(as.character(md[names(mc@mc), "EB_day"]))

mc_3ab_d = table(mc@mc, day)

mc_3ab_dscore = apply(mc_3ab_d, 1, function(x) mean(rep(c(0:7), times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & lfp["Dppa3",]>0.2 & lfp["Klf4",]<1, "darkkhaki", mc@colors)
mc@colors = ifelse(mc_3ab_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[lfp["Zscan4c",]>0.5] = "antiquewhite1"
mc@colors[legc["Sox1",]> -15 & legc["Tdgf1",] < -12] = "chartreuse3"

scdb_add_mc("eseb_07_3ab_bs500f_got_s2m", mc)

cut_edges = data.frame(mc1=c(456), mc2=c(382))

mcell_mc2d_force_knn("eseb_07_3ab_bs500f_got_s2m", "eseb_07_3ab_bs500f_got_s2m", 
				graph_id="eseb_07_3ab", 
				ignore_edges=cut_edges)

cut_mcs = unique(c(cut_edges$mc1,cut_edges$mc2))
mcell_mc_plot_submc_marks("eseb_07_3ab_bs500f_got_s2m", "eseb_07_3ab", 
					foc_mcs=cut_mcs,
					fig_fn="paper_figs/3ab_merge_edge.png",n_max_marks=15, w=900)

mcell_mc_plot_marks("eseb_07_3ab_bs500f_got_s2m", "eseb_07_3ab_bs500f", "eseb_07_3ab", add_metadata="EB_day")

mcell_mc2d_rotate("eseb_07_3ab_bs500f_got_s2m", alpha=90)

mcell_mc2d_plot("eseb_07_3ab_bs500f_got_s2m")
mcell_mc_plot_marks("eseb_07_3ab_bs500f_got_s2m", "eseb_07_3ab_bs500f", "eseb_07_3ab", add_metadata="EB_day")

mcell_mc2d_plot_by_factor("eseb_07_3ab_bs500f_got_s2m","eseb_07_3ab", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_3ab_bs500f_got_s2m","eseb_07_3ab", "EB_day",
						ncols=4, single_plot=F,neto_points=T)

md_vals = paste(md[names(mc@mc),"rep"], md[names(mc@mc),"EB_day"], sep="_")
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
names(md_vals) = names(mc@mc)
mcell_mc2d_plot_by_factor("eseb_07_3ab_bs500f_got_s2m","eseb_07_3ab", meta_field="rep_and_day",
								meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)


}

if(do_D1) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.025, "metacell")
tgconfig::set_param("mcell_mc2d_K", 50, "metacell")
tgconfig::set_param("mcell_mc2d_max_confu_deg", 5, "metacell")

mc = scdb_mc("eseb_07_D1_bs500f_got_s2m")
lfp = log2(mc@mc_fp)

m = scdb_mat("eseb_07_D1")

md = m@cell_metadata

day = as.integer(as.character(md[names(mc@mc), "EB_day"]))

mc_D1_d = table(mc@mc, day)

mc_D1_dscore = apply(mc_D1_d, 1, function(x) mean(rep(c(0:7), times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & lfp["Dppa3",]>0.2 & lfp["Klf4",]<1, "darkkhaki", mc@colors)
mc@colors = ifelse(mc_D1_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[lfp["Zscan4c",]>0.5] = "antiquewhite1"
#mc@colors[lfp["Lama1",]>1] = "#1A1A1A"

scdb_add_mc("eseb_07_D1_bs500f_got_s2m", mc)

cut_edges = data.frame(mc1=c(227,226,103), mc2=c(225,259,256))

mcell_mc2d_force_knn("eseb_07_D1_bs500f_got_s2m", "eseb_07_D1_bs500f_got_s2m", 
				graph_id="eseb_07_D1",
				ignore_edges=cut_edges)

cut_mcs = unique(c(cut_edges$mc1,cut_edges$mc2))
mcell_mc_plot_submc_marks("eseb_07_D1_bs500f_got_s2m", "eseb_07_D1", 
					foc_mcs=cut_mcs,
					fig_fn="paper_figs/D1_merge_edge.png",n_max_marks=15, w=900)

mcell_mc2d_rotate("eseb_07_D1_bs500f_got_s2m", alpha=-90)

mcell_mc2d_plot("eseb_07_D1_bs500f_got_s2m")
mcell_mc_plot_marks("eseb_07_D1_bs500f_got_s2m", "eseb_07_D1_bs500f", "eseb_07_D1", add_metadata="EB_day")

mcell_mc2d_plot_by_factor("eseb_07_D1_bs500f_got_s2m","eseb_07_D1", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_D1_bs500f_got_s2m","eseb_07_D1", "EB_day",
						ncols=4, single_plot=F,neto_points=T)
md_vals = paste(md[names(mc@mc),"rep"], md[names(mc@mc),"EB_day"], sep="_")
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
names(md_vals) = names(mc@mc)
mcell_mc2d_plot_by_factor("eseb_07_D1_bs500f_got_s2m","eseb_07_D1", meta_field="rep_and_day",
								meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)

}

if(do_TKO) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.015, "metacell")
tgconfig::set_param("mcell_mc2d_K", 50, "metacell")
tgconfig::set_param("mcell_mc2d_max_confu_deg", 6, "metacell")

mc = scdb_mc("eseb_07_TKO_bs500f_got_s2m")
lfp = log2(mc@mc_fp)
legc = log2(mc@e_gc+1e-5)

m = scdb_mat("eseb_07_TKO")

md = m@cell_metadata

day = as.integer(as.character(md[names(mc@mc), "EB_day"]))

mc_TKO_d = table(mc@mc, day)

mc_TKO_dscore = apply(mc_TKO_d, 1, function(x) mean(rep(c(0:7), times=x)))

mc@colors = ifelse(mc@colors!= "#FACB12" & lfp["Dppa3",]>0.2 & lfp["Klf4",]<1, "darkkhaki", mc@colors)
mc@colors = ifelse(mc_TKO_dscore < 0.4, "antiquewhite3", mc@colors)
mc@colors[lfp["Zscan4c",]>0.5] = "antiquewhite1"
mc@colors[legc["Sox1",]> -15 & legc["Tdgf1",] < -12] = "chartreuse3"

scdb_add_mc("eseb_07_TKO_bs500f_got_s2m", mc)

cut_edges = data.frame(mc1=c(), mc2=c())

mcell_mc2d_force_knn("eseb_07_TKO_bs500f_got_s2m", "eseb_07_TKO_bs500f_got_s2m", 
				graph_id="eseb_07_TKO")

mcell_mc_plot_marks("eseb_07_TKO_bs500f_got_s2m", "eseb_07_TKO_bs500f", "eseb_07_TKO", add_metadata="EB_day")

mcell_mc2d_rotate("eseb_07_TKO_bs500f_got_s2m", alpha=-90)

mcell_mc2d_plot("eseb_07_TKO_bs500f_got_s2m")

mcell_mc2d_plot_by_factor("eseb_07_TKO_bs500f_got_s2m","eseb_07_TKO", "EB_day",
						ncols=4, single_plot=T, neto_points=T)
mcell_mc2d_plot_by_factor("eseb_07_TKO_bs500f_got_s2m","eseb_07_TKO", "EB_day",
						ncols=4, single_plot=F,neto_points=T)
md_vals = paste(md[names(mc@mc),"rep"], md[names(mc@mc),"EB_day"], sep="_")
md_vals = paste(md[names(mc@mc),"line"], md_vals, sep="_")
names(md_vals) = names(mc@mc)
mcell_mc2d_plot_by_factor("eseb_07_TKO_bs500f_got_s2m","eseb_07_TKO", meta_field="rep_and_day",
								meta_data_vals = md_vals, ncols=8, single_plot=F, neto_points=T)


}

if(do_wt_sing) {

tgconfig::set_param("mcell_mc2d_T_edge", 0.02, "metacell")
tgconfig::set_param("mcell_mc2d_K", 50, "metacell")
tgconfig::set_param("mcell_mc2d_max_confu_deg", 4, "metacell")

mc = scdb_mc("eseb_07_wt_bs500f_got_s2m")
lfp = log2(mc@mc_fp)
mcell_mc2d_force_knn("eseb_sEB8_wt_bs500f_got_s2m", "eseb_sEB8_wt_bs500f_got_s2m", 
				graph_id="eseb_sEB8_wt")

mcell_mc_plot_marks("eseb_sEB8_wt_bs500f_got_s2m", "eseb_sEB8_wt_bs500f", "eseb_sEB8_wt", add_metadata="rep")
mcell_mc2d_plot("eseb_sEB8_wt_bs500f_got_s2m")

mcell_mc2d_plot_by_factor("eseb_sEB8_wt_bs500f_got_s2m","eseb_sEB8_wt", "rep",
						ncols=4, single_plot=T, neto_points=T)
}
