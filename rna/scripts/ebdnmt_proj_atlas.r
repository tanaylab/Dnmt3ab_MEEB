
library(metacell)
library(Matrix)
scdb_init("scrna_db", force_reinit=T)
scfigs_init("figs/")

source("scripts/ebdnmt_proj_atlas_util.r")

tgconfig::override_params("config/eseb.yaml",package="metacell")

recomp=F
if(1) {
	atlas_annot_w_all("eseb_07_wt_bs500f", qmat_id = "eseb_07_wt", recomp_knn=recomp, max_entropy=2.5)
}
if(1) {
	atlas_annot_w_all("eseb_07_3a_bs500f", qmat_id = "eseb_07_3a", recomp_knn=recomp, max_entropy=2.5)
}
if(1) {
	atlas_annot_w_all("eseb_07_3b_bs500f", qmat_id = "eseb_07_3b", recomp_knn=recomp, max_entropy=2.5)
}
if(0) {
	atlas_annot_w_all("eseb_07_3ab_bs500f", qmat_id = "eseb_07_3ab", recomp_knn=recomp, max_entropy=2.5)
}
if(0) {
	atlas_annot_w_all("eseb_07_TKO_bs500f", qmat_id = "eseb_07_TKO", recomp_knn=recomp, max_entropy=2.5)
}

