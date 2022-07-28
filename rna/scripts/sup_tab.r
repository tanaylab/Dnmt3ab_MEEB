
; 
library(metacell)
scdb_init("scrna_db/", force_reinit=T); 
scfigs_init("figs/")

mats = c("eseb_07_wt", "eseb_07_3a", "eseb_07_3b", "eseb_07_3ab", "eseb_07_TKO")

all = NULL
for(nm in mats) {
	cur_type = sub("eseb_07_", "", nm)
	mat = scdb_mat(nm)
	csize = colSums(mat@mat)
	md = mat@cell_metadata[colnames(mat@mat),]
	batches = unique(md[,c("amp_batch_id", "line", "EB_day", "rep")])
	rownames(batches) = batches[,"amp_batch_id"]
	batch_csize_med = tapply(csize, md$amp_batch_id, median)
	batch_n = tapply(csize, md$amp_batch_id, length)
	batches$type = cur_type
	batches$umi_med = batch_csize_med[batches$amp_batch_id]
	batches$n_qc_cell = batch_n[batches$amp_batch_id]
	if(is.null(all)) {
		all = batches
	} else {
		all = rbind(all, batches)
	}
}

write.table(all,"paper_figs/table_s1.txt", quote=F, sep="\t")
