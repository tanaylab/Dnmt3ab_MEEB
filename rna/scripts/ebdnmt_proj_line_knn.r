
if(0) {
	ref_mat_id = "eseb_07_wt"
	ref_mc_id = "eseb_07_wt_bs500f_got_s2m"
	q_mat_id = "eseb_07_3a"
	q_mc_id = "eseb_07_3a_bs500f_got_s2m"
	gset_id = "eseb_wt3a3b"
}
proj_knn_line = function(ref_mat_id, ref_mc_id, q_mat_id, q_mc_id, gset_id)
{
	rmat = scdb_mat(ref_mat_id)
	qmat = scdb_mat(q_mat_id)
	qmc = scdb_mc(q_mc_id)

	qmat_dsamp_n = quantile(colSums(qmat@mat),0.1)
	rmat_dsamp_n = quantile(colSums(rmat@mat),0.1)

	dsamp_n = floor(min(rmat_dsamp_n, qmat_dsamp_n))

	q_feats = gset_get_feat_mat(gset_id=gset_id, mat_id=q_mat_id, 
							downsamp =T, add_non_dsamp=T,
							downsample_n = dsamp_n)

	a_feats = gset_get_feat_mat(gset_id=gset_id, mat_id=ref_mat_id,
							downsamp = T, add_non_dsamp=T,
							downsample_n = dsamp_n)

	k_nonz_exp = 7
	K=30

	gr = tgs_cor_knn(log2(1+k_nonz_exp*as.matrix(q_feats)),
									log2(1+k_nonz_exp*as.matrix(a_feats)), K)

	q_line = qmat@cell_metadata[colnames(qmat@mat),"line"]
	r_line = rmat@cell_metadata[colnames(rmat@mat),"line"]

	gr$col1 = as.character(gr$col1)
	gr$col2 = as.character(gr$col2)

	f_j1 = q_line[gr$col1]=="J1" & r_line[gr$col2] == "J1"
	f_n15 = q_line[gr$col1]=="N15" & r_line[gr$col2] == "N15"

	match_j1 = gr[f_j1,]
	match_j1 = match_j1[!duplicated(match_j1$col1),]
	match_j1 = match_j1[match_j1$col1 %in% names(qmc@mc),]
	match_n15 = gr[f_n15,]
	match_n15 = match_n15[!duplicated(match_n15$col1),]
	match_n15 = match_n15[match_n15$col1 %in% names(qmc@mc),]

	tot_j1_mc = t(tgs_matrix_tapply(rmat@mat[,match_j1$col2], qmc@mc[match_j1$col1], sum))
	rownames(tot_j1_mc) = rownames(rmat@mat)
	tot_n15_mc = t(tgs_matrix_tapply(rmat@mat[,match_n15$col2], qmc@mc[match_n15$col1], sum))
	rownames(tot_n15_mc) = rownames(rmat@mat)

	f = colSums(tot_j1_mc) > 1e+5
	j1_p_egc = t(t(tot_j1_mc[,f])/(colSums(tot_j1_mc)[f]))
	j1_p_legc = log2(j1_p_egc + 1e-5)

	f = colSums(tot_n15_mc) > 1e+5
	n15_p_egc = t(t(tot_n15_mc[,f])/(colSums(tot_n15_mc)[f]))
	n15_p_legc = log2(n15_p_egc + 1e-5)

	tot_j1_q = t(tgs_matrix_tapply(qmat@mat[,match_j1$col1], 
												qmc@mc[match_j1$col1], sum))
	rownames(tot_j1_q) = rownames(qmat@mat)
	f = colSums(tot_j1_q) > 1e+5
	j1_egc = t(t(tot_j1_q[,f])/(colSums(tot_j1_q)[f]))
	j1_legc = log2(j1_egc + 1e-5)

	tot_n15_q = t(tgs_matrix_tapply(qmat@mat[,match_n15$col1], 
												qmc@mc[match_n15$col1], sum))
	rownames(tot_n15_q) = rownames(qmat@mat)
	f = colSums(tot_n15_q) > 1e+5
	n15_egc = t(t(tot_n15_q[,f])/(colSums(tot_n15_q)[f]))
	n15_legc = log2(n15_egc + 1e-5)

	return(list(j1_legc = j1_legc, n15_legc = n15_legc, j1_p_legc = j1_p_legc, n15_p_legc = n15_p_legc))
#for each q cell
}

plt_cmp = function(gnm, res1, res2, qmc1, qmc2) {
	layout(matrix(1:4,nrow=2, byrow=T));
	n15_mcs = intersect(colnames(res1$n15_p_legc),colnames(res1$n15_legc))
	j1_mcs = intersect(colnames(res1$j1_p_legc),colnames(res1$j1_legc))
	plot(res1$n15_legc[gnm, n15_mcs], res1$n15_p_legc[gnm,n15_mcs], pch=19, col=qmc1@colors[n15_mcs]); 
	abline(a=0,b=1)
	plot(res1$j1_legc[gnm, j1_mcs], res1$j1_p_legc[gnm,j1_mcs], pch=19, col=qmc1@colors[j1_mcs]); 
	abline(a=0,b=1)

	n15_mcs = intersect(colnames(res2$n15_p_legc),colnames(res2$n15_legc))
	j1_mcs = intersect(colnames(res2$j1_p_legc),colnames(res2$j1_legc))
	plot(res2$n15_legc[gnm, n15_mcs], res2$n15_p_legc[gnm,n15_mcs], pch=19, col=qmc2@colors[n15_mcs]); 
	abline(a=0,b=1)
	plot(res2$j1_legc[gnm, j1_mcs], res2$j1_p_legc[gnm,j1_mcs], pch=19, col=qmc2@colors[j1_mcs]); 
	abline(a=0,b=1)
}
bxplt_cmp = function(gnm, res1, res2, qmc1, qmc2) {
	layout(matrix(1:4,nrow=2, byrow=T));
	n15_mcs = intersect(colnames(res1$n15_p_legc),colnames(res1$n15_legc))
	j1_mcs = intersect(colnames(res1$j1_p_legc),colnames(res1$j1_legc))

	n15_splt = split(res1$n15_legc[gnm, n15_mcs] - res1$n15_p_legc[gnm,n15_mcs], 
											qmc1@colors[n15_mcs]); 
	boxplot(n15_splt, col=names(n15_splt))
	j1_splt = split(res1$j1_legc[gnm, j1_mcs] - res1$j1_p_legc[gnm,j1_mcs], 
											qmc1@colors[j1_mcs]); 
	boxplot(j1_splt, col=names(j1_splt))

	n15_mcs = intersect(colnames(res2$n15_p_legc),colnames(res2$n15_legc))
	j1_mcs = intersect(colnames(res2$j1_p_legc),colnames(res2$j1_legc))

	n15_splt = split(res2$n15_legc[gnm, n15_mcs] - res2$n15_p_legc[gnm,n15_mcs], 
											qmc1@colors[n15_mcs]); 
	boxplot(n15_splt, col=names(n15_splt))
	j1_splt = split(res2$j1_legc[gnm, j1_mcs] - res2$j1_p_legc[gnm,j1_mcs], 
											qmc1@colors[j1_mcs]); 
	boxplot(j1_splt, col=names(j1_splt))
}
if(0) {
	n15_mcs_a = intersect(colnames(res1$n15_p_legc),colnames(res1$n15_legc))
	j1_mcs_a = intersect(colnames(res1$j1_p_legc),colnames(res1$j1_legc))
	n15_mcs_b = intersect(colnames(res2$n15_p_legc),colnames(res2$n15_legc))
	j1_mcs_b = intersect(colnames(res2$j1_p_legc),colnames(res2$j1_legc))
	d3a_j1 = res1$j1_legc[, j1_mcs_a] - 
						res1$j1_p_legc[rownames(res1$j1_legc),j1_mcs_a]
	d3b_j1 = res2$j1_legc[, j1_mcs_b] - res2$j1_p_legc[rownames(res2$j1_legc),j1_mcs_b]
	d3a_n15 = res1$n15_legc[, n15_mcs_a] - 
					res1$n15_p_legc[rownames(res1$n15_legc),n15_mcs_a]
	d3b_n15 = res2$n15_legc[, n15_mcs_b] - res2$n15_p_legc[rownames(res2$n15_legc),n15_mcs_b]
}
