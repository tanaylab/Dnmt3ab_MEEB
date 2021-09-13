
gen_sfig1 = function()
{

#QCs on eseb

	mat = scdb_mat("eseb_07_wt")
	md = mat@cell_metadata
	
	bad_genes07 = c(as.character(read.table("data/filtered_genes_eseb07_wt.txt", sep="\t")$x))
	bad_genes07 = union(bad_genes07, c(as.character(read.table("data/filter_gene_line_batch_wt.txt", h=T)$x)))

	mat_ds = scm_downsamp(mat@mat, 2500)

	bad_cors = tgs_cor(as.matrix(t(mat_ds[bad_genes07,])))
	diag(bad_cors) = 0
	cmax = apply(bad_cors,1,max)
	bad_cors = bad_cors[cmax > 0.1, cmax > 0.1]

	bad_cors2 = tgs_cor(bad_cors)
	km = tglkmeans::TGL_kmeans(bad_cors2, 10, id_column=F)
	
#	hc = hclust(tgs_dist(bad_cors), "ward.D2")
	shades = colorRampPalette(c("brown", "darkblue", "blue", "white", "red", "darkred", "gold"))(1000)
	
	pheatmap::pheatmap(
			pmin(pmax(bad_cors[order(km$cluster), order(km$cluster)],-0.3),0.3),cluster_cols=F, cluster_rows=F,  color=shades, filename="paper_figs/sfig1/bad_genes_wt_cors.png", width=25,heigh=25, fonsize=3)

	names(km$cluster) = rownames(bad_cors2)
	gm = km$cluster
	gm_tot = tgs_matrix_tapply(t(mat_ds[names(gm),]), gm, sum)

	cnms = colnames(mat_ds)

	write.table(sort(gm), file="paper_figs/sfig1/bat_mods_tab.txt", quote=F)
	for(i in 1:10) {
		png(sprintf("paper_figs/sfig1/batch_mods_wt_%d.png", i), w=1200, h=600)
		par(mar = c(14,4,2,2))
		boxplot(split(gm_tot[i,], paste(md[cnms, "exp_state"], md[cnms,"rep"], sep="_R")),las=2, col="darkblue", ylab=NA, main=paste("module", i),cex.main=2, cex.axis=2 )
		dev.off()
	}
}
