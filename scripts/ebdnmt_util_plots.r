
ebdnmt_plot_kit = function(mc_id, mat_id, mc2d_id, base_dir, do_gene_maps=T)
{
	atlas_key = read.table("config/atlas_type_order.txt", h=T, stringsAsFactors=F)
	color_ord = 1:nrow(atlas_key)
	names(color_ord) = atlas_key$colour

	mat = scdb_mat(mat_id)
	mc = scdb_mc(mc_id)

	if(is.null(mc)) {
		stop("not mc id ", mc_id)
	}

	mc_day = table(mc@mc, mat@cell_metadata[names(mc@mc), "EB_day"])
	mc_day_max = apply(mc_day,1, function(x) sum(x*0:7)/sum(x))
	mc_ord = order(color_ord[mc@colors] + 0.01*mc_day_max)


	mcell_mc_plot_marks(mc_id, "eseb_key_tfs", 
			mat=mat_id, plot_cells=F, 
			mc_ord = mc_ord,
			fig_fn = sprintf("%s/tfs_heat.png", base_dir),
			reorder_marks=F, fold_burn=1.5)

#plot 2D without ESCs and singletons
	f_mc = mc@colors!="antiquewhite3" & mc@colors!="antiquewhite1"

	mgr = scdb_mc2d(mc2d_id)
	f_gr = f_mc[mgr@graph$mc1] & f_mc[mgr@graph$mc2] & mgr@graph$mc1!=mgr@graph$mc2
	deg = tabulate(c(mgr@graph$mc1[f_gr], mgr@graph$mc2[f_gr]), 
														nbins=length(mc@colors))

	f_mc = f_mc & deg!=0

	tgconfig::set_param("mcell_mc2d_cex",0.7,"metacell")

	mcell_mc2d_plot(mc2d_id, show_mcid=F,cell_outline=T, 
					sc_cex=0.8,edge_w=0.4, 
					fig_fn = sprintf("%s/proj_2d.png", base_dir), 
					filt_mc=f_mc)

	mcell_mc2d_plot_by_factor(mc2d_id = mc2d_id, mat_id = mat_id, meta_field="EB_day", single_plot=F, neto_points=T, filt_mc = f_mc, base_dir = base_dir)

# plot type time trace

	plot_type_trace(mat_id, mc_id, fig_fn = sprintf("%s/type_trace.png", base_dir))

# plot gene 2D projections

	if(do_gene_maps) {
		map_dir = sprintf("%s/genemaps", base_dir)
		if(!dir.exists(map_dir)) {
			dir.create(map_dir)
		}
		plot_gene_maps(mat_id = mat_id,
					mc_id = mc_id,
					mc2d_id = mc2d_id,
					base_dir = map_dir,
					filt_mc = f_mc)
	}
}

plot_gene_maps = function(mat_id, mc_id, mc2d_id, base_dir, filt_mc = NULL)
{
	tgconfig::set_param("mcell_mc2d_gene_shades", colorRampPalette(c("darkblue", "blue", "lightblue", "white", "lightgoldenrod1","orange1", "orange4"))(1000), "metacell")
	tgconfig::set_param("mcell_mc2d_gene_cell_cex",0.5,"metacell")
	tgconfig::set_param("mcell_mc2d_gene_mc_cex",2.5,"metacell")
	tgconfig::set_param("mcell_mc2d_gene_width",500,"metacell")
	tgconfig::set_param("mcell_mc2d_gene_height",500,"metacell")

	mat = scdb_mat(mat_id)
	mc = scdb_mc(mc_id)

	gs = read.table("data/ebdnmt_focus_genes",h=F, stringsAsFactors=F)

	message("will downsamp")
	mat_ds = scm_downsamp(mat@mat, n = round(quantile(colSums(mat@mat),0.1)))

	g_nms = intersect(gs$V1, rownames(mat_ds))
	g_nms = intersect(g_nms, rownames(mc@mc_fp))

	for(g in g_nms) {
		message("plot ", g)
		mcell_mc2d_plot_gene(mc2d_id=mc2d_id, base_dir = base_dir, g,
				max_lfp=2, min_lfp=-2, mat_ds=mat_ds, color_cells=T, neto_points=T,
				filt_mc = filt_mc)
	}
}

plot_type_trace = function(mat_id, mc_id, rep_i = NA, fig_fn=NULL)
{
	color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
	mc = scdb_mc(mc_id)
	md = scdb_mat(mat_id)@cell_metadata
	day = as.numeric(as.character(md[names(mc@mc),"EB_day"]))
	rep = as.numeric(as.character(md[names(mc@mc),"rep"]))
	names(day) = names(mc@mc)
	names(rep) = names(mc@mc)
	if(!is.na(rep_i)) {
		f = md[names(mc@mc), "rep"]==rep_i
		type_day = table(mc@colors[mc@mc[f]], day[names(mc@mc[f])])
		col_ord = intersect(color_ord$colour, mc@colors[f])
	} else {
		type_day = table(mc@colors[mc@mc], day[names(mc@mc)])
		col_ord = intersect(color_ord$colour, mc@colors)
	}
	samp_days = colnames(type_day)
	type_day_n = t(type_day)/colSums(type_day)

	poly_y = apply(type_day_n[,col_ord],1, cumsum)
	poly_y = rbind(rep(0,length(samp_days)), poly_y)
	if(is.null(fig_fn)) {
		fig_fn = sprintf("figs/day_type_%s.png", mc_id)
	}
	png(fig_fn, w=600, h=400)
	plot(NA, xlim=c(0,7), ylim=c(0,1))
	for(i in 2:nrow(poly_y)) {
		polygon(x=c(samp_days,rev(samp_days)), 
					y=c(poly_y[i-1,],rev(poly_y[i,])),
					col=rownames(poly_y)[i])
	}
	dev.off()

}
