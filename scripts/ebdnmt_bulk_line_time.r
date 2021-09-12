

#mat, mc, focus_colors

bulk_time_trace = function(mat_id, mc_id, foc_colors, 
							min_day, max_day, 
							ignore_line = F, ignore_rep = F, layer = NULL,
							min_umis = 100000)
{
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("cannot find mat id ", mat_id)
	}
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("cannot find mc id ", mc_id)
	}
	md = mat@cell_metadata
	if(is.null(foc_colors)) {
		foc_colors = unique(mc@colors)
	}

	cells = names(mc@mc)[mc@colors[mc@mc] %in% foc_colors & md[names(mc@mc),"line"]!="N15S"]
	
	cell_day = as.numeric(as.character( md[cells, "EB_day"]))
	if(ignore_line) {
		cell_line = md[cells, "perturb"]
	} else {
		if(ignore_rep) {
			cell_line = paste(md[cells, "perturb"],md[cells, "line"],sep=".")
		} else {
			cell_line = paste(md[cells, "perturb"],md[cells, "line"],sep=".")
			cell_line = paste(cell_line, md[cells, "rep"], sep=".")
		}
	}
	if(!is.null(layer)) {
		cell_line = paste(cell_line, layer[cells])
	}

	f_cells =  cell_day >= min_day & cell_day <= max_day
	umis = mat@mat[, cells[f_cells]]

	
	bulks = t(tgs_matrix_tapply(umis, paste(cell_day[f_cells], cell_line[f_cells],sep="."), sum))
	rownames(bulks) = rownames(umis)

	bulks_n = t(t(bulks)/colSums(bulks))
	bulks_n = bulks_n[,colSums(bulks)>min_umis]

	return(bulks_n)
}

