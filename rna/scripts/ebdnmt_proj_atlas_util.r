
build_atlases = function()
{
type_annot = read.table("data/gotg_type_annot.txt")
temp_type_annot = read.table("data/temp_type_annot.txt")

mcell_new_mcatlas("emb_gotg_75", mat_id = "emb_gotg_75", naming_type="tenx", mc_id = "emb_gotg_75_bs500f", graph_id= "emb_gotg_75", gset_id = "emb_gotg_75", type_annot=type_annot)

mcell_new_mcatlas("sing_emb_wt10", mat_id = "sing_emb_wt10", naming_type="mars", mc_id = "sing_emb_wt10_recolored", graph_id= "sing_emb_wt10", gset_id = "sing_emb_wt10", type_annot=temp_type_annot)

}

atlas_annot_w_all = function(mc_id, qmat_id, recomp_knn=T, max_entropy=2) 
{
	message("will annotate ", mc_id, " using all")

	for(atlas_id in c("emb_gotg_75", "sing_emb_wt10")) {
		message("doing atlas ", atlas_id)
		if(atlas_id == "emb_gotg_75") {
			atlas_short = "got"
		} else {
			atlas_short = "tem"
		}
		atlas = scdb_mcatlas(atlas_id)
		atlas@atlas_name = atlas_id

#annots = mcatlas_annotate_sc_by_projection(atlas, "eseb_07_wt", qmat_naming_type="mars", all_vs_all=F)

#t_annots = mcatlas_annotate_sc_by_projection(t_atlas, "eseb_07_wt", qmat_naming_type="mars", all_vs_all=F)


		qmc_id = mc_id
		new_qmc_id = paste(mc_id,atlas_short,sep="_")
		new_qmc_id_mc2mc = paste(new_qmc_id,"m2m",sep="_")
		new_qmc_id_sc2mc = paste(new_qmc_id,"s2m",sep="_")
		new_qmc_id_sc2sc = paste(new_qmc_id,"s2s",sep="_")

		obj_fn = sprintf("scrna_db/annots/%s.%s.Rda", mc_id, atlas_short)

		if(recomp_knn | !file.exists(obj_fn)) {
			message("will create new annots")
	      sc_annots = metacell:::mcatlas_annotate_sc_by_projection(atlas, qmat_id, 
										qmat_naming_type = "mars", all_vs_all = F)
  	    	sc_annots$color = as.character(sc_annots$color)
			save(sc_annots, file=obj_fn)
		} else {
			message("loading annots object")
			load(obj_fn)
		}
		s = mcatlas_annotate_mc_by_sc2sc_projection(
				atlas_id=atlas_id, 
				qmc_id = qmc_id, 
				qmat_id = qmat_id,
				qmat_naming_type = "mars", all_vs_all =F,
				new_qmc_id = new_qmc_id_sc2sc, 
				T_cor_gap=1, fig_cmp_dir=NULL,
				sc_annots = sc_annots)

		mc_annots = mcatlas_annotate_mc_by_mc2mc_projection(
				atlas_id=atlas_id, 
				qmc_id = qmc_id, 
				qmat_naming_type = "mars",
				new_qmc_id = new_qmc_id_mc2mc, 
				T_cor_gap=1, fig_cmp_dir=NULL)

		s2m_proj_egc = mcatlas_annotate_mc_by_sc2mc_projection(
				atlas_id=atlas_id, 
				qmat_id = qmat_id, 
				qmc_id = qmc_id, 
				qmat_naming_type = "mars",
				new_qmc_id = new_qmc_id_sc2mc, 
				max_entropy = max_entropy,
				fig_cmp_dir=NULL)

		save(s2m_proj_egc, mc_annots, sc_annots, file=obj_fn)
	}
}


