
library("devtools")
library("ggplot2"); 
library("ggrepel")
library("qvalue")

load_all("metacell")
scdb_init("scrna_db", force_reinit=T)

cmp_exclude_genes = c("Dnmt3a", "Dnmt3b", "AK145379;H19", "Igf2","Meg3;Mir1906-1;Mir770", "Grb10", "Nnat;Peg5")

df_gnames_to_correct = fread("config/gname_to_correct.txt", sep="\t")
gnames_to_correct = df_gnames_to_correct$short_name
names(gnames_to_correct) = df_gnames_to_correct$long_name

load_data = function(do_dsamp_3ab = T, do_dsamp_dtko=T, do_dsamp_iv =T)
{
	mc_wt = scdb_mc("eseb_07_wt_bs500f_got_s2m")
	mc_3a = scdb_mc("eseb_07_3a_bs500f_got_s2m")
	mc_3b = scdb_mc("eseb_07_3b_bs500f_got_s2m")
	mc_dko = scdb_mc("eseb_07_3ab_bs500f_got_s2m")
	mc_tko = scdb_mc("eseb_07_TKO_bs500f_got_s2m")
	mc_iv = scdb_mc("sing_emb_wt10_recolored")

	mat_wt = scdb_mat("eseb_07_wt")
	mat_3a = scdb_mat("eseb_07_3a")
	mat_3b = scdb_mat("eseb_07_3b")
	mat_dko = scdb_mat("eseb_07_3ab")
	mat_tko = scdb_mat("eseb_07_TKO")
	mat_iv = scdb_mat("sing_emb_wt10")

	legc_wt = log2(mc_wt@e_gc+1e-5)
	legc_3a = log2(mc_3a@e_gc+1e-5)
	legc_3b = log2(mc_3b@e_gc+1e-5)
	legc_dko = log2(mc_dko@e_gc+1e-5)
	legc_tko = log2(mc_tko@e_gc+1e-5)
	legc_iv = log2(mc_iv@e_gc+1e-5)

	tot_wt = colSums(mat_wt@mat)
	tot_3a = colSums(mat_3a@mat)
	tot_3b = colSums(mat_3b@mat)
	tot_dko = colSums(mat_dko@mat)
	tot_tko = colSums(mat_tko@mat)
	tot_iv = colSums(mat_iv@mat)

	md_wt = mat_wt@cell_metadata[names(mc_wt@mc),]
	md_3a = mat_3a@cell_metadata[names(mc_3a@mc),]
	md_3b = mat_3b@cell_metadata[names(mc_3b@mc),]
	md_dko = mat_dko@cell_metadata[names(mc_dko@mc),]
	md_tko = mat_tko@cell_metadata[names(mc_tko@mc),]
	md_iv = mat_iv@cell_metadata[names(mc_tko@mc),]

	day_line_wt = paste(md_wt$EB_day, md_wt$line, sep=".")
	day_line_3a = paste(md_3a$EB_day, md_3a$line, sep=".")
	day_line_3b = paste(md_3b$EB_day, md_3b$line, sep=".")
	day_line_dko = paste(md_dko$EB_day, md_dko$line, sep=".")
	day_line_tko = paste(md_tko$EB_day, md_tko$line, sep=".")

	mat_ds_wt = scm_downsamp(mat_wt@mat, 3000)
	if(do_dsamp_dtko) {
		mat_ds_dko = scm_downsamp(mat_dko@mat, 3000)
		mat_ds_tko = scm_downsamp(mat_tko@mat, 3000)
	} else {
		mat_ds_dko = NULL
		mat_ds_tko = NULL
	}
	if(do_dsamp_3ab) {
		mat_ds_3a = scm_downsamp(mat_3a@mat, 3000)
		mat_ds_3b = scm_downsamp(mat_3b@mat, 3000)
	} else {
		mat_ds_3a= NULL
		mat_ds_3b= NULL
	}
	if(do_dsamp_iv) {
		mat_ds_iv = scm_downsamp(mat_iv@mat, 3000)
	} else {
		mat_ds_iv = NULL
	}

	return(list(mc_wt = mc_wt, mc_3a = mc_3a, mc_3b = mc_3b, mc_dko = mc_dko, mc_tko = mc_tko, mc_iv = mc_iv,
				   mat_wt = mat_wt, mat_3a = mat_3a, mat_3b = mat_3b, mat_dko = mat_dko, mat_tko = mat_tko, mat_iv = mat_iv,
				   mat_ds_wt = mat_ds_wt, mat_ds_3a = mat_ds_3a, mat_ds_3b = mat_ds_3b, mat_ds_dko = mat_ds_dko, mat_ds_tko = mat_ds_tko, mat_ds_iv = mat_ds_iv,
				   md_wt = md_wt, md_3a = md_3a, md_3b = md_3b, md_dko = md_dko, md_tko = md_tko, md_iv = md_iv,
				   legc_wt = legc_wt, legc_3a = legc_3a, legc_3b = legc_3b, legc_dko = legc_dko, legc_tko = legc_tko, legc_iv = legc_iv,
				   day_line_wt = day_line_wt, day_line_3a = day_line_3a, day_line_3b = day_line_3b, day_line_dko = day_line_dko, day_line_tko = day_line_tko))
}

gen_gmods = function(d)
{

	mc_day_wt = table(d$mc_wt@mc, as.character(d$md_wt[,"EB_day"]))
	mc_day_wt_n = mc_day_wt/rowSums(mc_day_wt)

	f_not_esc = mc_day_wt_n[,1]<0.2

	c_foxa1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Foxa1", f_not_esc]))
	c_mesp1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Mesp1", f_not_esc]))
	c_eomes = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Eomes", f_not_esc]))
	c_utf1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Utf1", f_not_esc]))
	c_hand1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Hand1", f_not_esc]))
	c_foxc1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Foxc1", f_not_esc]))
	c_tal1 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Tal1", f_not_esc]))
	c_dppa3 = sort(apply(d$legc_wt[,f_not_esc], 1, cor, d$legc_wt["Dppa3", f_not_esc]))

	gmods = list(endo = names(tail(c_foxa1,30)), 
					 meso = names(tail(c_mesp1,30)), 
					 ps = names(tail(c_eomes,30)), 
					 epi = names(tail(c_utf1,30)),
					 heme = names(tail(c_tal1,30)),
					 dppa3 = names(tail(c_dppa3,30)),
					 emeso = names(tail(c_hand1,30)),
					 rmeso = names(tail(c_foxc1,30)))
	all_gs = unlist(gmods)
	dups = all_gs[duplicated(all_gs)]
	gmods = lapply(gmods, setdiff, dups)
	return(gmods)
}


plot_gmod_day_score = function(legc, day_line, mc, gs, bw=10, 
							legc_wt = d$legc_wt, day_line_wt = d$day_line_wt, 
							mc_wt = d$mc_wt, 
							thresh = NA,
							max_y=0.04, max_day=7, 
							line_nms)
{
	mc_dl = table(mc@mc, as.character(day_line))
	mc_dl_wt = table(mc_wt@mc, as.character(day_line_wt))
	n_miss = length(setdiff(gs,rownames(legc)))
		
	tot_gm = colSums(legc[intersect(gs,rownames(legc)),])
	tot_gm_wt = colSums(legc_wt[intersect(gs,rownames(legc_wt)),])
	if(is.na(thresh)) {
		min_gm_on_type = tapply(tot_gm_wt, mc_wt@colors, min)
		thresh = max(min_gm_on_type)
	}
	xlim = c(min(tot_gm_wt),max(tot_gm_wt))
	if(n_miss > 0) {
		tot_gm = tot_gm + n_miss*log2(1e-5)
	}
	n_l = length(line_nms)
	layout(matrix(1:(max_day*n_l), nrow=n_l, byrow=T))
	lcol = c("darkblue","cyan")
	cnt= 1
	hits=list()
	for(ln_nm in line_nms) {
		wt_ln_nm = ifelse(cnt==1, "J1$", "N15$")
		par(mar=c(ifelse(cnt==1,0,3),2,ifelse(cnt==1,2,0),0))
		frac_hits = c()
		frac_wt_hits = c()
		for(t in 1:max_day) {
			fracs = "f="
			mut_i = grep(sprintf("^%d.%s",t, ln_nm), colnames(mc_dl))
			wt_i = grep(sprintf("^%d.%s",t, wt_ln_nm), colnames(mc_dl_wt))
			if(length(mut_i) == 0 | length(wt_i) == 0) {
				frac_hits[t] = NA
				if(length(wt_i)!=0) {
					x_wt = mean(rep(tot_gm_wt, times=mc_dl_wt[, wt_i])>thresh)
					frac_wt_hits[t] = x_wt
				} else {
					frac_wt_hits[t] = 0
				}
				plot.new()
			} else {
#				message("plot ", ln_nm, " vs ", wt_ln_nm, " at ", mut_i, " and ", wt_i)
				x = mean(rep(tot_gm, times=mc_dl[, mut_i])>thresh)
				x_wt = mean(rep(tot_gm_wt, times=mc_dl_wt[, wt_i])>thresh)
				fracs = paste(ln_nm, paste("frac: ", round(x,3), sep=" ")," ")
				frac_hits[t] = x
				frac_wt_hits[t] = x_wt

				plot(density(rep(tot_gm_wt, times=mc_dl_wt[,wt_i]),bw=bw), 
								xlim=xlim, lty=2, lwd=2, col="black", 
								ylim=c(0, max_y), main=ifelse(cnt==1, fracs, NA),
								yaxt='n',xaxt=ifelse(cnt==1,'n', 's'))
				abline(v=thresh, lty=2, col="gray")
				lines(density(rep(tot_gm, times=mc_dl[,mut_i]),bw=bw),col=lcol[cnt], lwd=3)
			}
			if(t == 6) {
				par(mar=c(ifelse(cnt==1,0,3),0,ifelse(cnt==1,2,0),2))
			} else {
				par(mar=c(ifelse(cnt==1,0,3),2,ifelse(cnt==1,2,0),0))
			}
		}
		hits[[cnt]] = matrix(c(frac_wt_hits,frac_hits),nrow=2, byrow=T)
		cnt = cnt+1
	}
	return(hits)
}

plot_all_gmod_dists = function(d, gs)
{
	for(prog in names(gmods)) {
		all_hits = list()
		for(nm in c("dko","tko","3a", "3b")) {
			ln_nms = c("J1$", "N15$")
			if(nm == "dko") {
				ln_nms = c("DKO$", "DKO16$")
			}
			if(nm == "tko") {
				ln_nms = c("J1$")
			}
			nms = grep(nm, names(d), v=T)
			png(sprintf("paper_figs/prog_dist/%s.%s.png", prog, nm),w=700,h=ifelse(nm == "tko", 125, 250))
			hits = plot_gmod_day_score(legc= d[[sprintf("legc_%s",nm)]],
									  day_line = d[[sprintf("day_line_%s",nm)]], 
									  mc = d[[sprintf("mc_%s",nm)]], 
									  gs = gmods[[prog]],
									  thresh = ifelse(prog == "heme", -480, NA),
									  line_nms = ln_nms,
									  bw=5, max_y=0.06)
			dev.off()
			png(sprintf("paper_figs/prog_dist/dist.%s.%s.png", prog, nm),w=400,h=ifelse(nm == "tko", 200, 400))
			layout(matrix(1:length(hits),ncol=1))
			ymax = max(sapply(hits,max,na.rm=T))*1.1
			par(mar=c(2,2,2,2))
			for(i in 1:length(hits)) {
				barplot(hits[[i]], beside=T, 
					col=c("darkgray",ifelse(i==1, "darkblue", "cyan")), 
					ylim=c(0,ymax),xaxt='n')
				par(mar=c(2,2,1,2))
			}
			dev.off()
		}
	}
}
plot_3a3b_gmod_dists = function(d, gs)
{
	ln_nms = c("J1$", "N15$")
	for(prog in names(gmods)) {
		png("/tmp/safta.png")
		hits3a = plot_gmod_day_score(legc= d[["legc_3a"]],
				  day_line = d[["day_line_3a"]], 
				  mc = d[["mc_3a"]], 
				  gs = gmods[[prog]],
				  line_nms = ln_nms,
				  bw=5, max_y=0.06)
		hits3b = plot_gmod_day_score(legc= d[["legc_3b"]],
				  day_line = d[["day_line_3b"]], 
				  mc = d[["mc_3b"]], 
				  gs = gmods[[prog]],
				  line_nms = ln_nms,
				  bw=5, max_y=0.06)
		dev.off()

		ab_dat1 = rbind(hits3a[[1]][1,], hits3a[[1]][2,], hits3b[[1]][2,])
		ab_dat2 = rbind(hits3a[[2]][1,], hits3a[[2]][2,], hits3b[[2]][2,])
		png(sprintf("paper_figs/prog_dist/dist%s.ab.png", prog),w=400,h=400)
		layout(matrix(1:2,ncol=1))
		ymax = max(max(ab_dat1,na.rm=T), max(ab_dat2,na.rm=T))
		par(mar=c(2,2,2,2))
		barplot(ab_dat1, beside=T, 
				col=c("darkgray", "darkblue","darkred"),
				ylim=c(0,ymax),xaxt='n')
		par(mar=c(2,2,1,2))
		barplot(ab_dat2, beside=T, 
				col=c("darkgray", "darkblue","darkred"),
				ylim=c(0,ymax),xaxt='n')
		dev.off()
	}
}


compare_top_gmod_cells = function(d, gs, mut_nm, 
							back_nm = "wt", gs_back = gs,  K_cells = 50,
							focus_line = NA, focus_line_back = NA,
							exclude_genes=cmp_exclude_genes)
{
	mut_mat_ds = d[[sprintf("mat_ds_%s",mut_nm)]] 
	back_mat_ds = d[[sprintf("mat_ds_%s",back_nm)]] 

	if(!is.na(focus_line)) {
		mut_line = d[[sprintf("mat_%s",mut_nm)]]@cell_metadata[colnames(mut_mat_ds),"line"]
		mut_cells = colnames(mut_mat_ds)[mut_line==focus_line]
		mut_mat_ds = mut_mat_ds[,mut_cells]
	} else {
		mut_cells = colnames(mut_mat_ds)
	}
	if(!is.na(focus_line_back)) {
		back_line = d[[sprintf("mat_%s",back_nm)]]@cell_metadata[colnames(back_mat_ds),"line"]
		back_cells = colnames(back_mat_ds)[back_line==focus_line_back]
		back_mat_ds = back_mat_ds[,back_cells]
	} else {
		back_cells = colnames(back_mat_ds)
	}

	mut_mat = d[[sprintf("mat_%s",mut_nm)]]@mat[,mut_cells]
	back_mat = d[[sprintf("mat_%s",back_nm)]]@mat[,back_cells]

	gs = intersect(rownames(mut_mat), gs)
	gs_back = intersect(rownames(back_mat), gs_back)

	mut_us = mut_mat_ds[gs,]	
	back_us = back_mat_ds[gs_back,]	
	mut_score = colSums(log2(1+mut_us*3))
	back_score = colSums(log2(1+back_us*3))

	top_mut_cells = mut_score > -sort(-mut_score,partial=K_cells)[K_cells]
	top_back_cells = back_score > -sort(-back_score,partial=K_cells)[K_cells]
	mut_top_e_g = rowSums(mut_mat[,top_mut_cells])
	back_top_e_g = rowSums(back_mat[,top_back_cells])
	all_gs = intersect(names(mut_top_e_g), names(back_top_e_g))
	a = data.frame(legc_mut = log2(1e-5+mut_top_e_g[all_gs]/sum(mut_top_e_g)),
							legc_back = log2(1e-5+back_top_e_g[all_gs]/sum(back_top_e_g)),
							mut = mut_top_e_g[all_gs],
							back = back_top_e_g[all_gs])

	N = min(sum(a$mut), sum(a$back))
	tot1 = a$mut*N/sum(a$mut)
	names(tot1) = rownames(a)
	tot2 = round(a$back*N/sum(a$back))
	names(tot2) = rownames(a)
	exclude_genes = intersect(exclude_genes, names(tot1))
	tot1[exclude_genes]=0
	tot2[exclude_genes]=0
	df = data.frame(tot1=tot1, tot2=tot2, 
									le1=log2(3e-5+tot1/sum(tot1)),
									le2=log2(3e-5+tot2/sum(tot2)))
	df$mean = (df$le1+df$le2)/2
	df$diff= (df$le1-df$le2)

	df$lab = names(tot1)
	trim_label = df$lab
	names(trim_label) = trim_label
	trim_label[names(gnames_to_correct)] = gnames_to_correct
	df$lab = trim_label

	df = df[(df$tot1+df$tot2)>15,]
	ns = c(sum(tot1), sum(tot2))
	df$pv = apply(cbind(df$tot1,df$tot2),1,function(x) { mat=matrix(c(x,ns),nrow=2); chisq.test(mat)$p.value }) 
	df$qv = qvalue(df$pv)$qvalues
	return(df)
}

plot_cmp_top_gmods = function(d, gmods, n_genes=40, exclude_genes=cmp_exclude_genes, text_size=14, force=30, w=1000, h=1000, K_cells=100)
{
	for(nm in names(gmods)) {
		for(mut in c("dko", "tko", "3a","3b", "iv")) {
			message("plot ", nm, " mut ", mut)
			a = compare_top_gmod_cells(d, gmods[[nm]], mut, exclude_genes=exclude_genes, K_cells=K_cells)
			plt_cmp_gmod(a, n_genes, nm=nm, mut=mut, text_size= text_size, force=force,  w=w, h=h)
		}
	}
}
plt_cmp_gmod = function(df, n_genes, nm, mut, text_size = 14, force = 30, w=1000, h=1000)
{	
	df$rank = rank(df$pv)
	df$col = ifelse(df$diff < -log2(1.5) & df$pv < 1e-2 & df$rank<n_genes,"darkblue","gray")
	df$col = ifelse(df$diff > log2(1.5) & df$pv < 1e-2 & df$rank<n_genes,"darkred",df$col)

	lab_theme = element_text(size=32)
	p = ggplot(df, aes(le1, le2)) + 
				geom_point(col=df$col, size=text_size/2) +
				labs(x=sprintf("%s log2(%s RNA perc.)", nm, mut)) +
				labs(y=sprintf("%s log2(wt RNA perc.)", nm, mut)) +
				theme(axis.text.x = lab_theme, axis.text.y = lab_theme,
						axis.title.x = lab_theme, axis.title.y = lab_theme) +
				geom_text_repel(data = filter(df, df$col!="gray"), 
							size=text_size,force=force, max.iter=1e+4,aes(label=lab))
			
	png(sprintf("paper_figs/cmp_gmods/%s.%s.png", nm, mut),w=w,h=h)
	print(p)
	dev.off()
}

plt_cpm_3a3b_lines = function(d, gmods, prog, K_cells=50, text_size=12, force=30, w=1000,h=1000, exclude_genes=cmp_exclude_genes) 
{
	gs = gmods[[prog]]
	dfj1 = compare_top_gmod_cells(d, gs, mut_nm="3a", back_nm="3b", focus_line="J1", focus_line_back="J1", K_cells=K_cells, exclude_genes=exclude_genes)
	dfn15 = compare_top_gmod_cells(d, gs, mut_nm="3a", back_nm="3b", focus_line="N15", focus_line_back="N15", K_cells=K_cells, exclude_genes=exclude_genes)

	ns = union((dfj1$lab)[dfj1$qv<0.05], (dfn15$lab)[dfn15$qv<0.05])

	df = data.frame(line1 = dfj1[ns,"diff"], line2 = dfn15[ns,"diff"], lab=ns)

	ctest = cor.test(df$line1, df$line2)

	df$min_qv = log10(pmax(dfj1[ns,"qv"], dfn15[ns,"qv"]))
	df$min_qv = ifelse(df$line1 > 0 & df$line2 > 0, -df$min_qv, df$min_qv)
	df$min_qv[is.na(df$min_qv)]=0
	df$min_qv = pmax(pmin(df$min_qv, 4),-4)
	message(paste(quantile(df$min_qv),collapse=" "))

	lab_theme = element_text(size=32)
	p = ggplot(df, aes(line1, line2)) + 
				geom_point(aes(fill=min_qv), size=text_size/2,pch=21) +
				scale_fill_gradientn(colors=c("darkblue","blue","lightgray","lightgray", "red", "darkred"),
									breaks=c(-4,-1,-0.8,0.8,1,4)) +
				labs(x=sprintf("%s J1 3a-3b log2(RNA perc)", prog)) +
				labs(y=sprintf("%s N15 3a-3b log2(RNA perc)", prog)) +
				theme(axis.text.x = lab_theme, axis.text.y = lab_theme,
						axis.title.x = lab_theme, axis.title.y = lab_theme) +
				geom_text_repel(data = 
							filter(df, abs(df$line1+df$line2)>1 & abs(df$min_qv) > 1), 
							size=text_size, force=force, max.iter=1e+4,aes(label=lab))
	p = p + xlim(-2,2) + ylim(-2, 2)
	p = p + ggtitle(sprintf("%s repli cor %.3f pv %e", prog, ctest$estimate, ctest$p.value))
	png(sprintf("paper_figs/cmp_gmods/lines3a3b.%s.png", prog),w=w,h=h)
	print(p)
	dev.off()
}

plot_g_gmod_compare = function(d, gmods, gmod, gene, legc1, legc2, col1, col2, base_dir, suf="")
{
	f1 = legc1["Klf4",] < -13
	f2 = legc2["Klf4",] < -13

	gmod_gs = gmods[[gmod]]
	gmod_gs = setdiff(gmod_gs, gene)

	x1 = colSums(legc1[gmod_gs,f1])
	x2 = colSums(legc2[gmod_gs,f2])
	xlim = c(min(c(x1,x2)), max(c(x1,x2)))
	ylim = c(min(c(legc1[gene,f1],legc2[gene,f2])),
			   max(c(legc1[gene,f1],legc2[gene,f2])))
	
   lo = loess(legc1[gene,f1] ~ x1,span=0.3)$fitted
	png(sprintf("%s/%s_%s%s.png", base_dir, gmod, gene, suf),w=500,h=1000)
	layout(matrix(1:2,nrow=2))
	par(mar=c(0,2,2,2))
	plot(x1, legc1[gene,f1], pch=21,bg=col1[f1], xlab=gmod, ylab=gene, ylim=ylim,xlim=xlim, xaxt='n', main=gene, cex=1.6)
	lines(sort(x1), lo[order(x1)],col="darkred",lwd=4,lty=2)
	par(mar=c(2,2,0,2))
	plot(x2, legc2[gene,f2], pch=21,bg=col2[f2], xlab=gmod, ylab=gene, ylim=ylim, xlim=xlim, cex=1.6)
	lines(sort(x1), lo[order(x1)],col="darkred",lwd=4,lty=2)
	dev.off()
}

comp_egc_on_lines = function(d, type="wt", line1, line2)
{
	mat = d[[sprintf("mat_%s", type)]]
	md = d[[sprintf("md_%s", type)]]
	mc = d[[sprintf("mc_%s", type)]]

	f_line1 = md[names(mc@mc), "line"] == line1
	f_line2 = md[names(mc@mc), "line"] == line2
	cells1 = names(mc@mc)[f_line1]
	cells2 = names(mc@mc)[f_line2]

	egc1 = t(tgs_matrix_tapply(mat@mat[,cells1], mc@mc[cells1], sum))
	rownames(egc1) = rownames(mat@mat)
	egc2 = t(tgs_matrix_tapply(mat@mat[,cells2], mc@mc[cells2], sum))
	rownames(egc2) = rownames(mat@mat)

	egc1_n = t(tgs_matrix_tapply(mat@mat[,cells1], mc@mc[cells1], function(y) {exp(mean(log(1+y)))-1}))
	rownames(egc1_n) = rownames(mat@mat)
	mc_meansize1 = tapply(colSums(mat@mat[,cells1]), mc@mc[cells1], mean)
	egc1_n = t(t(egc1_n)/as.vector(mc_meansize1))

	egc2_n = t(tgs_matrix_tapply(mat@mat[,cells2], mc@mc[cells2], function(y) {exp(mean(log(1+y)))-1}))
	rownames(egc2_n) = rownames(mat@mat)
	mc_meansize2 = tapply(colSums(mat@mat[,cells2]), mc@mc[cells2], mean)
	egc2_n = t(t(egc2_n)/as.vector(mc_meansize2))

	return(list(egc1=egc1, egc2=egc2, 
					egc1_n = egc1_n, egc2_n = egc2_n
				))
}


strat_gmod_on_lines = function(d, type="wt", line1, line2, gmod, bin_vals = seq(0,100,10))
{
	mat = d[[sprintf("mat_%s", type)]]
	mat_ds = d[[sprintf("mat_ds_%s", type)]]
	mat = mat@mat[,colnames(mat_ds)]
	md = d[[sprintf("md_%s", type)]]
	mc = d[[sprintf("mc_%s", type)]]

	f_line1 = md[names(mc@mc), "line"] == line1
	f_line2 = md[names(mc@mc), "line"] == line2
	cells1 = intersect(names(mc@mc)[f_line1], colnames(mat))
	cells2 = intersect(names(mc@mc)[f_line2], colnames(mat))

	gmod_score = colSums(mat_ds[gmod,])
	bin_vals[length(bin_vals)] = max(c(bin_vals, max(gmod_score)))+1
	bin_vals[1] = min(c(bin_vals, min(gmod_score)))-1
	gmod_bin = as.numeric(cut(gmod_score, bin_vals))
	names(gmod_bin) = names(gmod_score)

	strati1 = t(tgs_matrix_tapply(mat[,cells1], gmod_bin[cells1], sum))
	rownames(strati1) = rownames(mat)
	strati2 = t(tgs_matrix_tapply(mat[,cells2], gmod_bin[cells2], sum))
	rownames(strati2) = rownames(mat)

	return(list(starti1= strati1, strati2=strati2,
					strati1_n = t(t(strati1)/colSums(strati1)),
					strati2_n = t(t(strati2)/colSums(strati2)),
					gmod_score = gmod_score
				))
}

compare_endomeso_3ab = function(d, gmods)
{
	endo3a = colSums(d$mat_ds_3a[gmods$endo,])	
	endo3b = colSums(d$mat_ds_3b[gmods$endo,])	

	meso3a = colSums(d$mat_ds_3a[gmods$meso,])	
	meso3b = colSums(d$mat_ds_3b[gmods$meso,])	
	emeso3a = colSums(d$mat_ds_3a[gmods$emeso,])	
	emeso3b = colSums(d$mat_ds_3b[gmods$emeso,])	
	rmeso3a = colSums(d$mat_ds_3a[gmods$rmeso,])	
	rmeso3b = colSums(d$mat_ds_3b[gmods$rmeso,])	
	epi3a = colSums(d$mat_ds_3a[gmods$epi,])	
	epi3b = colSums(d$mat_ds_3b[gmods$epi,])	

	rep_3a = as.character(d$md_3a[names(endo3a),"rep"])
	rep_3b = as.character(d$md_3b[names(endo3b),"rep"])
	day_3a = as.character(d$md_3a[names(endo3a),"EB_day"])
	day_3b = as.character(d$md_3b[names(endo3b),"EB_day"])
	line_3a = d$md_3a[names(endo3a),"line"]
	line_3b = d$md_3b[names(endo3b),"line"]

	rep_line_3a = paste(rep_3a, line_3a,sep="_")
	rep_line_3b = paste(rep_3b, line_3b,sep="_")

	plt = function(line, day, r1, r2, v3a, v3b, pref) {
		png(sprintf("paper_figs/fig3/%s_score_3ab_d%d_%s.png", pref,day,line), w=900,h=300)
		layout(matrix(1:length(r1),nrow=1))
		for(i in 1:length(r1)) {
			rl1 = paste(r1[i],line,sep="_")
			rl2 = paste(r2[i],line,sep="_")
			ks = ks.test(v3a[day_3a==day & rep_line_3a==rl1],v3b[day_3b==day & rep_line_3b==rl2])
			D = ks$statistic
			pv = ks$p.value
			plot(ecdf(log2(1+v3a[day_3a==day & rep_line_3a==rl1])), 
									verticals=T, do.points=F, lwd=3,
									main=sprintf("d%s, %s - %s %.2f %e", day,rl1,rl2,D,pv))
			lines(ecdf(log2(1+v3b[day_3b==day & rep_line_3b==rl2])), 
											col="red", verticals=T, do.points=F, lwd=3)
		}
		dev.off()
	}

#J1: d4 4 and 5, d5 only 4
	plt("J1", 5, c(2,3,4), c(1,3,4), endo3a, endo3b, "endo")
	plt("J1", 4, c(2,3,4), c(2,3,4), endo3a, endo3b, "endo")
	plt("J1", 5, c(2,3,4), c(1,3,4), meso3a, meso3b, "meso")
	plt("J1", 4, c(2,3,4), c(2,3,4), meso3a, meso3b, "meso")
	plt("N15", 5, c(1,2,3), c(1,2,3), endo3a, endo3b, "endo")
	plt("N15", 4, c(1,2,3), c(1,2,3), endo3a, endo3b, "endo")
	plt("N15", 5, c(1,2,3), c(1,2,3), meso3a, meso3b, "meso")
	plt("N15", 4, c(1,2,3), c(1,2,3), meso3a, meso3b, "meso")
	plt("J1", 5, c(2,3,4), c(1,3,4), meso3a+emeso3a+rmeso3a, meso3b+emeso3b+rmeso3b, "allmeso")
	plt("J1", 4, c(2,3,4), c(2,3,4), meso3a+emeso3a+rmeso3a, meso3b+emeso3b+rmeso3b, "allmeso")
	plt("N15", 5, c(1,2,3), c(1,2,3), meso3a+emeso3a+rmeso3a, meso3b+emeso3b+rmeso3b, "allmeso")
	plt("N15", 4, c(1,2,3), c(1,2,3), meso3a+emeso3a+rmeso3a, meso3b+emeso3b+rmeso3b, "allmeso")

	plt("J1", 5, c(2,3,4), c(1,3,4), emeso3a, emeso3b, "emeso")
	plt("J1", 4, c(2,3,4), c(2,3,4), emeso3a, emeso3b, "emeso")
	plt("N15", 5, c(1,2,3), c(1,2,3), emeso3a, emeso3b, "emeso")
	plt("N15", 4, c(1,2,3), c(1,2,3), emeso3a, emeso3b, "emeso")
	plt("J1", 5, c(2,3,4), c(1,3,4), rmeso3a, rmeso3b, "rmeso")
	plt("J1", 4, c(2,3,4), c(2,3,4), rmeso3a, rmeso3b, "rmeso")
	plt("N15", 5, c(1,2,3), c(1,2,3), rmeso3a, rmeso3b, "rmeso")
	plt("N15", 4, c(1,2,3), c(1,2,3), rmeso3a, rmeso3b, "rmeso")

	plt("J1", 5, c(2,3,4), c(1,3,4), epi3a, epi3b, "epi")
	plt("J1", 4, c(2,3,4), c(2,3,4), epi3a, epi3b, "epi")
	plt("N15", 5, c(1,2,3), c(1,2,3), epi3a, epi3b, "epi")
	plt("N15", 4, c(1,2,3), c(1,2,3), epi3a, epi3b, "epi")

	png("paper_figs/fig3/endo_all_J1.png", w=400,h=400)
	ks = ks.test(endo3a[day_3a==4 & line_3a=="J1"],
							endo3b[day_3b==4 & line_3b=="J1"])
	x = ecdf(log2(1+endo3a[day_3a==4 & line_3a=="J1"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("J1 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+endo3b[day_3b==4 & line_3b=="J1"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

	png("paper_figs/fig3/endo_all_N15.png", w=400,h=400)
	ks = ks.test(endo3a[day_3a==4 & line_3a=="N15"],
							endo3b[day_3b==4 & line_3b=="N15"])
	x = ecdf(log2(1+endo3a[day_3a==4 & line_3a=="N15"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("N15 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+endo3b[day_3b==4 & line_3b=="N15"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

	png("paper_figs/fig3/meso_all_J1.png", w=400,h=400)
	ks = ks.test(meso3a[day_3a==4 & line_3a=="J1"],
							meso3b[day_3b==4 & line_3b=="J1"])
	x = ecdf(log2(1+meso3a[day_3a==4 & line_3a=="J1"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("J1 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+meso3b[day_3b==4 & line_3b=="J1"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

	png("paper_figs/fig3/meso_all_N15.png", w=400,h=400)
	ks = ks.test(meso3a[day_3a==4 & line_3a=="N15"],
							meso3b[day_3b==4 & line_3b=="N15"])
	x = ecdf(log2(1+meso3a[day_3a==4 & line_3a=="N15"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("N15 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+meso3b[day_3b==4 & line_3b=="N15"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

	png("paper_figs/fig3/meso_all5_J1.png", w=400,h=400)
	ks = ks.test(meso3a[day_3a==5 & line_3a=="J1"],
							meso3b[day_3b==5 & line_3b=="J1"])
	x = ecdf(log2(1+meso3a[day_3a==5 & line_3a=="J1"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("J1 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+meso3b[day_3b==5 & line_3b=="J1"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

	png("paper_figs/fig3/meso_all5_N15.png", w=400,h=400)
	ks = ks.test(meso3a[day_3a==5 & line_3a=="N15"],
							meso3b[day_3b==5 & line_3b=="N15"])
	x = ecdf(log2(1+meso3a[day_3a==5 & line_3a=="N15"]))
	plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
						col="darkred", lwd=3, 
						main=sprintf("N15 %.2f %e", ks$statistic, ks$p.value),
						xlim=c(1,6))
	x = ecdf(log2(1+meso3b[day_3b==5 & line_3b=="N15"]))
	lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
	dev.off()

   png("paper_figs/fig3/epi_all_J1.png", w=400,h=400)
   ks = ks.test(epi3a[day_3a==4 & line_3a=="J1"],
                     epi3b[day_3b==4 & line_3b=="J1"])
   x = ecdf(log2(1+epi3a[day_3a==4 & line_3a=="J1"]))
   plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
                  col="darkred", lwd=3,
                  main=sprintf("J1 %.2f %e", ks$statistic, ks$p.value),
                  xlim=c(1,6))
   x = ecdf(log2(1+epi3b[day_3b==4 & line_3b=="J1"]))
   lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
   dev.off()

   png("paper_figs/fig3/epi_all_N15.png", w=400,h=400)
   ks = ks.test(epi3a[day_3a==4 & line_3a=="N15"],
                     epi3b[day_3b==4 & line_3b=="N15"])
   x = ecdf(log2(1+epi3a[day_3a==4 & line_3a=="N15"]))
   plot(seq(1,6,0.01),1-x(seq(1,6,0.01)), type="l",
                  col="darkred", lwd=3,
                  main=sprintf("N15 %.2f %e", ks$statistic, ks$p.value),
                  xlim=c(1,6))
   x = ecdf(log2(1+epi3b[day_3b==4 & line_3b=="N15"]))
   lines(seq(1,6,0.01),1-x(seq(1,6,0.01)), col="darkblue", lwd=3)
   dev.off()
}

plot_gmods_cors = function(d, gmods, lead_g,
						cols=rep("darkblue", length(gmods)))
{
	i = 1
	for(nm in names(gmods)) {
		gcor = tgs_cor(t(d$legc_wt[gmods[[nm]],]))
		png(sprintf("paper_figs/fig1/%s.png", nm), w=250, h=600)
		par(mar=c(3,12,2,1))
		barplot(sort(gcor[lead_g[i],setdiff(gmods[[nm]],lead_g[i])]), horiz=T, las=2, col=cols[i])
		dev.off()
		i = i + 1
	}
}
