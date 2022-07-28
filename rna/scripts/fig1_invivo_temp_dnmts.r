
library(metacell)
scdb_init("scrna_db/", force_reinit=T)

max_t = 13
mat = scdb_mat("sing_emb_wt10")
md = mat@cell_metadata
mc = scdb_mc("sing_emb_wt10_recolored")
cell_t = md$age_group
names(cell_t) = rownames(md)

ctot = colSums(mat@mat)

t_cols = mc@color_key$color
names(t_cols) = mc@color_key$group
t_nms = mc@color_key$group
names(t_nms) = mc@color_key$color

n_mt = table(mc@mc, md[names(mc@mc),"age_group"])
p_mt = t(t(n_mt)/colSums(n_mt))

base_3a = max(log2(1e-5+mc@e_gc["Dnmt3a",]))
base_3b = max(log2(1e-5+mc@e_gc["Dnmt3b",]))

mct = scdb_mctnetwork("sing_emb_wt10")
plot_traceback = function(type, max_t=13, add=T)
{
	col = t_cols[type]
	f_type = mc@colors==col
	p_type = p_mt[,max_t] * f_type
	p_type = p_type/sum(p_type)

	probs = mctnetwork_propogate_from_t(mct, t=max_t, p_type)
	mc_prob = t(t(probs$probs)/colSums(probs$probs))

	prof_a = log2(1e-5+mc@e_gc["Dnmt3a",] %*% mc_prob)
	prof_b = log2(1e-5+mc@e_gc["Dnmt3b",] %*% mc_prob)
	if(add) {
		lines(1:max_t, prof_b-base_3b, lwd=3, col=col, type="l")
	} else {
		plot(1:max_t, prof_b-base_3b, lwd=3, col=col, ylim=c(-5,0), type="l", ylab=NA,xlab=NA, xaxt='n', yaxt='n')
	}
	lines(1:max_t, prof_a-base_3a, lwd=3, col=col,lty=2)
	abline(h=c(-1,-2,-3,-4), lty=2, col="gray")
	abline(v=seq(2,12,l=6), lty=2, col="gray")
}

png("paper_figs/fig1/dnmt_invivo_temporal.png", w=1500,h=300)
layout(matrix(1:5, nrow=1))

plot_traceback("Foregut", add=F)
plot_traceback("Node/Notochord", add=F)
plot_traceback("Amnion/Chorion", add=F)
plot_traceback("ExE mesoderm",add=T)
plot_traceback("Allantois",add=T)
plot_traceback("Erythroid 2",add=T)
plot_traceback("Cardiomyocytes",add=F)
plot_traceback("Rostral mesoderm",add=T)
plot_traceback("Paraxial mesoderm",add=T)
plot_traceback("Caudal mesoderm",add=T)
plot_traceback("Definitive ectoderm",add=F)
plot_traceback("Rostral neural plate",add=T)
plot_traceback("Caudal neural plate",add=T)
plot_traceback("Caudal epiblast",add=T)
dev.off()

mat = scdb_mat("eseb_07_wt")
md = mat@cell_metadata
mc = scdb_mc("eseb_07_wt_bs500f_s2m")
legc = log2(1e-5+mc@e_gc)
mc_t_n = table(mc@mc, as.character(md[names(mc@mc),"EB_day"]))
mc_t_e = apply(mc_t_n,1,function(x) return(sum(x*0:7)/sum(x)))

f_endo = mc@colors=="#EF5A9D" | mc@colors=="#c19f70" | mc@colors=="#F397C0"

f_meso = mc@colors == "#cc7818" | mc@colors=="#B51D8D" | mc@colors=="#C594BF" | mc@colors=="#408DA1" | mc@colors=="#C9EBFB" | mc@colors=="#8DB5CE" |

f_bad = mc@colors == "antiquewhite3" | mc@colors=="darkkhaki" | mc@colors=="gray"
f_ps = mc@colors == "#DABE99"
f_epi = mc@colors == "#635547"

base_3a = max(log2(1e-5+mc@e_gc["Dnmt3a",]))
base_3b = max(log2(1e-5+mc@e_gc["Dnmt3b",]))
fmes = f_meso | f_ps | f_epi
fend = f_endo | f_ps | f_epi
png("paper_figs/fig1/eb_3a3b_rna.png",w=800,h=600)
layout(matrix(1:4,nrow=2))
par(mar=c(2,2,1,2))
plot(mc_t_e[fmes], legc["Dnmt3a",fmes]-base_3a, pch=21,bg=mc@colors[fmes],cex=1.2, ylim=c(-5,0), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(v=6, lty=2,lwd=2)
abline(h=seq(-1,-5,l=5),col="gray", lty=2)
plot(mc_t_e[fmes], legc["Dnmt3b",fmes]-base_3b, pch=21,bg=mc@colors[fmes],cex=1.2, ylim=c(-5,0), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(v=6, lty=2,lwd=2)
abline(h=seq(-1,-5,l=5),col="gray", lty=2)
plot(mc_t_e[fend], legc["Dnmt3a",fend]-base_3a, pch=21,bg=mc@colors[fend],cex=1.2, ylim=c(-5,0), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(v=6, lty=2,lwd=2)
abline(h=seq(-1,-5,l=5),col="gray", lty=2)
plot(mc_t_e[fend], legc["Dnmt3b",fend]-base_3b, pch=21,bg=mc@colors[fend],cex=1.2, ylim=c(-5,0), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(v=6, lty=2,lwd=2)
abline(h=seq(-1,-5,l=5),col="gray", lty=2)
dev.off()
