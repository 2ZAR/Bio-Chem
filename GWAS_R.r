#Practice#

#Quantitative trait locus, QTL - GWAS#

> install.packages("qtl2")
> library(qtl2)
> iron <- read_cross2(file=system.file("extdata", "iron.zip", package="qtl2"))
> summary(iron)

> names(iron)

> head(iron$gmap)

> ironmap <- insert_pseudomarkers(map=iron$gmap, step=1)
> head(ironmap, n=2)

> ironpr <- calc_genoprob(cross=iron, map=ironmap, error_prob=0.002)
> names(ironpr)

> dimnames(ironpr$'19')

> (ironpr$'19')[1:3,,"D19Mit68"]

> (ironpr$'19')[1:3,,"c19.loc4"]

> (ironpr$'19')[1:3,,"c19.loc5"]

> plot_genoprob(ironpr, ironmap, ind=1, chr=19)
> grav <- read_cross2(file=system.file('extdata'. 'grav2.zip', package='qtl2'))
> summary(grav)

> names(grav)

> head(grav$gmap)

> gravmap <- insert_psudomarkers(map=grav$gmap, step=1)
> head(gravmap, n=2)

> gravpr <- calc_genoprob(cross=grav, map=gravmap)
> (gravpr$'1')[1:5,,"PVV4"]

> (gravpr$'1')[1:5,,"c1.loc1"]

> (gravpr$'1')[1:5,,"c1.loc2"]

> xcovar <- get_x_covar(iron)
> head(xcovar)

> out <- scan1(genoprobs=ironpr, pheno=iron$pheno, xcovar=xcovar)
> head(out, n=10)

> plot_scan1(out, map=ironmap, lodcolumn="liver")
> plot_scan1(out, map=ironmap, lodcolumn="spleen")

> operm <- scan1perm(genoprobs=ironpr, pheno=iron$pheno, xcovar=xcovar, n_perm=1000)
> hist(operm[,'liver'], breaks=50, xlab="LOD", main="LOD scores for liver scan with threshold in red")
> abline(v=summary(operm)[,'liver'], col='red', lwd=2)
> summary(operm)

> summary(operm, alpha=c(0.2, 0.05))

> operm2 <- scan1perm(ironpr, iron$pheno, xcovar=xcovar, n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(ironmap))
> summary(operm2, alpha=c(0.2, 0.05))

> shuffled_order <- sample(rownames(iron$pheno))
> pheno_permuted <- iron$pheno
> rownames(pheno_permuted) <- shuffled_order
> xcovar_permuted <- xcovar
> rownames(xcovar_permuted) <- shuffled_order
> out_permuted <- scan1(genoprobs=ironpr, pheno=pheno_permuted, xcovar=xcovar_permuted)
> plot(out_permuted, ironmap)
> head(shuffled_order)
> operm <- scan1perm(ironpr, iron$pheno, xcovar=xcovar, n_perm=1000)
> thr <- summary(operm)
> find_peaks(scan1_output_out, map=ironmap, threshold=thr, prob=0.95, expand2markers=FALSE)
> find_peaks(scan1_output_out, map=ironmap, threshold=thr, peakdrop=1.8, prob=0.95, expand2markers=FALSE)

> find_peaks(out, ironmap, threshold=3, drop=2)

> find_peaks(out, ironmap, prob=0.90, expand2markers=FALSE)

> find_peaks(out, ironmap, prob=0.95, expand2markers=FALSE)

> n_samples <- 25
> heatmap(kinship[1:n_samples, 1:n_samples], symm=TRUE)
> grid <- calc_grid(map=iron$gmap, step=1)
> pr_grid <- probs_to_grid(probs=ironpr, grid=grid)
> kinship_grid <- calc_kinship(probs=pr_grid)

> out_pg <- scan(ironpr, iron$pheno, kinship=kinship, xcovar=xcovar)
> kinship_loco <- calc_kinship(ironpr, "loco")
> out_pg_loco <- scan1(ironpr, iron$pheno, kinship_loco, xcovar=xcovar)
> plot_scan1(out_pg_loco, map=ironmap, lodcolumn="liver", col="black")
> plot_scan1(out_pg, map=ironmap, lodcolumn="liver", col="blue", add=TRUE)
> plot_scan1(out, map=ironmap, lodcolumn="liver", col="green" add=TRUE)

> bin_pheno <- apply(iron$pheno, 2, function(a) as.numeric(a > median(a)))
> rownames(bin_pheno) <- rownames(iron$pheno)
> out_bin <- scan1(ironpr, bin_pheno, xcovar=xcovar, model="binary")
> par(mar=c(5.1, 4.1, 1.1, 1.1))
> ymx <- maxlod(out_bin)
> plot(out_bin, ironmap, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
> plot(out_bin, ironmap, lodcolumn=2, col="violetred", add=TRUE)
> legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out_bin), bg="gray90")
> find_peaks(out_bin, ironmap, threshold=3.5, drop=1.5)

> grav_kinship <- calc_kinship(gravpr, "loco")
> out_grav <- scan1(genoprobs=gravpr, pheno=grav$pheno, kinship=grav_kinship)
> which(out_grav == maxlod(out_grav), arr.ind=TRUE)

> plot(ou_grav, lodcolumn=133, map=gravmap)

> c2eff <- scan1coef(ironpr[,"2"], iron$pheno[,"liver"])
> dim(c2eff)

> head(c2eff)

> plot_coef(c2eff, ironmap, legend="topright")
> plot_coef(c2eff, ironmap, columns=1:3, scan1_output=out, main="Chromosome 2 QTL effects and LOD scores", legend="topright")
> c2effB <- scan1coef(ironpr[,"2"], iron$pheno[,"liver"], contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
> dim(c2effB)

> head(c2effB)

> plot_coef(c2effB, ironmap["2"], columns=2:3, col=c("green", "pink"))
> c2eff_pg <- scan1coef(ironpr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])
> dim(c2eff_pg)

> head(c2eff_pg)

> plot_coef(c2eff_pg, ironmap, columns=1:3, scan1_output=out_pg_loco, main="Chromosome 2 QTL effects and LOD scores", legend="topright")
> c2effB_pg <- scan1coef(ironpr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]], contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
> plot(c2effB_pg, ironmap["2"], columns=2:3, col=c("yellow", "purple"))
> c2blup <- scan1blup(ironpr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])
> plot_coef(c2eff, ironmap["2"], columns=1:3, add=TRUE, lty=2, legend="topright")
> c2eff_bin <- scan1coef(ironpr[,"2"], bin_pheno[,"liver"], model="binary")
> g <- maxmarg(ironpr, ironmap, chr=2, pos=28.6, return_char=TRUE)
> par(mar=c(4.1, 4.1, 0.6, 0.6))
> plot_pxg(g, iron$pheno[,"liver"], ylab="Liver phenotype")
> plot_pxg(g, iron$pheno[,"liver"], SEmult=2, swap,axes=TRUE, xlab="Liver phenotype")
> c16eff <- scan1coef(ironpr[,"16"], iron$pheno[,"liver"])
> plot_coef(c16eff, ironmap, legend="topright", scan1_output=out)

> DOex <- read_cross2(file = "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
> mpr <- calc_genoprob(DOex, error_prob=0.002)
> ampr <- genoprob_to_alleleprob(mpr)
> k <- calc_kinship(ampr, "loco"
> sex <- (DOex$covar$Sex == "male")*1
> names(sex) <- rownames(DOex$covar)
> sex <- setNames((DOex$covar$sex == "male")*1, rownames(DOex$covar))
> out <- scan1(ampr, DOex$pheno, k, sex)
> par(mar=c(4.1, 4.1, 0.6, 0.6))
> plot(out, DOex$gmap)
> coef_c2 <- scan1coef(ampr[,"2"], DOex$pheno, k[["2"]], sex)
> par(mar=c(4.1, 4.1, 0.6, 0.6))
> plot_coefCC(coef_c2, DOex$gmap["2"], bgcolor="gray95", legend="bottomleft")
> plot_coefCC(coef_c2, DOex$gmap["2"], scan1_output=out, bgcolor="gray95", legend="bottomleft")

> query_variants <- create_variant_query_func(경로1)
> variants_2_97.5 <- query_variants(2, 97, 98)
> query_genes <- create_gene_query_func(경로2)
> genes_2_97.5 <- query_genes(2, 97, 98)
> peak_Mbp <- max(out, DOex$pmap)$pos
> variants <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
> out_snps <- scan1snps(mpr, Doex$pmap, Doex$pheno, k[["2"]], sex, query_func=query_variants, chr=2, start=peak_Mbp - 1, end=peak_Mbp + 1, keep_all_snps=TRUE)
> par(mar=c(4.1, 4.1, 0.6, 0.6))
> plot_snpasso(out_snps$lod, out_snps$snpinfo)
> plot(out_snps$lod, out_snps$snpinfo)
> genes <- query_genes(2, peak_Mbp - 1, peak_Mbp + 1)
> plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)
> top <- top_snps(out_snps$lod, out_snps$snpinfo)
> print(top[,c(1, 8:15, 20)], row.names=FALSE)

> out_gwas <- scan1snps(mpr, DOex$pmap, DOex$pheno, k[[2]], sex, query_func=query_variants)
> par(mar=c(4.1, 4.1, 0.6, 0.6))
> plot(out_gwas$lod, out_gwas$snpinfo, altcol="green4", gap=0)
