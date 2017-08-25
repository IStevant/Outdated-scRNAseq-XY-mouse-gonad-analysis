source("170824_functions.R")

###########################################
#                                         #
#          Load and prepare files         #
#                                         #
###########################################

load(
	file="../data/male_rpkm.Robj"
)

load(
	file="../data/male_count.Robj"
)

colnames(male_count) <- colnames(male_rpkm)
males <- male_rpkm

male_stages <- sapply(strsplit(colnames(males), "_"), `[`, 1)
names(male_stages) <- colnames(males)

male_captures <- sapply(strsplit(colnames(males), "_"), `[`, 3)
male_stagePalette <- c(
	"#2754b5", 
	"#8a00b0", 
	"#d20e0f", 
	"#f77f05", 
	"#f9db21"
)

males_data <- log(males+1)


###########################################
#                                         #
#          Plot stats about cells         #
#                                         #
###########################################


rle(male_captures)
rle(male_stages)
rle(paste(male_stages, male_captures))


gene_per_cell_males <- males
gene_per_cell_males[gene_per_cell_males>0] <- 1
gene_per_cell_males <- colSums(gene_per_cell_males)

gene_per_cell_males <- data.frame(
	cells=names(gene_per_cell_males),
	geneNb = gene_per_cell_males
)
sf1_vs_gfp <- as.data.frame(
	t(
		log(males[c("eGFP", "Nr5a1"),]+1)
	)
)

cor(sf1_vs_gfp, method="spearman")
cor(sf1_vs_gfp, method="pearson")

pdf("../graph/male_gene_per_cell_gfp_sf1.pdf", width=5)

ggplot(gene_per_cell_males, aes(geneNb)) +
	geom_histogram(color="black", fill="grey") +
	theme_bw() +
	labs(x = "Detected genes", y="Cell count")+
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)

ggplot(sf1_vs_gfp, aes(x=Nr5a1, y=eGFP)) +
	geom_point(color="black") +
	geom_smooth(method=lm, color="red", se = FALSE) +
	theme_bw() +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)


dev.off()



###########################################
#                                         #
#            Var Gene Selection           #
#                                         #
###########################################

males_data <- log(males+1)

male_var_genes <- getMostVarGenes(
	data=males,
	binNb=60,
	zScore=1.5
)

males_data <- males_data[rownames(males_data) %in% male_var_genes,]

###########################################
#                                         #
#              RtSNE Analysis             #
#                                         #
###########################################

male_sub_pca <- prcomp(
	t(males_data), 
	center=TRUE, 
	scale=TRUE
)

plot_percent_var(
	pca=male_sub_pca, 
	pc=30
)


plot_pca(
	pca=male_sub_pca, 
	pc=2, 
	conditions=male_stages, 
	colours=male_stagePalette
)

male_t_sne <- plot_tSNE(
	pca=male_sub_pca,
	pc=9,
	iter=5000,
	conditions=male_stages,
	colours=male_stagePalette
)

male_clustering <- plot_clusters(
	tsne=male_t_sne, 
	dist=3.5
	)

male_clustering <- paste("XY", male_clustering, sep="_")
names(male_clustering) <- names(males)

male_clusterPalette <- c(
	"#457cff", 
	"#8ae400", 
	"#00c5ec", 
	"#009900",
	"#8dfaff"
)

male_cluster_name <- as.factor(male_clustering)
levels(male_cluster_name) <- c("Progenitors", "Pre-Sertoli", "Endothelial", "Sertoli", "Leydig")



cluster_stage <- paste(male_stages, male_cluster_name, sep="_")
cluster_stage <- cluster_stage[order(cluster_stage)]

cluster_stage <-rle(cluster_stage)

cluster_stage <-data.frame(
	cluster_stage$lengths,
	cluster_stage$values
	)


###########################################
#                                         #
#               DE Analysis               #
#                                         #
###########################################


DE_male <- prepare_for_DE (
	male_count, 
	male_cluster_name, 
	male_stages
)


male_DE_genes <- findDEgenes(
	DE_male, 
	qvalue=0.05
)

de_clusters <- get_up_reg_clusters(
	male_count, 
	male_cluster_name, 
	male_DE_genes
)

write.csv(
	de_clusters, 
	quote = FALSE, 
	file="170814_male_DE_genes_per_clusters.csv"
)

test_male_DE_genes <- get_top_up_reg_clusters(
	male_count, 
	male_clustering, de_clusters, 
	100
)

top_de_clusters <- de_clusters[de_clusters$genes %in% test_male_DE_genes,]

matrix <- log(males[rownames(males) %in% test_male_DE_genes,]+1)

plot_heatmap(
	matrix=matrix, 
	clusters=male_cluster_name,
	stages=male_stages
)


###########################################
#                                         #
#        GO terms of the DE genes         #
#                                         #
###########################################

library("clusterProfiler")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")


# GO terms

de_genes <- read.csv(file="170814_male_DE_genes_per_clusters.csv", row.names=1)
gene_names <- rownames(subset(de_genes, qval<0.01))

#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

de_gene_clusters <- de_gene_names[which(de_gene_names$genes==entrez_genes$SYMBOL),c("genes", "cluster")]

de_gene_clusters <- data.frame(
	ENTREZID=entrez_genes$ENTREZID[which(entrez_genes$SYMBOL==de_gene_clusters$genes)],
	cluster=de_gene_clusters$cluster
)

list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, de_gene_clusters$cluster)

formula_res <- compareCluster(
	ENTREZID~cluster, 
	data=de_gene_clusters, 
	fun="enrichGO", 
	OrgDb="org.Mm.eg.db",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
)

lineage1_ego <- simplify(
	formula_res, 
	cutoff=0.4, 
	by="p.adjust", 
	select_fun=min
)

dotplot(lineage1_ego, showCategory=8)


pdf("../graph/GO_term_DE_genes_clusters.pdf", width=10, height=12)
dotplot(formula_res, showCategory=10)
dev.off()

###########################################
#                                         #
#          heatmap_market_genes           #
#                                         #
###########################################

markerGenes <- c(
	"Nr2f1",
	"Lhx9",
	"Tcf21",
	"Pdgfra",
	"Wnt5a",
	"Arx",
	"Runx1",
	"Sry",
	"Gadd45g",
	"Mro",
	"Aard",
	"Amh",
	"Dhh",
	"Hsd17b3",
	"Cyp11a1",
	"Cyp17a1",
	"Star",
	"Insl3",
	"Jack3",
	"Pecam1",
	"Cdh5",
	"Flt4",
	"Esam"
)

gene_subset <- as.matrix(log(males[rownames(males) %in% markerGenes,]+1))

cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(male_cluster_name[male_cluster_name=="Progenitors"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(male_cluster_name[male_cluster_name=="Pre-Sertoli"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(male_cluster_name[male_cluster_name=="Endothelial"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(male_cluster_name[male_cluster_name=="Sertoli"])]
cl5_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(male_cluster_name[male_cluster_name=="Leydig"])]

heatmap_gene_subset <- cbind(
	cl1_gene_subset, 
	cl2_gene_subset, 
	cl3_gene_subset,
	cl4_gene_subset,
	cl5_gene_subset
	)

heatmap_male_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)

heatmap_male_clusters <- c(
	rep(1,ncol(cl1_gene_subset)), 
	rep(2,ncol(cl2_gene_subset)),
	rep(3,ncol(cl3_gene_subset)),
	rep(4,ncol(cl4_gene_subset)),
	rep(5,ncol(cl5_gene_subset))
	)

plot_heatmap_2(heatmap_gene_subset, heatmap_male_clusters)


#####################################################
#                                                   #
#  Prepare start and ends for lineage construction  #
#                                                   #
#####################################################

male_clustering[
	which(
		sapply(strsplit(names(male_clustering), "_"), `[`, 1)=="E10.5" & male_clustering=="XY_1"
	)] <- "XY_1_E10.5"


male_clustering[
	which(
		sapply(strsplit(names(male_clustering), "_"), `[`, 1)=="E16.5" & male_clustering=="XY_1"
	)] <- "XY_1_E16.5"



male_clustering[
	which(
		male_clustering=="XY_3"
	)] <- "XY_1"

###########################################
#                                         #
#            Diffusion maps               #
#                                         #
###########################################

male_dm <- run_diffMap(
	males_data, 
	male_clustering,
	sigma="global"
)

plot_eigenVal(
	dm=male_dm
)


male_clusterPalette <- c(
	"#6a7cff", 
	"#8ae400", 
	"#00ccff", 
	"#008000",
	"#55ffdd"
)


plot_dm_2D(
	dm=male_dm, 
	condition=male_cluster_name, 
	colour=male_clusterPalette
)


plot_dm_2D(
	dm=male_dm, 
	condition=male_stages, 
	colour=male_stagePalette
)

plot_dm_3D(
	dm=male_dm, 
	dc=c(1:3),
	condition=male_cluster_name, 
	colour=male_clusterPalette
)

plot_dm_3D(
	dm=male_dm, 
	dc=c(1:3),
	condition=male_stages, 
	colour=male_stagePalette
)

male_lineage <- get_lineage(
	dm=male_dm, 
	dim=c(1:6), 
	condition=factor(male_clustering),
	start="XY_1_E10.5", 
	end=c("XY_1_E16.5", "XY_4", "XY_5")
)


plot_dm_3D(
	dm=male_dm, 
	dc=c(1:3),
	condition=male_stages, 
	colour=male_stagePalette
)
plot3d(male_lineage, dim=c(1:3), add=TRUE, lwd=5)
rgl.postscript( "../graph/male_dm_stages_new_cols.svg", fmt = "svg", drawText = TRUE )

plot_dm_3D(
	dm=male_dm, 
	dc=c(1:3),
	condition=male_cluster_name, 
	colour=male_clusterPalette
)
plot3d(male_lineage, dim=c(1:3), add=TRUE, lwd=5)
rgl.postscript( "../graph/male_dm_lineages_new_cols.svg", fmt = "svg", drawText = TRUE )



male_pseudotime <- get_pseudotime(male_lineage, wthres=0.9)
rownames(male_pseudotime) <- colnames(males)



plot_gene_per_lineage(
	rpkm_matrix=males, 
	pseudotime=male_pseudotime, 
	gene="Nr2f2", 
	stages=male_stages, 
	clusters=male_cluster_name, 
	stage_colors=male_stagePalette,
	cluster_colors=male_clusterPalette2
)


plot_gene_per_lineage(
	rpkm_matrix=males, 
	pseudotime=male_pseudotime, 
	gene="Sox9", 
	stages=male_stages, 
	clusters=male_cluster_name, 
	stage_colors=male_stagePalette,
	cluster_colors=male_clusterPalette2
)


plot_gene_per_lineage(
	rpkm_matrix=males, 
	pseudotime=male_pseudotime, 
	gene="Hsd3b1", 
	stages=male_stages, 
	clusters=male_cluster_name, 
	stage_colors=male_stagePalette,
	cluster_colors=male_clusterPalette2
)


###########################################
#                                         #
#            Loess smoothing              #
#                                         #
###########################################

male_lineage1_sig_gene_pseudoT <- get_var_genes_pseudotime(
	males, 
	male_count, 
	male_pseudotime, 
	lineageNb=1, 
	male_cluster_name
)

write.csv(male_lineage1_sig_gene_pseudoT, file="male_lineage1_gene_var_genes_pseudotime_df_4.csv")

male_lineage1_sig_gene_pseudoT <- read.csv(file="male_lineage1_gene_var_genes_pseudotime_df_4.csv", row.names=1)


tiff("../graph/male_lineage1_var_genes_qval_0.05_df_4_new_cols.tif", res = 300, height = 20, width = 17, units = 'cm')

male_lineage1_clustering <- get_gene_clustering(
	male_lineage1_sig_gene_pseudoT, 
	males, 
	male_pseudotime, 
	lineageNb=1, 
	male_cluster_name,
	clusterNb=13,
	qvalue=0.05
)

dev.off()

write.csv(male_lineage1_clustering, file="male_lineage1_gene_clustering_qval_0.05_df_4.csv")
male_lineage1_clustering <- data.frame(genes=rownames(male_lineage1_clustering), clusters=male_lineage1_clustering)
tf_dynamics <- male_lineage1_clustering[rownames(male_lineage1_clustering) %in% tf_list$V1,]
write.csv(tf_dynamics, file="male_lineage1_gene_clustering_qval_0.05_TF_df_4.csv")


#############################################################
male_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
	males, 
	male_count, 
	male_pseudotime, 
	lineageNb=2, 
	male_cluster_name
)

tiff("../graph/male_lineage2_var_genes_qval_0.05_df_4_new_cols.tif", res = 300, height = 20, width = 17, units = 'cm')
male_lineage2_clustering <- get_gene_clustering(
	male_lineage2_sig_gene_pseudoT, 
	males, 
	male_pseudotime, 
	lineageNb=2, 
	male_cluster_name,
	clusterNb=10,
	qvalue=0.05
)
dev.off()

write.csv(male_lineage2_clustering, file="male_lineage2_gene_clustering_qval_0.05_df_4.csv")

male_lineage2_clustering <- data.frame(genes=rownames(male_lineage2_clustering), clusters=male_lineage2_clustering)
tf_dynamics <- male_lineage2_clustering[rownames(male_lineage2_clustering) %in% tf_list$V1,]

write.csv(tf_dynamics, file="male_lineage2_gene_clustering_qval_0.05_df_4.csv")


#############################################################
male_lineage3_sig_gene_pseudoT <- get_var_genes_pseudotime(
	males, 
	male_count, 
	male_pseudotime, 
	lineageNb=3, 
	male_cluster_name
)

tiff("../graph/male_lineage3_var_genes_qval_0.05_df_4_new_cols.tif", res = 300, height = 20, width = 17, units = 'cm')
male_lineage3_clustering <- get_gene_clustering(
	male_lineage3_sig_gene_pseudoT, 
	males, 
	male_pseudotime, 
	lineageNb=3, 
	male_cluster_name,
	clusterNb=7,
	qvalue=0.05
)
dev.off()

write.csv(male_lineage3_clustering, file="male_lineage3_gene_clustering_qval_0.05_df_4.csv")

male_lineage3_clustering <- data.frame(genes=rownames(male_lineage3_clustering), clusters=male_lineage3_clustering)

tf_dynamics <- male_lineage3_clustering[rownames(male_lineage3_clustering) %in% tf_list$V1,]

write.csv(tf_dynamics, file="male_lineage3_gene_clustering_qval_0.05_df_4.csv")



library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")


# Progenitors lineage

male_lineage3_clustering <- read.csv(file="male_lineage3_gene_clustering_qval_0.05_df_4.csv", row.names=1)
male_lineage3_clustering <- data.frame(clusters=male_lineage3_clustering, genes=rownames(male_lineage3_clustering))

lineage3_genes <- rownames(male_lineage3_clustering)

#convert gene ID into entrez genes
lineage3_entrez_genes <- bitr(lineage3_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
male_lineage3_clustering <- male_lineage3_clustering[male_lineage3_clustering$genes %in% lineage3_entrez_genes$SYMBOL,]
male_lineage3_clustering <- data.frame(entrez_genes=lineage3_entrez_genes$ENTREZID, clusters=male_lineage3_clustering)

male_lineage3_clustering[male_lineage3_clustering==3] <- "early"
male_lineage3_clustering[male_lineage3_clustering==4] <- "bimodal"

male_lineage3_clustering[male_lineage3_clustering==5] <- "intermediate"
male_lineage3_clustering[male_lineage3_clustering==6] <- "intermediate"
male_lineage3_clustering[male_lineage3_clustering==7] <- "intermediate"

male_lineage3_clustering[male_lineage3_clustering==1] <- "late"
male_lineage3_clustering[male_lineage3_clustering==2] <- "late"



lineage3_gene_clusters <- split(male_lineage3_clustering$entrez_genes, male_lineage3_clustering$clusters.Gene_Clusters)

lineage3_formula_res <- compareCluster(
	entrez_genes~clusters.Gene_Clusters, 
	data=male_lineage3_clustering, fun="enrichGO", 
	OrgDb="org.Mm.eg.db",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
)

lineage3_ego <- simplify(
	lineage3_formula_res, 
	cutoff=0.5, 
	# by="p.adjust", 
	select_fun=min
)

dotplot(lineage3_ego, showCategory=6)

dotplot(lineage3_formula_res, showCategory=10)

pdf("../graph/go_term_progenitor_pseudotime.pdf", width=10)
dotplot(lineage3_ego, showCategory=6)
dev.off()


###########################################
#                                         #
#     DE genes lineage 1 vs lineage 3     #
#                                         #
###########################################


male_lineage1_clustering <- read.csv(file="male_lineage1_gene_clustering_qval_0.05_df_4.csv", row.names=1)
male_lineage3_clustering <- read.csv(file="male_lineage3_gene_clustering_qval_0.05_df_4.csv", row.names=1)

de_genes_pseudotime <- unique(c(
	rownames(male_lineage1_clustering),
	rownames(male_lineage3_clustering)
	))


male_pseudotime_lineage1 <- male_pseudotime[,"curve1"]
male_pseudotime_lineage2 <- male_pseudotime[,"curve2"]
male_pseudotime_lineage3 <- male_pseudotime[,"curve3"]

male_pseudotime_lineage1[!is.na(male_pseudotime_lineage1)] <- 1
male_pseudotime_lineage2[!is.na(male_pseudotime_lineage2)] <- 2
male_pseudotime_lineage3[!is.na(male_pseudotime_lineage3)] <- 3

lineages <- paste(
	male_pseudotime_lineage1, 
	male_pseudotime_lineage2, 
	male_pseudotime_lineage3, 
	sep="_"
)

lineages <- data.frame(
	rownames(male_pseudotime),
	lineages,
	stringsAsFactors=FALSE
)


lineages[lineages=="1_NA_NA"] <- 1
lineages[lineages=="NA_2_NA"] <- 2
lineages[lineages=="NA_NA_3"] <- 3
lineages[lineages=="1_2_3"] <- "1_3"
lineages[lineages=="NA_2_3"] <- 3
lineages[lineages=="1_NA_3"] <- "1_3"
lineages[lineages=="1_2_NA"] <- 1

lineages_all <- paste("l", lineages$lineages, sep="")
names(lineages_all) <- lineages[,1]


lineages <- lineages_all

genes_pseudoT <- prepare_for_DE (
	male_count, 
	clustering=lineages_all, 
	stages=male_stages
)

lineages <- lineages[- grep("lNA_NA_NA", lineages)]

lineages_1_3 <- lineages[- grep("l2", lineages)]
lineages_1_3 <- lineages_1_3 [- grep("l1_3", lineages_1_3)]

genes_pseudoT <- genes_pseudoT[, names(lineages_1_3)]

genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)

genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 30, ]


DE_genes_pseudoT <- differentialGeneTest(
	genes_pseudoT, 
	fullModelFormulaStr="~cellType",
	cores = 3
)

DE_genes_pseudoT <- subset(DE_genes_pseudoT, qval< 0.05)


count <- male_count[,names(lineages_1_3)]

DE_genes_pseudoT <- get_up_reg_clusters(
	count, 
	lineages_all, 
	DE_genes_pseudoT
)



write.csv(DE_genes_pseudoT, file="male_lineage3_lineage1_DE_gene_pseudotime_qval_0.05_2.csv")


DE_genes_pseudoT <- read.csv(file="male_lineage3_lineage1_DE_gene_pseudotime_qval_0.05_2.csv", row.names=1)
DE_genes_pseudoT <- subset(DE_genes_pseudoT, qval<0.0001)

de_genes <- rownames(subset(DE_genes_pseudoT, cluster %in% c("l1","l3")))

lineage1_male <- males[de_genes ,names(lineages[lineages %in% c("l1_3", "l1")])]

lineage3_male <- males[de_genes ,names(lineages[lineages %in% c("l1_3", "l3")])]

lineage1_3_male <- males[de_genes ,names(lineages[lineages %in% c("l1_3", "l1", "l3")])]


# Get cells from each lineages
common_cells <- names(lineages[lineages %in% c("l1_3")])
l1_cells <- names(lineages[lineages %in% c("l1")])
l3_cells <- names(lineages[lineages %in% c("l3")])


# order cells with pseudo-time
common_cells_pseudotime <- male_pseudotime[common_cells,"curve1"]
common_cells_pseudotime <- names(common_cells_pseudotime[order(common_cells_pseudotime)])
common_cells <- common_cells[order(match(common_cells, common_cells_pseudotime))]

# order cells with pseudo-time
l1_cells_pseudotime <- male_pseudotime[l1_cells,"curve1"]
l1_cells_pseudotime <- names(l1_cells_pseudotime[order(l1_cells_pseudotime)])
l1_cells <- l1_cells[order(match(l1_cells, l1_cells_pseudotime))]

# order cells with pseudo-time
l3_cells_pseudotime <- male_pseudotime[l3_cells,"curve3"]
l3_cells_pseudotime <- names(l3_cells_pseudotime[order(l3_cells_pseudotime)])
l3_cells <- l3_cells[order(match(l3_cells, l3_cells_pseudotime))]


common_cells_matrix <- males[de_genes, common_cells]
common_cells_matrix <- common_cells_matrix[,order(match(colnames(common_cells_matrix), common_cells))]

common_cells_matrix <- males[de_genes, common_cells]
common_cells_matrix <- common_cells_matrix[,order(match(colnames(common_cells_matrix), common_cells))]


l1_common_cells_matrix <- common_cells_matrix
colnames(l1_common_cells_matrix) <- paste(colnames(l1_common_cells_matrix), "l1", sep="_")

l3_common_cells_matrix <- common_cells_matrix[,ncol(common_cells_matrix):1]

l1_cells_matrix <- males[de_genes, l1_cells]
l1_cells_matrix <- l1_cells_matrix[,order(match(colnames(l1_cells_matrix), l1_cells))]

l3_cells_matrix <- males[de_genes, l3_cells]
l3_cells_matrix <- l3_cells_matrix[,order(match(colnames(l3_cells_matrix), l3_cells))]
l3_cells_matrix <- l3_cells_matrix[,ncol(l3_cells_matrix):1]


data_heatmap <- log(data.frame(
	l3_cells_matrix,
	l3_common_cells_matrix,
	l1_common_cells_matrix,
	l1_cells_matrix
	)+1)


cellType <- c(
	as.vector(male_cluster_name[colnames(l3_cells_matrix)]),
	as.vector(rev(male_cluster_name[colnames(l3_common_cells_matrix)])),
	as.vector(male_cluster_name[colnames(l3_common_cells_matrix)]),
	as.vector(male_cluster_name[colnames(l1_cells_matrix)])
)



annotation_col <- data.frame(
	cellType=cellType,
	Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)

rownames(annotation_col) <- colnames(data_heatmap)

cellTypeCol = c(
	"#bd66d4", 
	"#8dcf38", 
	"#ffd65c",
	"#77b0f3"
)

names(cellTypeCol) <- unique(cellType)

annotation_colors <- list(
	cellType=cellTypeCol,
		Stages=c(
		E16.5="#fc8d59", 
		E13.5="#fee090", 
		E12.5="#e0f3f8", 
		E11.5="#91bfdb", 
		E10.5="#4575b4"
	)
)

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(15)), warm(16))
breaksList = seq(-3, 3, by = 0.2)



gene_clustering <- pheatmap(
	data_heatmap, 
 	scale="row",
 	kmeans_k=15,
 	gaps_col=length(c(colnames(l3_cells_matrix), colnames(l3_common_cells_matrix))),
	show_colnames=FALSE, 
	cluster_cols=FALSE,
	clustering_method="ward.D2",
	annotation_col=annotation_col,
	annotation_colors=annotation_colors,
	color=mypalette,
	breaks=breaksList
)


write.csv(gene_clustering$kmeans$cluster, file="male_lineage3_lineage1_DE_gene_pseudotime_qval_0.0001_gene_clustering_kmeans_k15_scaled.csv")
