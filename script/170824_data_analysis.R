source("handcraft_seurat-like_20161123.R")

###########################################
#                                          #
#          Load and prepare files          #
#                                          #
###########################################

load(
    file="../data/male_rpkm.Robj"
)

prot_coding_genes <- read.csv(file="../data/prot_coding.csv", row.names=1)
males <- male_rpkm[rownames(male_rpkm) %in% as.vector(prot_coding_genes$x),]

males <- males[,!colnames(males) %in% "E16.5_XY_20150202_C94_150331_8"]
males <- males[rowSums(males)>0,]


load(
    file="../data/male_count.Robj"
)

male_count <- male_count[rownames(male_count) %in% rownames(males),]
colnames(male_count) <- colnames(male_rpkm)
male_count <- male_count[,!colnames(male_count) %in% "E16.5_XY_20150202_C94_150331_8"]


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

tf_list <- read.csv(file="../data/mart_export_TF_GO0003700.txt", header=FALSE)

###########################################
#                                          #
#          Plot stats about cells          #
#                                          #
###########################################

pdf("../graph/20170908_graphs_qc_pcq.pdf", width=6)
male_pca <- prcomp(
    t(log(males+1)), 
    center=TRUE, 
    scale=TRUE
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=male_stages, 
    colours=male_stagePalette
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=male_clustering, 
    colours=male_clusterPalette
)

dev.off()


# rle(male_captures)
# rle(male_stages)
# rle(paste(male_stages, male_captures))

gene_per_cell_males <- males
gene_per_cell_males <- rowSums(gene_per_cell_males)
gene_per_cell_males[gene_per_cell_males>0] <- 1

sum(gene_per_cell_males)



gene_per_cell_males <- data.frame(
    cells=names(gene_per_cell_males),
    geneNb = gene_per_cell_males
)


gene_per_cell_males <- males
gene_per_cell_males[gene_per_cell_males>0] <- 1
gene_per_cell_males <- colSums(gene_per_cell_males)

gene_per_cell_males <- data.frame(
    cells=names(gene_per_cell_males),
    geneNb = gene_per_cell_males
)

median(gene_per_cell_males$geneNb)



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
    # ggtitle(paste("Median=", median(gene_per_cell_males$geneNb)," genes per cell", sep="")) +
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
    # geom_smooth(method=lm, color="red", se = FALSE) +
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



pdf("../graph/20170908_graphs_qc_pca_batch.pdf", width=6)


e11.5 <- males[,colnames(males) %in% names(male_stages[male_stages=="E11.5"])]
e11.5 <- e11.5[rowSums(e11.5)>0,]
e11.5_captures <- as.factor(sapply(strsplit(colnames(e11.5), "_"), `[`, 3))
levels(e11.5_captures) <- c("Capture 1", "Capture 2")
e11.5_palette <- c("#0000ff", "#8a00b0")


male_pca <- prcomp(
    t(log(e11.5+1)), 
    center=TRUE, 
    scale=TRUE
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=e11.5_captures, 
    colours=e11.5_palette
)


e12.5 <- males[,colnames(males) %in% names(male_stages[male_stages=="E12.5"])]
e12.5 <- e12.5[rowSums(e12.5)>0,]
e12.5_captures <- as.factor(sapply(strsplit(colnames(e12.5), "_"), `[`, 3))
levels(e12.5_captures) <- c("Capture 1", "Capture 2")
e12.5_palette <- c("#ff0000", "#ff7674")


male_pca <- prcomp(
    t(log(e12.5+1)), 
    center=TRUE, 
    scale=TRUE
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=e12.5_captures, 
    colours=e12.5_palette
)


e13.5 <- males[,colnames(males) %in% names(male_stages[male_stages=="E13.5"])]
e13.5 <- e13.5[rowSums(e13.5)>0,]
e13.5_captures <- as.factor(sapply(strsplit(colnames(e13.5), "_"), `[`, 3))
levels(e13.5_captures) <- c("Capture 1", "Capture 2")
e13.5_palette <- c("#803c00", "#f77f05")


male_pca <- prcomp(
    t(log(e13.5+1)), 
    center=TRUE, 
    scale=TRUE
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=e13.5_captures, 
    colours=e13.5_palette
)


e16.5 <- males[,colnames(males) %in% names(male_stages[male_stages=="E16.5"])]
e16.5 <- e16.5[rowSums(e16.5)>0,]
e16.5_captures <- as.factor(sapply(strsplit(colnames(e16.5), "_"), `[`, 3))
levels(e16.5_captures) <- c("Capture 1", "Capture 2")
e16.5_palette <- c("#7ca328", "#f9db21")


male_pca <- prcomp(
    t(log(e16.5+1)), 
    center=TRUE, 
    scale=TRUE
)

plot_pca(
    pca=male_pca, 
    pc=2, 
    conditions=e16.5_captures, 
    colours=e16.5_palette
)



dev.off()




###########################################
#                                          #
#            Var Gene Selection              #
#                                          #
###########################################

males_data <- getMostVarGenes(males, fitThr=2)
males_data <- log(males_data+1)

###########################################
#                                          #
#              RtSNE Analysis              #
#                                          #
###########################################

male_sub_pca <- PCA(
    t(males_data), 
    ncp = ncol(males_data), 
    graph=FALSE
)

significant_pcs <- permutationPA(
    male_sub_pca$ind$coord, 
    B = 100, 
    threshold = 0.05, 
    verbose = TRUE, 
    seed = NULL
)

male_t_sne <- plot_tSNE(
    pca=male_sub_pca,
    pc=significant_pcs$r,
    iter=5000,
    conditions=male_stages,
    colours=male_stagePalette
)

res.pca <- PCA(
    t(males_data), 
    ncp = significant_pcs$r, 
    graph=FALSE
)

res.hcpc <- HCPC(
    res.pca, 
    graph = FALSE
    )

plot(res.hcpc, choice ="tree", cex = 0.6)

male_clustering <- res.hcpc$data.clust$clust
male_clustering <- paste("C", male_clustering, sep="")
names(male_clustering) <- colnames(males)


male_clusterPalette <- c(
        "#457cff", 
        "#aeff01", 
        "#00c5ec", 
        "#009900",  
        "#a38cff",
        "#8dfaff"
    )

plot_tSNE(
    tsne=male_t_sne, 
    conditions=male_clustering, 
    colours= male_clusterPalette
)

names(male_clustering) <- names(males)
male_cluster_name <- as.factor(male_clustering)


###########################################
#                                          #
#                 DE Analysis              #
#                                          #
###########################################


DE_male <- prepare_for_DE (
    male_count, 
    male_clustering, 
    male_stages
)


male_DE_genes <- findDEgenes(
    DE_male, 
    qvalue=0.05
)

de_clusters <- get_up_reg_clusters(
    males, 
    male_cluster_name, 
    male_DE_genes
)

write.csv(
    de_clusters, 
    quote = FALSE, 
    file="170905_male_DE_genes_per_clusters.csv"
)

###########################################
#                                          #
#               GO term DE genes              #
#                                          #
###########################################

de_genes <- read.csv(file="170814_male_DE_genes_per_clusters.csv", row.names=1)

de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes

#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,c("genes", "cluster")]

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
    cutoff=0.3, 
    by="p.adjust", 
    select_fun=min
)


write.csv(formula_res@compareClusterResult, file="compared_GO_term_DE_cluster.csv")

dotplot(formula_res, showCategory=5)
dotplot(lineage1_ego, showCategory=5)



###########################################
#                                          #
#          Marker gene enrichment          #
#                                          #
###########################################

marker_gene_list <- read.csv(file="../data/marker_gene_list.txt", header=TRUE)
subset <- as.matrix(log(males[rownames(males) %in% marker_gene_list$gene,]+1))
mean_exp_marker_cT <- data.frame(matrix(ncol = 4))
colnames(mean_exp_marker_cT) <- c("median", "genes","cellType", "cluster")

for(cluster in unique(mean_exp_marker_cT$cluster)){
    data <- mean_exp_marker_cT[mean_exp_marker_cT$cluster==cluster,]
    test <- compare_means(median ~ cellType, data = data, ref.group = ".all.", method = "t.test")
    p <- ggboxplot(data, x = "cellType", y = "median", color = "white")+
      geom_violin(scale="width", aes(fill=cellType))+
      geom_boxplot(fill="white", outlier.shape = NA, width = 0.2)+
      geom_jitter(size=0.8, width = 0.2)+
      theme_bw()+
      ggtitle(cluster)+
    stat_compare_means(method = "anova", label.y = 10)+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", label.y = 8.5)+ # Pairwise comparison against all
    theme(
        text=element_text(size=45),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        legend.position= "none",
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        aspect.ratio=0.5
    )

    print(p)
}

###########################################
#                                          #
#         Diffusion map Analysis                 #
#                                          #
###########################################

male_dm <- run_diffMap(
    males_data, 
    male_clustering,
    sigma=19
)

plot_eigenVal(
    dm=male_dm
)


plot_dm_3D(
    dm=male_dm, 
    dc=c(1:3),
    condition=male_clustering, 
    colour=male_clusterPalette
)


plot_dm_3D(
    dm=male_dm, 
    dc=c(1:3),
    condition=male_stages, 
    colour=male_stagePalette
)


male_clusterPalette <- c(
        "#457cff", 
        "#aeff01", 
        "#00c5ec", 
        "#009900",  
        "#a38cff",
        "#8dfaff"
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
    dim=c(1:5), 
    condition=factor(male_clustering),
    start="C2",
    end=c("C3", "C5", "C6")
)


plot_dm_3D(
    dm=male_dm, 
    dc=c(1:3),
    condition=male_stages, 
    colour=male_stagePalette
)
plot3d(male_lineage, dim=c(1:3), add=TRUE, lwd=5)


plot_dm_3D(
    dm=male_dm, 
    dc=c(1:3),
    condition=male_clustering, 
    colour=male_clusterPalette
)
plot3d(male_lineage, dim=c(1:3), add=TRUE, lwd=5)

male_pseudotime <- get_pseudotime(male_lineage, wthres=0.90)
rownames(male_pseudotime) <- colnames(males)

male_clusterPalette2 <- c(
    "#518ecf", 
    "#c98b79", 
    "#457cff"
)

plot_smoothed_gene_per_lineage(
    rpkm_matrix=males, 
    pseudotime=male_pseudotime, 
    lin=c(1:3),
    gene="Wnt5a", 
    stages=male_stages, 
    clusters=male_clustering, 
    stage_colors=male_stagePalette,
    cluster_colors=male_clusterPalette,
    lineage_colors=male_clusterPalette2
)


plot_smoothed_genes <- function(genes, lin){
    male_clusterPalette2 <- c("#518ecf", "#c98b79", "#457cff")
    for (gene in genes){
        plot_smoothed_gene_per_lineage(
            rpkm_matrix=males, 
            pseudotime=male_pseudotime, 
            lin=lin,
            gene=gene, 
            stages=male_stages, 
            clusters=male_cluster_name, 
            stage_colors=male_stagePalette,
            cluster_colors=male_clusterPalette,
            lineage_colors=male_clusterPalette2
        )
    }
}


gene_list_lineage1_3 <- c(
    "Pdlim3",
    "Fblim1",
    "Pbx1",
    "Dmrt1",
    "Sfrp1",
    "Ctnnal1",
    "Ptgr1",
    "Dmrt1",
    "Pdgfa"
)

gene_list_lineage3 <- c(
    "Cdk1",
    "Cbx2",
    "Gata3",
    "Meis2",
    "Ren1",
    "Robo2",
    "Ar",
    "Wt1",
    "Arx",
    "Pdgfra",
    "Gli2",
    "Tcf21",
    "Ptch1"
)

gene_list_lineage1 <- c(
    "Nr0b1",
    "Sry",
    "Fgf9",
    "Sox9",
    "Cst8",
    "Dhh",
    "Hsd17b1"
)


plot_smoothed_genes(gene_list_lineage1_3, c(1,2))
plot_smoothed_genes(gene_list_lineage1, c(1))
plot_smoothed_genes(gene_list_lineage3, c(2))

###########################################
#                                          #
#            DE genes lineages              #
#                                          #
###########################################

de_clusters <- read.csv(
    file="170905_male_DE_genes_per_clusters.csv"
)

de_clusters <- de_clusters[de_clusters$qval<0.00001,]
de_clusters_l1_l2 <- de_clusters[de_clusters$Associated.cluster..max.mean. %in% c("C2", "C3", "C4", "C6"),]

de_matrix <- males[rownames(males) %in% de_clusters_l1_l2$genes,]

L1_lineage <- male_pseudotime[!is.na(male_pseudotime[,1]),1]
L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), L1_rpkm_exp_lineage))]


L2_lineage <- male_pseudotime[!is.na(male_pseudotime[,2]),2]
L2_ordered_lineage <- L2_lineage[order(L2_lineage, decreasing = TRUE)]
L2_rpkm_exp_lineage <- de_matrix[,names(L2_ordered_lineage)]
L2_cells <- L2_rpkm_exp_lineage[,order(match(names(L2_ordered_lineage), L2_rpkm_exp_lineage))]

L1_lineage_cells <- names(L1_lineage)
L2_lineage_cells <- names(L2_lineage)

comp_list <- comparelists(L1_lineage_cells, L2_lineage_cells)

common_cells <- comp_list$intersect
L1_spe_cells <- L1_lineage_cells[!L1_lineage_cells %in% comp_list$intersect]
L2_spe_cells <- L2_lineage_cells[!L2_lineage_cells %in% comp_list$intersect]

L1_cellLin <- c(
    rep_along("common cells", common_cells), 
    rep_along("L1 cells", L1_spe_cells)
)
names(L1_cellLin) <- c(common_cells, L1_spe_cells)

L1_cellLin <- L1_cellLin[order(match(names(L1_cellLin), colnames(L1_cells)))]

L2_cellLin <- c(
    rep_along("common cells", common_cells), 
    rep_along("L2 cells", L2_spe_cells)
)
names(L2_cellLin) <- c(common_cells, L2_spe_cells)
L2_cellLin <- L2_cellLin[order(match(names(L2_cellLin), colnames(L2_cells)))]


cellType_L1 <- male_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)

cellType_L2 <- male_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)


cellLin <- c(
    L2_cellLin,
    L1_cellLin
)

data_heatmap <- log(data.frame(
    L2_cells,
    L1_cells
)+1)


cellType <- c(
    cellType_L2,
    cellType_L1
)

annotation_col <- data.frame(
    cellLineages=cellLin,
    cellType=cellType,
    Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)

rownames(annotation_col) <- colnames(data_heatmap)

cellTypeCol <- c(
    C1="#a38cff" , # leyd
    C2="#aeff01", # preS
    C3="#00c5ec", # endo
    C4="#009900", # sert
    C5="#457cff", # prog
    C6="#8dfaff" # leyd
)

names(cellTypeCol) <- unique(cellType)

annotation_colors <- list(
    cellType=cellTypeCol,
        Stages=c(
            E10.5="#2754b5", 
            E11.5="#8a00b0", 
            E12.5="#d20e0f", 
            E13.5="#f77f05", 
            E16.5="#f9db21"
        )
)

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(11)), warm(14))
breaksList = seq(-2.2, 2.5, by = 0.2)

gene_clustering <- pheatmap(
    data_heatmap, 
    scale="row",
    kmeans_k=15,
    # breaks=c(-0.5, 0, 0.5),
    gaps_col=length(cellType_L2),
    show_colnames=FALSE, 
    # show_rownames=FALSE, 
    cluster_cols=FALSE,
    # cluster_rows=FALSE,
    # cutree_rows=6,
    clustering_method="ward.D",
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    # color=viridis(8)
    color=mypalette,
    breaks=breaksList
)

write.csv(gene_clustering$kmeans$cluster, file="male_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k15_scaled.csv")

###########################################
#                                          #
#            Loess smoothing                  #
#                                          #
###########################################

male_lineage1_sig_gene_pseudoT <- get_var_genes_pseudotime(
    males, 
    male_count, 
    male_pseudotime, 
    lineageNb=1, 
    male_cluster_name
)

write.csv(male_lineage1_sig_gene_pseudoT, file="male_lineage1_gene_var_genes_pseudotime_df_3_new_psdt.csv")

male_lineage1_clustering <- get_gene_clustering(
    male_lineage1_sig_gene_pseudoT, 
    males, 
    male_pseudotime, 
    lineageNb=1, 
    male_cluster_name,
    clusterNb=13,
    qvalue=0.05
)

write.csv(male_lineage1_clustering, file="male_lineage1_gene_clustering_qval_0.05_df_3_new_psdt.csv")

male_lineage1_clustering <- data.frame(genes=rownames(male_lineage1_clustering), clusters=male_lineage1_clustering)
tf_dynamics <- male_lineage1_clustering[rownames(male_lineage1_clustering) %in% tf_list$V1,]
write.csv(tf_dynamics, file="male_lineage1_gene_clustering_qval_0.05_TF_df_3_new_psdt.csv")

#############################################################
male_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
    males, 
    male_count, 
    male_pseudotime, 
    lineageNb=2, 
    male_cluster_name
)

male_lineage2_clustering <- get_gene_clustering(
    male_lineage2_sig_gene_pseudoT, 
    males, 
    male_pseudotime, 
    lineageNb=2, 
    male_cluster_name,
    clusterNb=9,
    qvalue=0.05
)

write.csv(male_lineage2_clustering, file="male_lineage2_gene_clustering_qval_0.05_df_3_new_psdt.csv")

male_lineage2_clustering <- data.frame(genes=rownames(male_lineage2_clustering), clusters=male_lineage2_clustering)
tf_dynamics <- male_lineage2_clustering[rownames(male_lineage2_clustering) %in% tf_list$V1,]

write.csv(tf_dynamics, file="male_lineage2_gene_clustering_qval_0.05_df_3_new_psdt_TF.csv")

#############################################################
male_lineage3_sig_gene_pseudoT <- get_var_genes_pseudotime(
    males, 
    male_count, 
    male_pseudotime, 
    lineageNb=3, 
    male_cluster_name
)
write.csv(male_lineage3_sig_gene_pseudoT, file="male_lineage3_gene_var_genes_pseudotime_df_3_new_psdt.csv")

male_lineage3_clustering <- get_gene_clustering(
    male_lineage3_sig_gene_pseudoT, 
    males, 
    male_pseudotime, 
    lineageNb=3, 
    male_cluster_name,
    clusterNb=10,
    qvalue=0.05
)

write.csv(male_lineage3_clustering, file="male_lineage3_gene_clustering_qval_0.05_df_3_new_psdt.csv")

male_lineage3_clustering <- data.frame(genes=rownames(male_lineage3_clustering), clusters=male_lineage3_clustering)
tf_dynamics <- male_lineage3_clustering[rownames(male_lineage3_clustering) %in% tf_list$V1,]
write.csv(tf_dynamics, file="male_lineage3_gene_clustering_qval_0.05_df_3_new_psdt_TF.csv")
