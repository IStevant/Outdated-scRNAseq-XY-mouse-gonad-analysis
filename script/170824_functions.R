###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("taRifx")
library("matrixStats")
library("ggplot2")
library("Rtsne")
library("fpc")
library("factoextra")
library("monocle")
library("viridis")
library("gplots")
library("RColorBrewer")
library("destiny")
library("slingshot")
library("rgl")
library("scatterplot3d")
library("made4")
library("pheatmap")

r3dDefaults$windowRect <- c(0,50, 700, 700)

###########################################
#                                         #
#            Var Gene Selection           #
#                                         #
###########################################

evenbins <- function(x, bin.count=10, order=T) {
    bin.size <- rep(length(x) %/% bin.count, bin.count)
    bin.size <- bin.size + ifelse(1:bin.count <= length(x) %% bin.count, 1, 0)
    bin <- rep(1:bin.count, bin.size)
    return(factor(bin, levels=1:bin.count, ordered=order))
}


getMostVarGenes <- function(data=data,binNb=binNb,zScore=zScore){
	data_no0 <- as.matrix(
		data[log(rowSums(data)+1)>1,]
		)
	medianGeneExp <- rowMedians(data_no0)
	names(medianGeneExp)<- rownames(data_no0)
	orderedMedianGeneExp <- medianGeneExp[order(medianGeneExp)]
	geneBin <- evenbins(orderedMedianGeneExp, binNb)
	varGenes <- vector('character')

	for (bin in 1:binNb){
		binGeneName <- names(orderedMedianGeneExp[geneBin==bin])
		geneBinData <- data_no0[binGeneName,]
		meanBinExp <- rowMeans(geneBinData)
		varBinExp <- rowVars(geneBinData)
		dispBinVar <- log((varBinExp/meanBinExp)+1)
		z.scores <- lapply(dispBinVar,function(x)(x-mean(dispBinVar))/sd(dispBinVar))
		varGenes <- c(varGenes, names(z.scores[z.scores>zScore]))
	}

	varGenes <- unique(varGenes)
	print(
		paste(
			length(varGenes),
			"selected genes",
			sep=" "
		)
	)

	meanGeneExp <- log(rowMeans(data_no0)+1)
	sqrtCoefVar <- sqrt(
		rowSds(data_no0)/rowMeans(data_no0)
	)

	gene_color <- meanGeneExp
	gene_color[gene_color>=0] <- "#595959"
	gene_color[varGenes] <- "#5a9ca9"

	plotData <-  data.frame(
		x=meanGeneExp, 
		y=sqrtCoefVar,
		gene_col=gene_color 
		)
	p <- ggplot(plotData) + 
	geom_point(aes(x,y, color=(gene_col))) + 
	theme_bw() +
	labs(x = "log(mean expression)", y="squared coefficient of variation")+
	ggtitle(paste(length(varGenes), " selected genes (bin=", binNb,", zScore=", zScore,")", sep="")) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		legend.position= "none",
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=0.5
	)+
	scale_color_manual(
		values=c("#595959","#5a9ca9")
	)
	print(p)
	return(varGenes)
}



###########################################
#                                         #
#              RtSNE Analysis             #
#                                         #
###########################################

run_tSNE <- function(pca=pca, pc=pc, iter=2000, perplexity=0){
	set.seed(1)
	if (perplexity>0){
			rtsne_out <- Rtsne(pca$x[,pc] , pca=FALSE, max_iter=iter, verbose=TRUE, perplexity=perplexity)

		} else {
			rtsne_out <- Rtsne(pca$x[,pc] , pca=FALSE, max_iter=iter, verbose=TRUE)

		}
	tSNE <- data.frame(
		rtsne_out$Y
	)
	return(tSNE)
}

plot_tSNE <- function(pca=pca, pc=pc, iter=2000, conditions=conditions, colours=colours, perplexity=0){
	tsne <- run_tSNE(
		pca, 
		1:pc, 
		iter,  
		perplexity
	)
	tSNE <- data.frame(
		tsne,
		conditions
	)
	colnames(tSNE)<- c("tSNE_1", "tSNE_2", "cond")
	g<-ggplot(tSNE, aes(tSNE_1,tSNE_2)) +
	geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
	theme_bw() +
	scale_fill_manual(
		values=colours,
		name=""
		) +
		scale_color_manual(
			values=colours,
			name=""
		) +
	ggtitle(paste("t-SNE plot (", pc," PCs)", sep="")) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title = element_text(size =16 ,face="bold"),
		plot.title = element_text(size=18, face="bold", hjust = 0.5),
		aspect.ratio=1
	)
	print(g)
	return(tsne)
}


# t-SNE colored by gene expression level (Sup. Fig. 6)

tsne_gene_exp <- function(tsne_result, gene, rpkm){
	colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
	mypalette <- c(warm(20))


	p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
	geom_point(shape = 21, stroke=0.25, aes(fill=as.numeric(log(rpkm[gene,]+1))), color="white", size = 2) +
	# scale_fill_gradient2(high="darkred", low="yellow")+
	scale_fill_gradientn(colours = mypalette)+
	theme_bw() +
	ggtitle(gene) +
	theme(
		plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
		axis.text=element_text(size=12),
		axis.title=element_text(size=16),
		legend.text = element_text(size =12),
		legend.title=element_blank(),
		aspect.ratio=1,
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
	)
	print(p)

}




###########################################
#                                         #
#                 Plot PCA                #
#                                         #
###########################################

plot_pca <- function(pca=pca, pc=pc, conditions=conditions, colours=colours){
		# PCs <- 1:pc
		PCs <- paste("PC",1:pc, sep="")
		percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
		cond <- factor(conditions)
		col <- factor(conditions)
		levels(col) <- colours
		col <- as.vector(col)

		scores <- as.data.frame(pca$x)
		PCs.combinations <- combn(PCs,2)
		g <- apply(
			PCs.combinations,
			2,
			function(combination)
			{
				p1 <- ggplot(scores, aes_string(x=combination[1], y=combination[2])) +
				geom_point(shape = 21, size = 3, stroke=0.5, aes(fill=cond, colour=cond), alpha=6/10) +
				theme_bw() +
				scale_fill_manual(
					values=colours,
					name=""
				) +
				scale_color_manual(
					values=colours,
					name=""
				) +
				xlab(paste(combination[1], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit=3),"%)", sep=""))+
				ylab(paste(combination[2], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit=3),"%)", sep=""))+
				# ggtitle("PCA")+
				theme(
					axis.text=element_text(size=16),
					axis.title=element_text(size=16),
					legend.text = element_text(size =16),
					legend.title = element_text(size =16 ,face="bold"),
					plot.title = element_text(size=18, face="bold", hjust = 0.5),
					aspect.ratio=1#,
					# legend.position="none"
				)
				print(p1)
			}
		)
}


plot_pca_3D <- function(pca=pca, pc=c(1:3), conditions=conditions, colours=colours){
	cond <- factor(conditions)
	col <- factor(conditions)
	levels(col) <- colours
	col <- as.vector(col)

	PCs <- paste("PC",pc, sep="")


	data <- data.frame(
		pca$x[,PCs[1]], 
		pca$x[,PCs[2]], 
		pca$x[,PCs[3]]
	)

	colnames(data) <- PCs

	plot3d(
		data,
		col=col,
		size=6.5,
	)

	plot3d(
		data,
		size=8,
		add = TRUE
	)

}



###########################################
#                                         #
#             Plot screeplot              #
#                                         #
###########################################

plot_percent_var <- function(pca=pca, pc=pc){
	percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
	percent_var_explained <- data.frame(
		PC=1:length(percent_var_explained),
		percent_Var=percent_var_explained
	)

	sub_percent_var_explained <- percent_var_explained[1:pc,]

	p <- ggplot(sub_percent_var_explained, aes(x=PC, y=percent_Var)) + 
		geom_col()+
		theme_bw() +
		xlab("PCs") +
		ylab("% Variance") +
		ggtitle("Screeplot")+
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=0.5
		)

	print(p)
}


###########################################
#                                         #
#        Clustering with DBSCAN           #
#                                         #
###########################################

plot_clusters <- function(tsne=tsne, dist=3.5){
	set.seed(123)
	db <-dbscan(tsne, dist, MinPts = 2)
	g <- fviz_cluster(
		db, 
		tsne, 
		stand = FALSE, 
		geom = "point", 
		ellipse=FALSE, 
		ggtheme=theme_bw(), 
		xlab="tSNE_1", 
		ylab="tSNE_2", 
		pointsize=2.5, 
		shape=19, 
		show.clust.cent=FALSE) +
		ggtitle(paste("DBSCAN clustering (dist=", dist, ")", sep="")) +
		theme(
			axis.text=element_text(size=16),
			axis.title=element_text(size=16),
			legend.text = element_text(size =16),
			legend.title = element_text(size =16 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			aspect.ratio=1
		)

	print(g)
	return(db$cluster)
}

###########################################
#                                         #
#                   DE                    #
#                                         #
###########################################


prepare_for_DE <- function(count_matrix=count_matrix, clustering=clustering, stages=stages){
	# Prepare tables for monocle object
	expr_matrix <- as.matrix(count_matrix)
	sample_sheet <- data.frame(cells=names(count_matrix), stages=stages, cellType=clustering)
	rownames(sample_sheet)<- names(count_matrix)
	gene_annotation <- as.data.frame(rownames(count_matrix))
	rownames(gene_annotation)<- rownames(count_matrix)
	colnames(gene_annotation)<- "genes"
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = gene_annotation)

	# Create a CellDataSet from the relative expression levels
	HSMM <- newCellDataSet(
		as(expr_matrix, "sparseMatrix"),
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit=0.5,
		expressionFamily=negbinomial.size()
	)
	
	HSMM <- detectGenes(HSMM, min_expr = 5)
	# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
	HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]

	HSMM <- estimateSizeFactors(HSMM)
	HSMM <- estimateDispersions(HSMM)

	return(HSMM)
}


findDEgenes <- function(HSMM=HSMM, qvalue=qvalue){
	diff_test_res <- differentialGeneTest(
		HSMM,
		fullModelFormulaStr="~cellType",
		cores = 3
	)

	sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
	sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)

	print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
	print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

	diff_test_res <- subset(diff_test_res, qval< qvalue)

	return(diff_test_res)
}


plot_DE_heatmap <- function(rpkm_matrix=rpkm_matrix, diff_test_res=diff_test_res, qvalue=qvalue, condition=condition) {
	sig_genes <- subset(diff_test_res, qval < qvalue)
	dim(sig_genes)
	de_genes <- as.matrix(log(rpkm_matrix[rownames(rpkm_matrix) %in% rownames(sig_genes),]+1))
	plot_heatmap(de_genes, condition)
}



get_up_reg_clusters <- function(count, clustering, DE_genes){
	cluster_nb <- unique(clustering)
	mean_per_cluster <- vector()

	DE_genes <- DE_genes[order(rownames(DE_genes)),]
	count <- count[order(rownames(count)),]
	count_de_genes <- count[rownames(count) %in% DE_genes$genes,]

	for (clusters in cluster_nb) {
		mean <- rowMeans(
				count_de_genes[,
				colnames(count_de_genes) %in% names(clustering[clustering==clusters])
			]
		)
		mean_per_cluster <- cbind(
			mean_per_cluster,
			mean
		)
	}
	colnames(mean_per_cluster) <- cluster_nb
	up_reg_cluster <- colnames(mean_per_cluster)[apply(mean_per_cluster,1,which.max)]
	de_genes_table <- data.frame(
		DE_genes,
		cluster=up_reg_cluster
	)

	return(de_genes_table)
}

get_top_up_reg_clusters <- function(count, clustering, DE_genes, gene_nb){
	cluster_nb <- unique(DE_genes$cluster)
	DE_genes <- DE_genes[order(DE_genes$qval),]
	DE_genes <- DE_genes[DE_genes$qval<1.10e-5,]

	top_DE_genes <- vector()
	for (clusters in cluster_nb) {
		genes <- rownames(DE_genes[DE_genes$cluster==clusters,])[1:gene_nb]
		top_DE_genes <- c(top_DE_genes, genes)
	}

	return(top_DE_genes)
}

###########################################
#                                         #
#              Draw heatmaps              #
#                                         #
###########################################


plot_heatmap <- function(matrix=matrix, clusters=clusters, stages=stages){
	annotation_col <- data.frame(
		Cell_Clusters=clusters,
		Stages=stages
	)

	cell_cluster_colors <- c("#bd66d4", "#77b0f3", "#8dcf38", "#ffd65c", "#fb7072")
	names(cell_cluster_colors) <- unique(clusters)

	annotation_colors <- list(
		Stages=c(E16.5="#fc8d59", E13.5="#fee090", E12.5="#e0f3f8", E11.5="#91bfdb", E10.5="#4575b4"),
		Cell_Clusters=cell_cluster_colors
	)

	pheatmap(
		matrix,  
		show_colnames=FALSE, 
		show_rownames=FALSE, 
		clustering_method="ward.D2",
		annotation_col=annotation_col,
		annotation_colors=annotation_colors,
		color=viridis(100)
	)

}


plot_heatmap_2 <- function(matrix=matrix, conditions=conditions){
	cond <- factor(conditions)
	colours <- all_clusterPalette <- rainbow(length(unique(conditions)))

	colours <- colours[1:length(unique(cond))]
	col <- factor(cond)
	levels(col) <- colours
	col <- as.vector(col)
	hmcol = viridis(100)
	heatmap.2(
		matrix, 
		scale='none', 
		trace="none",
		Colv=FALSE,
		keysize = 1, 
		col = hmcol, 
		ColSideColors = col, 
	)

	par(xpd=TRUE, mar=c (5.1, 4, 4, 3))
	legend("topleft", inset=c(-0.1,0.1),legend=levels(cond),fill=colours)
}


###########################################
#                                         #
#             Diffusion map               #
#                                         #
###########################################


run_diffMap <- function(data=data, condition=condition, sigma="local"){
	destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
	destinyObj$condition <- factor(condition)
	dm <- DiffusionMap(destinyObj, sigma)
	return(dm)
}

plot_eigenVal <- function(dm=dm){
	linepad <- .5
	plot(
		eigenvalues(dm), 
		ylim = 0:1, 
		pch = 20, 
		xlab ='Diffusion component (DC)', 
		ylab ='Eigenvalue'
	)
}

plot_dm_2D <- function(dm=dm, dc=2, condition=condition, colours=colours){
	DCs <- paste("DC",1:dc, sep="")

	dm_eigen <- data.frame(
			dm@eigenvectors
		)

	DCs.combinations <- combn(DCs,2)
	g <- apply(
		DCs.combinations,
		2,
		function(combination)
		{
			p1 <- ggplot(dm_eigen, aes_string(x=combination[1], y=combination[2])) +
			geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=condition)) +
			theme_bw() +
			scale_fill_manual(
				values=colours,
				name=""
			) +
			xlab(combination[1])+
			ylab(combination[2])+
			ggtitle("Diffusion Map")+
			theme(
				axis.text=element_text(size=16),
				axis.title=element_text(size=16),
				legend.text = element_text(size =16),
				legend.title = element_text(size =16 ,face="bold"),
				plot.title = element_text(size=18, face="bold", hjust = 0.5),
				aspect.ratio=1
			)
			print(p1)
		}
	)

}


plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
	cond <- factor(condition)
	col <- factor(condition)
	levels(col) <- colours
	col <- as.vector(col)
	DCs <- paste("DC",dc, sep="")

	data <- data.frame(
		dm@eigenvectors[,DCs[1]], 
		dm@eigenvectors[,DCs[2]], 
		dm@eigenvectors[,DCs[3]]
	)
	colnames(data) <- DCs

	plot3d(
		data,
		col=col,
		size=6.5,
		box = FALSE
	)

	legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")

	plot3d(
		data,
		size=7.5,
		add = TRUE
	)

}




plot_dm_3D_transparent <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
	cond <- factor(condition)
	col <- factor(condition)
	levels(col) <- colours
	col <- as.vector(col)
	DCs <- paste("DC",dc, sep="")

	data <- data.frame(
		dm@eigenvectors[,DCs[1]], 
		dm@eigenvectors[,DCs[2]], 
		dm@eigenvectors[,DCs[3]]
	)
	colnames(data) <- DCs

	plot3d(
		data,
		col=col,
		size=6.5,
		alpha=0,
		box = FALSE
	)

	legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")

}




###########################################
#                                         #
#          Lineage prediction             #
#                                         #
###########################################

get_lineage <- function(dm=dm, dim=c(0:20), condition=condition, start=start, end=end){

	data <- data.frame(
		dm@eigenvectors[,dim]
	)
	crv <- slingshot(
		data, 
		condition, 
		start.clus = start, 
		end.clus=end
	)

	return(crv)
}


# Thanks to Guillaume Devailly for this function
rankKeepNA <- function(x) {
    return(
        ifelse(
            is.na(x),
            NA,
            rank(
            	x, 
            	na.last = TRUE, 
            	ties.method="random"
            )
        )
    )
}


get_pseudotime <- function(pseudotime, wthres=wthres){
	pseudoT <- list()
	for(lineage in 1:length(pseudotime@curves))local({
		curve <- pseudotime@curves[[lineage]]
		lambda <- curve$lambda
		weight <- curve$w
		ps <- curve$lambda
		ps[weight < wthres] <- NA
		ps <- rankKeepNA(ps)
		pseudoT[[lineage]] <<- ps
	})
	df <- t(do.call("rbind",pseudoT))
	colnames(df) <- names(pseudotime@curves)
	return(df)
}



plot_gene_per_lineage <- function(
	rpkm_matrix=rpkm_matrix, 
	pseudotime=pseudotime, 
	geneName=geneName, 
	stages=stages, 
	clusters=clusters, 
	stage_colors=stage_colors,
	cluster_colors=cluster_colors
	){

	myplots <- list()
	total_pseudotime <- vector()
	for (lineages in 1:ncol(pseudotime)){
		lineage <- as.vector(pseudotime[,lineages])
		total_pseudotime <- c(total_pseudotime, lineage)
		total_pseudotime <- na.omit(total_pseudotime)
	}
	max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
	max_pseudotime <- max(total_pseudotime)

	pseudotime_data <- data.frame(
		pseudotime=numeric(),
		lineage=numeric(),
		stages=character(),
		clusters=character(),
		gene=numeric()
	)

	for (lineages in 1:ncol(pseudotime)){

		i <- lineages

		lineage <- pseudotime[,lineages]
		sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)

		data <- data.frame(
			pseudotime=lineage,
			lineage=paste("Lineage ",lineages, sep=""),
			stages=stages,
			clusters=clusters,
			gene=t(sub_data)
			)
		colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

		pseudotime_data <- rbind(pseudotime_data, data)

	}


	#
	# Old plot I keep for record in case I need it
	#

	# p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	# scale_shape_manual(values = 21:25) +
	# geom_point(size = 3, stroke=0.7, aes(shape=condition1, fill=condition2), color="white", na.rm = TRUE)+
	# geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
	# ylab("log(RPKM+1)") +
	# theme_bw() +
	# ggtitle(geneName)+
	# # scale_color_manual(
	# # 	values=colours
	# # ) +
	# scale_fill_manual(
	# 	values=colours
	# ) +
	# guides(
	# 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
 #    ) +
	# theme(
	# 	axis.text=element_text(size=16),
	# 	axis.title=element_text(size=16),
	# 	legend.text = element_text(size =16),
	# 	legend.title=element_blank(),
	# 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	# 	aspect.ratio=0.5#,
	# 	# legend.position="bottom"
	# ) +
	# facet_grid(.~lineage) +
	# theme(
	# 	strip.text.x = element_text(size = 16)
	# )

	# pseudotime_data[pseudotime_data=="Lineage 1"] <- "Sertoli lineage"
	# pseudotime_data[pseudotime_data=="Lineage 2"] <- "Leydig lineage"
	# pseudotime_data[pseudotime_data=="Lineage 3"] <- "Progenitor lineage"

	# p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	# # geom_point(shape=21, size = 3, stroke=0.7, aes(fill=clusters), color="white", na.rm = TRUE)+
	# # geom_point(shape=108, size = 8,  aes(y=-1, color=stages), na.rm = TRUE)+
	# # geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
	# geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
	# ylab("log(RPKM+1)") +
	# theme_bw() +
	# ggtitle(geneName)+
	# scale_color_manual(
	# 	values=cluster_colors
	# ) +
	# scale_fill_manual(
	# 	values=cluster_colors
	# ) +
	# # ylim(0, max(pseudotime_data$gene))+
	# # guides(
	# # 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
 # #    ) +
	# theme(
	# 	axis.text=element_text(size=16),
	# 	axis.title=element_text(size=16),
	# 	legend.text = element_text(size =16),
	# 	legend.title=element_blank(),
	# 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	# 	aspect.ratio=0.5,
	# 	legend.position="bottom"
	# ) +
	# # facet_grid(.~lineage) +
	# theme(
	# 	strip.text.x = element_text(size = 16)
	# )

	p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
	geom_point(shape=21, size = 2.5, stroke=0.5, aes(fill=clusters), color="white", na.rm = TRUE)+
	geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE)+
	geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
	ylab("log(RPKM+1)") +
	theme_bw() +
	ggtitle(geneName)+
	scale_color_manual(
		values=stage_colors
	) +
	scale_fill_manual(
		values=cluster_colors
	) +
	theme(
		axis.text=element_text(size=16),
		axis.title=element_text(size=16),
		legend.text = element_text(size =16),
		legend.title=element_blank(),
		plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
		aspect.ratio=0.8,
		legend.position="bottom"
	) +
	facet_grid(.~lineage) +
	theme(
		strip.text.x = element_text(size = 16)
	)
	print(p)
}


compare_lineage <- function(pseudotime=pseudotime, condition=condition){
	myLineages <- list()
	for (lineages in 1:length(pseudotime))local({
		i <- lineages
		lineage <- pseudotime[[lineages]]$pseudotime
		myLineages[[i]] <<- lineage
	})
	df <- t(do.call("rbind",myLineages))
	df[df>=0] <- 1
	df[is.na(df)] <- 0
	sdf <- as.character(rowSums(df))
	return(sdf)

}


###########################################
#                                         #
#            Loess smoothing              #
#                                         #
###########################################

smooth_gene_exp <- function(data=data, pseudotime=pseudotime, span=0.75){
	smooth_data <- data
	for (gene in 1:nrow(data)){
		gene_exp <- t(data[gene,])
		smooth <- loess(formula=gene_exp~pseudotime, span=span)
		smooth_data[gene,] <- predict(smooth, newdata=pseudotime)
	}
	return(smooth_data)
}


get_var_genes_pseudotime <- function(rpkm, count, pseudotime, lineageNb, clusters){
	lineage_pseudotime_all <- pseudotime[,lineageNb]
	lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb]
	rpkm_exp_lineage <- rpkm[,names(lineage)]
	count_exp_lineage <- count[,names(lineage)]
	clusters_in_lineage <- clusters[names(clusters) %in% names(lineage)]

	print("Prepare for DE...")
	genes_pseudoT <- prepare_for_DE (
		count, 
		lineage_pseudotime_all, 
		lineage_pseudotime_all
	)
	print("Done.")

	genes_pseudoT$Pseudotime <- lineage_pseudotime_all

	print("Compute DE genes...")

	#Select cells from the lineage of interest
	genes_pseudoT <- genes_pseudoT[,names(lineage)]
	# Re-detect how many cells express each genes in the subset of cells
	genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
	# Remove genes expressed in less than 10 cells
	genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]

	DE_genes_pseudoT <- differentialGeneTest(
		genes_pseudoT, 
		fullModelFormulaStr="~sm.ns(Pseudotime, df=4)",
		cores = 3
	)
	print("Done.")

	return(DE_genes_pseudoT)
}



get_gene_clustering <- function(DE_genes_pseudoT, rpkm, pseudotime, lineageNb, cell_clusters, clusterNb=10, qvalue=0.1){
	lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb]
	clusters_in_lineage <- cell_clusters[names(cell_clusters) %in% names(lineage)]

	print(paste("Select genes with q-value<", qvalue, "...", sep=""))
	sig_gene_pseudoT <- row.names(subset(DE_genes_pseudoT, qval < qvalue))
	sig_var_gene_pseudoT <- rpkm[sig_gene_pseudoT,names(lineage)]

	print(paste(length(sig_gene_pseudoT), " selected genes.", sep=""))

	print("Smooth gene expression...")
	smooth <- smooth_gene_exp(
		sig_var_gene_pseudoT, 
		lineage, 
		span=0.5
	)
	print("Done.")

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58', '#081d58'))
	warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
	mypalette <- c(rev(cold(20)), warm(20))

	smooth <- smooth[ , order(lineage)]

	print("Cluster the genes...")

	cell_cluster_palette <- c(
		"#457cff", 
		# "#8ae400", 
		"#00c5ec", 
		"#009900", 
		"#8dfaff"
	)

	cell_cluster_colors <- cell_cluster_palette[1:length(unique(clusters_in_lineage))]
	names(cell_cluster_colors) <- unique(clusters_in_lineage)

	gene_clustering <- pheatmap(
		smooth, 
		scale="row", 
		cutree_rows=clusterNb,
		clustering_method="ward.D2",
		silent=TRUE
	)

	print("Plot heatmap...")

	clusters <- cutree(gene_clustering$tree_row, k = clusterNb)
	clustering <- data.frame(clusters)
	clustering[,1] <- as.character(clustering[,1])
	colnames(clustering) <- "Gene_Clusters"

	gene_cluster_palette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#49beaa', '#611c35', '#2708a0')
	gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
	names(gene_cluster_colors) <- 1:max(clusters)


	annotation_row <- clustering

	annotation_col <- data.frame(
		Cell_Clusters=clusters_in_lineage,
		Stages=sapply(strsplit(names(clusters_in_lineage), "_"), `[`, 1)
	)

	annotation_colors <- list(
		Stages=c(
			# P6="#d73027", 
			E10.5="#2754b5", 
			E11.5="#8a00b0", 
			E12.5="#d20e0f", 
			E13.5="#f77f05", 
			E16.5="#f9db21"
		),
		Cell_Clusters=cell_cluster_colors,
		Gene_Clusters=gene_cluster_colors
	)

	pheatmap(
		smooth, 
		scale="row", 
		cluster_cols=FALSE, 
		show_colnames=FALSE, 
		show_rownames=FALSE, 
		cutree_rows=clusterNb, 
		clustering_method="ward.D2",
		annotation_col=annotation_col,
		annotation_row=annotation_row,
		annotation_colors=annotation_colors,
		annotation_names_row=FALSE,
		color=mypalette
	)

	return(clustering)
}
