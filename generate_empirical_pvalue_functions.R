### generate empirical pvalue functions  ###
## final update : 2021-11 -30
##=========================================================================

library(dplyr)
library(plyr)
library(parallel)
library(DESeq2)
library(edgeR)
library(reshape2)
library(data.table)

#zero count filtering 
filtering <- function(countmatrix, filtering_num){
	print(paste0("zero count: ", nrow(countmatrix[rowSums(countmatrix)==0,])))
	countmatrix <- filter(countmatrix, apply(countmatrix, 1, mean) > filtering_num)
	countmatrix <- countmatrix[order(apply(countmatrix, 1, mean)),,drop = FALSE]

	return(countmatrix)
}

### create pseudo-replicates in window
#only for DMSO model 
createMetagene <- function(window){
	
	#window_1 = unlist(window[,c(1,3,5)])
	#window_2 = unlist(window[,c(2,4,6)])

	metagene_for_each_gene <- lapply(1:iter, function(stage){
		
		window <- unlist(window)
		set.seed(stage)
		sample = t(data.frame(sample(window, sampling*2, replace = TRUE)))
		rownames(sample) <- stage
		colnames(sample) <- c(paste0(rep('s1_', sampling), seq(1,sampling)), paste0(rep('s2_', sampling), seq(1,sampling)))

		return(sample)})

		metagene_for_each_genedf = as.data.frame(do.call(rbind, metagene_for_each_gene))

		return(metagene_for_each_genedf)
}

### DEG analysis with pseudo-replicates of DMSO sample (control) - DESeq2 / edgeR
#DESEq code
deseq_metagenes <- function(metagene_for_each_genedf){
	coldata <- data.frame(name = colnames(metagene_for_each_genedf), type = as.factor(c(rep(1, sampling), rep(2, sampling))))
	dds <- DESeqDataSetFromMatrix(countData = metagene_for_each_genedf, colData = coldata, design = ~ type)
	dds <- estimateSizeFactors(dds)
	#dds <- estimateDispersions(dds, quiet = TRUE, fitType = 'local')
	dds <- estimateDispersionsGeneEst(dds)
	ddsfit <- try(estimateDispersionsFit(dds , fitType = 'mean', quiet = TRUE), silent = TRUE)
	if(class(ddsfit) == "try-error") { 		
		dispersions(dds) <- mcols(dds)$dispGeneEst
		dds <- nbinomWaldTest(dds)
		result <- results(dds)
		return (result)
	}

	else{
		dds <- estimateDispersionsMAP(ddsfit)
		dds <- nbinomWaldTest(dds)
		result <- results(dds)
		return(result)
	}
}

deseq_DMSO <- function(countmatrix, window_size){

	w = window_size

	metagene_pval_list <- mclapply(1:nrow(countmatrix),  mc.cores = core, function(xx){

		if (xx %% 1000 == 1){
			print(paste0((xx-1), "th gene processing"))
		}

		if (xx - w > 0 && xx + w <= nrow(countmatrix)){ window = countmatrix[(xx - w):(xx + w), ] }
		else if (xx - w <=0) {	window = countmatrix[1:9, ] }
		else if (xx + w > nrow(countmatrix)) { window = countmatrix[(nrow(countmatrix)-8):nrow(countmatrix), ] }
		
		metagene_for_each_genedf <- createMetagene(window)
		result <- deseq_metagenes(metagene_for_each_genedf)
		result$pvalue[is.na(result$pvalue)] <- 1
		deseqpval <- sort(result$pvalue)

		return(deseqpval)
	})
	return(metagene_pval_list)
}

#edgeR code
edger_metagenes <- function(metagene_for_each_genedf){
	d <- DGEList(counts=metagene_for_each_genedf, group=as.factor(c(rep(1, sampling), rep(2, sampling))))
	group <- as.factor(c(rep(1, sampling), rep(2, sampling)))
	design <- model.matrix(~0+group)
	d <- calcNormFactors(d, method='TMM')
	d1 <- estimateGLMCommonDisp(d,design)
	d1 <- estimateGLMTrendedDisp(d1,design, method="power")
	d1 <- estimateGLMTagwiseDisp(d1,design)
	et12 <- exactTest(d1)

	return(et12)
}

edger_DMSO <- function(countmatrix, window_size){

	w = window_size

	metagene_pval_list <- mclapply(1:nrow(countmatrix),  mc.cores = core, function(xx){

		if (xx %% 1000 == 1){ print(paste0((xx-1), "th gene processing")) }
		if (xx - w > 0 && xx + w <= nrow(countmatrix)){ window = countmatrix[(xx - w):(xx + w), ] }
		else if (xx - w <=0) {	window = countmatrix[1:9, ] }
		else if (xx + w > nrow(countmatrix)) { window = countmatrix[nrow(countmatrix)-8:nrow(countmatrix), ] }
		
		metagene_for_each_genedf <- createMetagene(window)
		et12 <- edger(metagene_for_each_genedf)
		edgerpval <-  sort(et12[[1]]$PValue)
		return(edgerpval)
	})
	return(metagene_pval_list)
}


# conversion of pvalue into empirical pvalue

nbin = 100

DMSO_pval_to_empiricalpval <- function(pval_list){

	empirical_pval_per_bindf <- mclapply(1:length(pval_list), mc.cores = core, function(xx){
		pval_df <- pval_list[[xx]]
		pval_df <- do.call(rbind, pval_df); colnames(pval_df) <- c(1:iter) 

		#bin by count 
		group <- rep(c(1:nbin), each = round(nrow(pval_df)/nbin)); group <- group[1:nrow(pval_df)]
		names(group) <- rownames(pval_df)
		group_quantile <- data.frame(group = group, percentile = seq(from = 1, to = nbin, length = nrow(pval_df)))
		pval_df <- cbind(pval_df, group_quantile)

		#split by group(bin)
		pval_perBin <- split(pval_df[,1:iter], pval_df$group)

		pval_perBin <- lapply(1:length(pval_perBin), function(xx){
			pval_perBin <- unname(unlist(pval_perBin[[xx]]))
			return(pval_perBin)
		})
		names(pval_perBin) <- c(1:nbin)

		#change into log2 scale bc of loess 
		empirical_pval_per_bin <- lapply(1:length(pval_perBin), function(xx){
			temp <- unique(pval_perBin[[xx]]) 
			temp <- data.frame(original_pval = log2(sort(temp, decreasing = FALSE)))
			temp$empirical_pval <- log2(c(1:nrow(temp)) / nrow(temp))
			temp$group <- xx
			return(temp)
		})
		names(empirical_pval_per_bin) <- c(1:nbin)

		empirical_pval_per_bindf <- do.call(rbind, empirical_pval_per_bin)
		return(empirical_pval_per_bindf)
		})
	
	return(empirical_pval_per_bindf)
}


# loess ; empirical pvalue ~ group_quantile + original pvalue

loess_empirical <- function(emppval_perbin_list){
	loesslist <- mclapply(1:length(emppval_perbin_list), mc.cores = core, function(xx){
		emppval_perbin <- emppval_perbin_list[[xx]]
		attach(emppval_perbin)
		loess <- loess(empirical_pval ~ original_pval + group)
		detach(emppval_perbin)

		return(loess)
	})
}

# Drug data 
get_drug_list <- function(list){
	temp <- strsplit(names(list),'[.]')
	drug_list <- unlist(lapply(temp, function(xx) return(xx[[3]])))
	drug_list <- drug_list[-grep('Niclo', drug_list)]
	return(drug_list)
}


deseq_drug <- function(drug_count, version, dose){
	
	#version = plate num
	id = paste0(version, '.', dose)
	DMSO_count <- DMSO_count_list[[id]]
	countdata0 <- cbind(DMSO_count, drug_count)
	#countdata0 <- countdata0[order(apply(countdata0, 1, mean)),,drop = FALSE]
	coldata <- data.frame(name = colnames(countdata0), strain = as.factor(c(rep(0, ncol(DMSO_count)), rep(1, ncol(drug_count)))))
	dds <- DESeqDataSetFromMatrix(countData = countdata0, colData = coldata, design = ~ strain)
	ddsres <- DESeq(dds, quiet = TRUE)
	result <- results(ddsres)
	result <- data.frame(result)
	colnames(result) <- c("baseMean", "logFC", "logFCSE", "stat", "pvalue", "padj")
	result <- result[,c("logFC", "pvalue", "padj", "stat")]

	#binning
	DMSO_count0 <- DMSO_count[order(apply(DMSO_count, 1, mean)),,drop = FALSE]
	gene_order <- rownames(DMSO_count0)
	group <- rep(c(1:nbin), each = round(length(gene_order)/nbin))
	group <- group[1:length(gene_order)]
	names(group) <- gene_order
	group_quantile <- data.frame(group = group, percentile = seq(from = 1, to = nbin, length = length(gene_order)))

	result <- result[gene_order,]
	result <- cbind(result, group_quantile)


	return(result)
}

edger_drug <- function(drug_count, version, dose){

	#version = plate num
	id = paste0(version, '.', dose)
	DMSO_count <- DMSO_count_list[[id]]
	countdata0 <- cbind(DMSO_count, drug_count)
	#countdata0 <- filtering(countdata0, 1)
	coldata <- data.frame(name = colnames(countdata0), strain = as.factor(c(rep(0, ncol(DMSO_count)), rep(1, ncol(drug_count)))))
	d <- DGEList(counts=countdata0, group = coldata$strain)
	d <- calcNormFactors(d, method='TMM')
	d1 <- estimateCommonDisp(d)
	d1 <- estimateTagwiseDisp(d1)
	et12 <- exactTest(d1)
	et12[[1]]$padj <- p.adjust(et12$table$PValue, method='BH')
	result <- et12[[1]]; colnames(result) <- c("logFC", "logCPM", "pvalue", "padj")
	result <- result[,c("logFC", "pvalue", "padj")]

	#binning
	DMSO_count0 <- DMSO_count[order(apply(DMSO_count, 1, mean)),,drop = FALSE]
	gene_order <- rownames(DMSO_count0)
	group <- rep(c(1:nbin), each = round(length(gene_order)/nbin))
	group <- group[1:length(gene_order)]
	names(group) <- gene_order
	group_quantile <- data.frame(group = group, percentile = seq(from = 1, to = nbin, length = length(gene_order)))

	result <- result[gene_order,]
	result <- cbind(result, group_quantile)

	return(result)
}


#degres_list = origianl deg analysis result DMSO vs drug
#loessres_list = DMSO loess result
#pval_list = input of loess
#group_quantile_list = group_quantile 

group_quantile <- function(pval_list) {
	group_quantile <- mclapply(1:length(pval_list), mc.cores = core, function(xx){
		#loess result per drugs 
		pval <- pval_list[[xx]]
		pval_df <- do.call(rbind, pval); colnames(pval_df) <- c(1:iter) 

		#nbin = bin size
		nbin = 100

		#bin by count 
		group <- rep(c(1:nbin), each = round(nrow(pval_df)/nbin)); group <- group[1:nrow(pval_df)]
		names(group) <- rownames(pval_df)
		group_quantile <- data.frame(group = group, percentile = seq(from = 1, to = nbin, length = nrow(pval_df)))
		
		return(group_quantile)
	})
}


predict_empiricalPvalue <- function(degres_list, loessres_list, pval_list, group_quantile_list){
	empres_list <- mclapply(1:length(degres_list), mc.cores = core, function(xx){
			degres <- degres_list[[xx]]
			plate <- paste0(strsplit(names(degres_list)[[xx]], '[.]')[[1]][1:2], collapse = '.')
			loessres <- loessres_list[[plate]]
			degres[setdiff(rownames(degres), names(pval_list[[plate]])),]$percentile <- 1
			degres[rownames(group_quantile_list[[plate]]),]$percentile <- group_quantile_list[[plate]]$percentile
			newdata <- data.frame(group = degres$percentile, original_pval = log2(degres$pvalue))	
			degres$empirical_pval <- 2^(predict(loessres, newdata = newdata))
			return(degres)
		})
	names(empres_list) <- names(degres_list)
	
	#extrapolation problem solution
	empres_list <- lapply(empres_list, function(xx){
		xx <- xx[!is.na(xx$pvalue),]
		if(sum(is.na(xx$empirical_pval)) > 0){
			xx[is.na(xx$empirical_pval),]$empirical_pval <- 1e-06
		xx[which(xx$empirical_pval > 1),]$empirical_pval <- 1
		}
		return(xx)
	})

	return(empres_list)


}