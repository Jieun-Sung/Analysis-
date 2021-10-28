library(DESeq2)
library(edgeR)

#DEG analysis
deseq <- function(countdata, coldata){
	dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ strain)
	ddsres <- DESeq(dds, quiet = TRUE)
	result <- results(ddsres)

	return(result)
}

edger <- function(countdata, coldata){
	group <- as.factor(coldata$strain)
	d <- DGEList(counts=countdata, group=group)
	d <- calcNormFactors(d, method='TMM')
	d1 <- estimateCommonDisp(d, design, verbose=FALSE)
	d1 <- estimateTagwiseDisp(d1)
	et12 <- exactTest(d1)
	result <- et12[[1]]
	result$padj <- p.adjust(result$PValue, method='BH')

	return(result)
}


# For strict samples
deseq <- function(metagene_for_each_genedf, coldata){
	dds <- DESeqDataSetFromMatrix(countData = metagene_for_each_genedf, colData = coldata, design = ~ type)
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersions(dds, quiet = TRUE, fitType = 'local')
	dds <- nbinomWaldTest(dds)
	result <- results(dds)

	#ddsres <- DESeq(dds, quiet = TRUE, fitType = 'mean')
	#result <- results(ddsres)

	#pvalue <- data.frame(DESeqPval = sort(result$pvalue))
	#pvalue[is.na(pvalue)] = 1

	return(result)
}

#edgeR code
edger <- function(metagene_for_each_genedf){
	d <- DGEList(counts=metagene_for_each_genedf, group=as.factor(c(rep(1, sampling), rep(2, sampling))))
	group <- as.factor(c(rep(1, sampling), rep(2, sampling)))
	design <- model.matrix(~0+group)
	d <- calcNormFactors(d, method='TMM')
	d1 <- estimateCommonDisp(d, design, verbose=FALSE)
	d1 <- estimateTagwiseDisp(d1)
	et12 <- exactTest(d1)

	return(et12)
}