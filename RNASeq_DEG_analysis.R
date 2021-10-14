library(GEOquery)
library(limma)
library(affy)
library(oligoClasses)
library(oligo)
library(EnhancedVolcano)
library(affyio)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(stringr)
library(VennDiagram)
library(gplots)
library("RColorBrewer")
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(devtools)

X11.options(type='cairo')


##################################################################################################################################################
# load data
##################################################################################################################################################

gtf_file = "/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf"

read_featurecounts <- function(output_dir, mc.cores){
	files <- list.files(output_dir)

	gtf = rtracklayer::import(gtf_file) %>% as.data.frame()
	pcg = gtf %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
	pcg$gene_id = gsub("\\.[0-9]+","",pcg$gene_id)  #ID만 가져오기 

	#check on path
	path = paste0(output_dir,"/",file,"/",file,".featureCounts.tsv")

	count.expr.list <- mclapply(1:length(files), mc.cores = mc.cores, function(xx){
    	file = files[xx]
    	countmatrix = read.table(path, header = TRUE)
    	countmatrix$Geneid <- gsub("\\.[0-9]+","",countmatrix$Geneid)
    	countmatrix.pcg <- filter(countmatrix, countmatrix$Geneid %in% pcg$gene_id)
    
   	 	return (countmatrix.pcg)}
		)

	names(count.expr.list) <- files
	countdf0 = do.call(cbind, count.expr.list)

	count <-countdf0[,grep('bam.files', colnames(countdf0))] 
	rownames(count) <- countdf0$SRR6080215.Geneid
	rownames(count) <- gsub("\\.", "", rownames(count))
	colnames(count) <- files

	#filtering 
	count0 <- nrow(count[!rowSums(count) == 0,])

	return(count0)
}


raw_to_tpm <- function(count, genelength){
	rpk <- count/(genelength / 1000)
	coef <- sum(rpk) / 1e6
	tpm <- rpk/coef
	return(tpm)
}
tpms <- apply(count, 2, function(x) raw_to_tpm(x, genelength))


##################################################################################################################################################
# DEG Analysis / Choose DESeq2 or edgeR
##################################################################################################################################################

#-------------------------1. DESeq2-------------------------#
dds <- DESeqDataSetFromMatrix(countData = count, colData = coldata, design = ~ patient + strain)
ddsres <- DESeq(dds)
res0 <- results(ddsres)

res_nona <- na.omit(res0)
up_DEG <- res_nona[res_nona$log2FoldChange > log2(3) & res_nona$padj < 0.05,]
nrow(up_DEG) 

dn_DEG <- res_nona[res_nona$log2FoldChange < -log2(3) & res_nona$padj < 0.05,]
nrow(dn_DEG) 

DEG <- rbind(up_DEG, dn_DEG)

#volcano plot
with(res0, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3), col="grey"))
with(subset(res0, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res0, padj<0.05 & abs(log2FoldChange)>log2(3)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA plot
vsdata <- vst(ddsres, blind=FALSE)
plotPCA(vsdata, intgroup="strain") + geom_label(aes(label = type))

#tSNE


#-------------------------2. edgeR-------------------------#
d <- DGEList(counts=count, group = as.factor(c(rep(1, 5), rep(2, 5))))
d <- calcNormFactors(d, method='TMM')
d1 <- estimateCommonDisp(d, verbose=FALSE)
d1 <- estimateTagwiseDisp(d1)
et12 <- exactTest(d1)


design <- model.matrix(~0+group)
colnames(design) <- gsub('group', '', colnames(design))
xxx <- estimateGLMCommonDisp(x, design)
xxx <- estimateGLMTrendedDisp(xxx, design, method='power') ## method according to the number of tags: "auto", "bin.spline"(defaluttag > 200), "power"(default tag <= 200), "spline", "bin.loess".
xxx <- estimateGLMTagwiseDisp(xxx, design)

#MDS plot
plotMDS(lcpm2, labels=group, col=col.group, main="MDS plot")

#volcano plot 
vol <- cbind(et12$table$logFC, fdr); colnames(vol) <- c('logFC', 'FDR'); vol <- as.data.frame(vol)
jpeg(filename='Volcano plot : control vs NAFLD with exact.jpg', units='px', bg='white', type='cairo')
EnhancedVolcano(vol,
	 lab = rownames(et12),
	 x = 'logFC',
	 y = 'FDR',
	 title='Control vs NAFLD',
	 pCutoff=0.05,
	 FCcutoff=2)


##################################################################################################################################################
# pathway analysis 
##################################################################################################################################################

#-------------------------1. hypergeometric test-------------------------#
#query = DEG set 
#refGMT
#gspace = common gene between query and gmt 

hypergeoTestForGeneset <- function(query, refGMT, gspace) {
  if(!all(query %in% gspace)) {
  stop(paste(length(setdiff(query, gspace)),'Query items were found outside of background space. Check inputs.'))
  }
  if(length(query) == 0) stop('Query length is zero.')

  N = length(gspace) # no of balls in urn
  k = length(query) # no of balls drawn from urn (DEG no)

  enrRes = lapply(refGMT, function(refgenes) {
    q = length(intersect(refgenes, query)) # no of white balls drawn
    m = length(intersect(gspace, refgenes)) # no of white balls in urn
    I = intersect(refgenes, query)

    pVal = phyper(q-1, m, N-m, k, lower.tail = FALSE)
    odds = (q / k) / (m / N)
    jacc = q / length(union(query, refgenes))

    return(data.frame(pVal = pVal, oddsRatio=odds, tan = jacc, int=q, bg=N))
  })

  enrRes = do.call(rbind, enrRes)
  enrRes$ID = names(refGMT)
  enrRes$adjpVal = p.adjust(enrRes$pVal, method = c("BH"))
  enrRes$logP = -log10(enrRes$adjpVal)
  enrRes = enrRes[,c('ID','adjpVal', 'pVal','logP','oddsRatio','tan','int','bg')]

  return(enrRes)
}

#make required files and test 

library(GSA)
cgpgmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cgp.v6.2.entrez.gmt")
cpgmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cp.v7.2.entrez.gmt")
kegggmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cp.kegg.v6.2.entrez.gmt") 

gmtdf <- cgpgmt[[1]]; names(gmtdf) <- cgpgmt[[2]]
cpdf <- cpgmt[[1]] ; names(cpdf) <- cpgmt[[2]]

gspace <- union(unlist(gmtdf), unlist(deglist))
result <- hypergeoTestForGeneset(deglist, gmtdf, gspace)
sigpathres <- result[which(result$adjpVal < 0.05 & result$int > 2),]
sigpathresdf = do.call(rbind, sigpathres)

#heatmap drawing 
prm = dcast(sigpathresdf, ID ~ type, value.var='adjpVal')
rownames(prm) <- prm$ID ; prm <- prm[-1]

pix = apply(prm, 2, function(x){

  ixix = which(x < 0.05)
  length(ixix) >=2 

  })
prm1 <- prm[,pix]

hmup <- Heatmap(-log10(as.matrix(prm1)), show_row_names = F,
        rect_gp = gpar(col = "grey80", lwd  = 0.5),
        col = colorRamp2(c(0, 2, 4), c("white", "red", 'darkred')), 
        cluster_columns = F, cluster_rows = T,
        column_order = colnames(grp2),
        row_title = "Pathway",
        # row_split = pathres$type,
        border=T,
        #clustering_distance_rows = function(x, y) 1 - cor(x, y),
        clustering_method_rows = "complete",
        heatmap_legend_param = list( col_fun = colorRamp2(c(0, 2, 4), c("white", "red", 'darkred')), 
        title = "-log10(adj.P.Val)", 
        legend_height = unit(5,"cm"),
        heatmap_legend_side = "left"))

#-------------------------2. GO analysis-------------------------#
# convert to ENTREZ id
res83596_2_E <- bitr(up_res83596_2$SystematicName, fromType="REFSEQ", toType="ENTREZID", OrgDb="org.Mm.eg.db")
go_83596_2 <- enrichGO(res83596_2_E$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all", pvalueCutoff = 0.05) 
nrow(go1_83596_2) 

barplot(go_83596_2, showCategory=10, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale = "free")
edox_83596 <- setReadable(go_83596_2, "org.Mm.eg.db", 'ENTREZID')
p_83596_2 <- cnetplot(edox_83596, categorySize="pvalue", foldChange=up_res83596_2$ENTREZID, circular = FALSE, colorEdge = TRUE)
h_83596_2 <- heatplot(edox_83596, foldChange=up_res83596_2$ENTREZID)


##################################################################################################################################################
#Chip-Seq analysis
##################################################################################################################################################

#Encode - Chip-Seq analysis
encodegmt <- GSA.read.gmt("/spstorage/USERS/hanbi/otherdata/ENCODE_ChIP-seq.gmt")
remapgmt <- GSA.read.gmt("/spstorage/USERS/hanbi/otherdata/ReMap_ChIP-seq.gmt")
encodedf <- encodegmt[[1]] ; names(encodedf) <- encodegmt[[2]]
remapdf <- remapgmt[[1]] ; names(remapdf) <- remapgmt[[2]]

syml = bitr(unique(unlist(entl)), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")

encodegspace = intersect(unique(unlist(encodedf)), syml$SYMBOL)
remapgspace = intersect(unique(unlist(remapdf)), syml$SYMBOL)

# mouse entrez -> human entrez -> gene symbol -> hypergeometric test
tfres1 <- lapply(names(cl_list), function(cl){

  print (paste(cl, "th cluster processing"))
  clx <- cl_list[cl][[1]][, c("REFSEQ", "cluster")]
  nmx = bitr(clx$REFSEQ, fromType="REFSEQ", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  clx <- merge(clx, nmx, by = "REFSEQ")
  symbol <- m2h[which (m2h$ENTREZID %in% clx$ENTREZID),c("entrezid.h", "ENTREZID")]
  clx1 <- merge(clx, symbol, by = "ENTREZID", all.x = T)
  nmx = bitr(clx1$entrezid.h, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  colnames(nmx) <- c("entrezid.h", "SYMBOL")
  clx2 <- merge(clx1, nmx, by = "entrezid.h", all.x = T)
  clxx <- intersect(clx2$SYMBOL, encodegspace)
  result <- hypergeoTestForGeneset(clxx, remapdf, encodegspace)
  result$type = cl
  return (result)
})
names(tfres1) <- names(cl_list)


#sigtfres : Final result 
sigtfres <- lapply(names(tfres1), function(cl){

  tfx <- tfres1[cl][[1]]
  tfx <- tfx[which(tfx$adjpVal < 0.05 & tfx$int > 2),]
 
  return(tfx)
  })
names(sigtfres) <- names(tfres1)

sigtfresdf = do.call(rbind, sigtfres)

tfhm = dcast(sigtfresdf, ID ~ type, value.var='adjpVal')
rownames(tfhm) <- tfhm$ID ; tfhm <- tfhm[-1]
tfhm <- tfhm[, colnames(grp2)]

pix1 = apply(tfhm, 1, function(x){

  ixix = which(x < 0.05)
  # ixix = setdiff(ixix, c(11,13))
  length(ixix) >= 1

  })
# pix1 = names(sort(pix1, decreasing=T))

tfhm1 <- tfhm[pix1,]

tfuphmm <- Heatmap(-log10(tfhm1), show_row_names = T,
        rect_gp = gpar(col = "grey80", lwd  = 0.5),
        col = colorRamp2(c(0, 2, 4), c("white", "red", 'darkred')),
        row_title_gp = gpar(fontsize = 0.1),
        column_dend_height = unit(0.5, "cm"),
        cluster_columns = F, cluster_rows = T,
        column_order = colnames(grp2),
        column_title_gp = gpar(fontsize = 0),
        row_title = "Pathway",
        # row_split = pathres$type,
        border=T,
        clustering_distance_rows = function(x, y) 1 - cor(x, y),
        clustering_method_rows = "complete",
        heatmap_legend_param = list( col_fun = colorRamp2(c(0, 2, 4), c("white", "red", 'darkred')), 
        title = "-log10(adj.P.Val)", 
        legend_height = unit(5,"cm"),
        heatmap_legend_side = "left"))

draw(tfuphmm, row_title_gp = gpar(fontsize = 0.2), heatmap_legend_side='left', padding = unit(c(0,20,10,50),'points'))