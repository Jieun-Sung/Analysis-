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

#1. DEG Analysis for each group (SS, NASH, FIBROSIS, HCC)

#1) File reading
#agilent microarray 

FileName <- list.files("data/", pattern = "GSM")
project2 <- read.maimages(FileName, source="agilent.median", 
                          columns=list(G="gMedianSignal", Gb="gBGMedianSignal"),
                          #path="/media/sung/J'S USB/SPL/GSE/DEGs/GSE83596/data",
                          path="/spstorage/INTERNSHIP/sung/GSE83596/data")

#Affymetrix microarray
getGEOSuppFiles("GSE45327")
setwd("D:/SPL/GSE/DEGs/GSE45327")
untar("GSE45327_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep = "/"), gunzip)
exp.cel_45327 <- read.celfiles(list.celfiles("data/", full.names = T))

GSE83596<-getGEO('GSE83596', GSEMatrix=TRUE)
GSE83596<-GSE83596[[1]]
eset<-exprs(GSE83596)

exp.cel_83596 <- project2 #raw data
exp.cel_83596 <- exp.cel_83596[,1:28]
exp.rma_83596 <- normalizeBetweenArrays(exp.cel_83596) #agilent file
exp.rma_45327 <- rma(exp.cel_45327) #affymetrix

#MDS / PCA dimension reduction to visualize the differentiation of the data
plotMDS(exp.rma_83596)

# data frame
eset83596 <- exprs(GSE83596)
eset83596 <- data.frame(eset83596)
eset83596 <- eset83596[, 1:28]
colnames(eset83596)
pset83596 <- pData(GSE83596)
pset83596 <- pset83596[1:28, ]
rownames(pset83596)
fset83596 <- fData(GSE83596)

#linear regression 
eset83596$REFSEQ<- fset83596$REFSEQ
eset83596 <- aggregate(. ~REFSEQ, data = eset83596, mean)
eset83596 <- eset83596[, -1]
colnames(eset83596)

grp_83596_1 <- pset83596$`disease state:ch1`
grp_83596_2 <- pset83596$`age:ch1`
grp_83596 <- paste(grp_83596_1, grp_83596_2)

design_83596 <- model.matrix(~0 + grp_83596)
nchar(colnames(design_83596))

fit_83596 <- lmFit(exp.rma_83596, design_83596)

#2. Find DEG

cont_83596 <- makeContrasts(steatosis6weeks - ctrl6weeks ,levels=design_83596)
fit.cont3 <- contrasts.fit(fit_83596,cont_83596)
fit.cont3 <- eBayes(fit.cont3)
res_83596_1 <- topTable(fit.cont3, number=Inf,lfc=0,adjust="BH")
up_res83596_1 <- res_83596_1[res_83596_1$logFC > 1.2 & res_83596_1$adj.P.Val < 0.05,]
nrow(up_res83596_1)  #133
down_res83596_1 <- res_83596_1[res_83596_1$adj.P.Val < 0.05 & res_83596_1$logFC < -1.2,]
nrow(down_res83596_1) #18

EnhancedVolcano(res_83596_1,
                lab = rownames(res_83596_1),
                x = 'logFC',
                y = 'P.Value',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)

#3. Pathway enrichment analysis using GO

# convert to ENTREZ id
up_res83596_1_E <- bitr(up_res83596_1$SystematicName, fromType="REFSEQ", toType="ENTREZID", OrgDb="org.Mm.eg.db")
down_res83596_1_E <- bitr(down_res83596_1$SystematicName, fromType="REFSEQ", toType="ENTREZID", OrgDb="org.Mm.eg.db")

write.csv(up_res83596_1_E,file="//media/sung/J'S USB/SPL/GSE/DEGs/GSE83596/up_83596_6wks_id.csv")
write.csv(down_res83596_1_E,file="//media/sung/J'S USB/SPLGSE/DEGs/GSE83596/down_83596_6wks_id.csv")

go_83596_1 <- enrichGO(up_res83596_1_E$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all", pvalueCutoff = 0.05) #keyType ='SYMBOL',pAdjustMethod = "BH",
go1_83596_1 <- enrichGO(down_res83596_1_E$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all", pvalueCutoff = 0.05)

nrow(go_83596_1)  #87
nrow(go1_83596_1) #133

barplot(go_83596_1, showCategory=10, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale = "free")
barplot(go1_83596_1, showCategory=10, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale = "free")

edox_83596 <- setReadable(go_83596_1, "org.Mm.eg.db", 'ENTREZID')
edox1_83596 <- setReadable(go1_83596_1, "org.Mm.eg.db", 'ENTREZID')

p_83596_1 <- cnetplot(edox_83596, categorySize="pvalue", foldChange=up_res83596_1$ENTREZID, circular = FALSE, colorEdge = TRUE)
p1_83596_1 <- cnetplot(edox1_83596, categorySize="pvalue", foldChange=down_res83596_1$ENTREZID, circular = FALSE, colorEdge = TRUE)

h_83596_1 <- heatplot(edox_83596, foldChange=up_res83596_1$ENTREZID)
h1_83596_1 <- heatplot(edox1_83596, foldChange=down_res83596_1$ENTREZID)

cowplot::plot_grid(p_83596_1)
cowplot::plot_grid(p1_83596_1)

cowplot::plot_grid(h_83596_1)
cowplot::plot_grid(h1_83596_1)


#3. Clustering 
#ordering dataframe
res_83596_1<- res_83596_1 %>% arrange(rownames(res_83596_1)) %>% arrange(SystematicName)
res_83596_2<- res_83596_2 %>% arrange(rownames(res_83596_2)) %>% arrange(SystematicName)
res_83596_3<- res_83596_3 %>% arrange(rownames(res_83596_3)) %>% arrange(SystematicName)
res_83596_4<- res_83596_4 %>% arrange(rownames(res_83596_4)) %>% arrange(SystematicName)

#grouping by logFC
dfl = lapply(list(res_83596_1,res_83596_2,res_83596_3,res_83596_4), function(xx){
  idx = abs(xx$logFC)> 1 & xx$adj.P.Val < 0.05
  xx1 = xx$logFC
  xx1[!idx] = 0
  xx1 = ifelse(xx1 > 0, 1, ifelse(xx1 < 0,-1, 0))
})

dfl = do.call(cbind, dfl)
dfl <- data.frame(dfl)
dfl$REFSEQ <- res_83596_1$SystematicName
dfl$Order <- rownames(res_83596_1)
colnames(dfl) <- c("SS", "NASH", "Fibrosis", "HCC", "REFSEQ", "Order")

# make groups and assign
dfg <- unique(dfl[,1:4])
dfg <- dfg %>% arrange(desc(HCC)) %>% arrange(desc(Fibrosis)) %>% arrange(desc(NASH)) %>% arrange(desc(SS))
    

#ncl <- df_cl %>% group_by(df_cl$cluster) %>% summarise(n = n())
ncl <- table(df_cl$cluster)
scl <- rownames(subset(ncl, ncl > 10))
s_cl <- rownames(subset(ncl, ncl <= 10))

df_cl[(which(df_cl$cluster %in% s_cl)),'cluster'] <- "class 20"  #class20 = 0000

df_hm <- subset(df_cl, cluster!="class 20") #4745
df_hm <- df_hm %>% arrange(REFSEQ) %>% arrange(Order)
resl <- unique(df_hm$Order)


# Heatmap drawing 
fc_cl <- data.frame(SS=subset(res_83596_1, rownames(res_83596_1) %in% resl)$logFC, 
                    NASH=subset(res_83596_2, rownames(res_83596_2) %in% resl)$logFC, 
                    Fibrosis=subset(res_83596_3, rownames(res_83596_3) %in% resl)$logFC, 
                    HCC=subset(res_83596_4, rownames(res_83596_4) %in% resl)$logFC, 
                    cluster=df_hm$cluster, REFSEQ=df_hm$REFSEQ)
 #only significant DEGs

# colnames(fc_hm) <- c("SS", "NASH", "Fibrosis", "HCC", "cluster", "REFSEQ")

Heatmap(as.matrix(fc_cl[1:4]), show_row_names = F,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", 'red')), 
        cluster_columns = F, cluster_rows = TRUE,
        column_order = c('SS', 'NASH', 'Fibrosis', 'HCC'),
        row_title = "Genes",
        row_split = fc_hm$cluster,
        border=T,
        clustering_distance_rows = function(x, y) 1 - cor(x, y),
        clustering_method_rows = "complete")


#4. Pathway analysis using human ID / cp, cgp. kegg database

#hypergeometric test
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



cl <- sort(unique(fc_cl$cluster)) #16개 -> 22개 

#cl_list : expression list with REFSEQ
cl_list <- list()
for ( i in 1:length(cl)){
    xx = cl[i]
    cll1 <- subset(fc_cl, fc_cl$cluster ==  cl[i])
    cl_list <- append(cl_list, list(xx = cll1))
}
names(cl_list) <- cl

#convert REFSEQ id into ENTREZ 
#cl_list : expression data with refseq and clustering
#entl : cl_list with entrez id 

load("/spstorage/USERS/gina/Project/2B4/m2h.RData")
m2h <- rename(m2h, "EntrezGene.ID" = "ENTREZID")

entl1 <- mclapply(names(cl_list), mc.cores = 20, function(xx){

  cl <- cl_list[xx]
  nmx <- bitr(cl[[1]]$REFSEQ, fromType="REFSEQ", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  cl <- merge(cl[[1]], nmx, by = "REFSEQ", all.x = T)
  cl <- na.omit(cl)
  ent_h <- m2h[which (m2h$ENTREZID %in% cl$ENTREZID),]
  cl <- merge(cl, ent_h, by = "ENTREZID", all.x = T)
  cl <- na.omit(cl)

  res <- data.frame( h.entrez = cl$entrezid.h)

  return(res)
  })

names(entl1) <- names(cl_list)


library(GSA)
cgpgmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cgp.v6.2.entrez.gmt")
cpgmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cp.v7.2.entrez.gmt")
kegggmt <- GSA.read.gmt("/spstorage/DB/MSigDB/c2.cp.kegg.v6.2.entrez.gmt") 

gmtdf <- cgpgmt[[1]]
names(gmtdf) <- cgpgmt[[2]]

cpdf <- cpgmt[[1]]
names(cpdf) <- cpgmt[[2]]


# cgp
gspace <- intersect(unlist(gmtdf), unlist(entl))

wholepathres <- mclapply(names(entl), mc.cores=20, function(cl){

  print (paste(cl, "th cluster processing"))
  clx <- entl[cl][[1]]
  clxx <- intersect(clx$h.entrez, gspace)
  result <- hypergeoTestForGeneset(clxx, gmtdf, gspace)
  result$type = cl
  return (result)
})
names(wholepathres) <- names(cl_list)

#cp
gspace <- intersect(unlist(cpdf), unlist(entl1))

cpwholepathres <- mclapply(names(entl1), mc.cores=20, function(cl){

  print (paste(cl, "th cluster processing"))
  clx <- entl[cl][[1]]
  clxx <- intersect(clx$h.entrez, gspace)
  result <- hypergeoTestForGeneset(clxx, cpdf, gspace)
  result$type = cl

  return (result)
})
names(cpwholepathres) <- names(cl_list)

sigpathres <- lapply(names(cpwholepathres), function(cl){

  pathx <- cpwholepathres[cl][[1]]
  pathx <- pathx[which(pathx$adjpVal < 0.05 & pathx$int > 2),]
 
  return(pathx)
  })
names(sigpathres) <- names(cl_list)

sigpathresdf = do.call(rbind, sigpathres)

pathhm<- lapply(names(cpwholepathres), function(cl){

  pathx <- cpwholepathres[cl][[1]]
  pathx <- pathx[which(pathx$ID %in% sigpathresdf$ID),]
 
  return(pathx)
})
names(pathhm) <- names(cl_list)

#nrow of sigpathresdf = 356

#heatmap drawing \
pathresdf = do.call(rbind, pathhm)
prm = dcast(pathresdf, ID ~ type, value.var='adjpVal')
rownames(prm) <- prm$ID ; prm <- prm[-1]

pix = apply(prm, 2, function(x){

  ixix = which(x < 0.05)
  length(ixix) >=2 

  })
prm1 <- prm[,pix]

pix1 = apply(prm, 1, function(x){

  ixix = which(x < 0.05)
  # ixix = setdiff(ixix, c(11,13))
  length(ixix) >= 2

  })
# pix1 = names(sort(pix1, decreasing=T))

prm1 <- prm[pix1,]

#prm : pvalue for each pathways 
# prm[is.na(prm)] = 1
grp2 <- dfg[which(rownames(dfg) %in% colnames(prm)),]
grp2 <- dfg[which(rownames(dfg) %in% colnames(prm1)),]
#ordering grp2 by -1,0,1
grp2[grp2==0] = 'none'
grp2[grp2==1] = 'up'
grp2[grp2==-1] = 'down'
grp2 <- t(grp2)

grp2 <- grp2[, c("class 4","class 1", "class 11", "class 10", "class 17","class 19","class 21", "class 24")]

prm2 <- prm1[,which(colnames(prm1) %in% colnames(grp2))]

antup <- Heatmap(grp2, 
        name='DEG',
        row_order = c('SS', 'NASH', 'Fibrosis', 'HCC'), 
        height = unit(2, "cm"),
        border=T,
        rect_gp = gpar(col = "black", lwd  = 0.2),
        cluster_columns = F, column_order = colnames(grp2), 
        col = structure(c('blue','white','red'), names=c('down','none','up'))
        )

hmup <- Heatmap(-log10(as.matrix(prm2)), show_row_names = F,
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

hmtotalup <- antup %v% hmup
draw(hmtotalup, heatmap_legend_side='left')



#5. Chip-Seq analysis using ENCODE / remap database
encodegmt <- GSA.read.gmt("/spstorage/USERS/hanbi/otherdata/ENCODE_ChIP-seq.gmt")
remapgmt <- GSA.read.gmt("/spstorage/USERS/hanbi/otherdata/ReMap_ChIP-seq.gmt")
encodedf <- encodegmt[[1]]
names(encodedf) <- encodegmt[[2]]

remapdf <- remapgmt[[1]]
names(remapdf) <- remapgmt[[2]]

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

tf<- mclapply(names(tfres1), mc.cores = 20, function(cl){

  tfx <- tfres1[cl][[1]]
  tfx <- tfx[which(tfx$ID %in% sigtfresdf$ID & tfx$type %in% colnames(grp2)),]
 
  return(tfx)
})
names(tf) <- names(tfres1)

tfdf = do.call(rbind, tf)
tfhm = dcast(tfdf, ID ~ type, value.var='adjpVal')
rownames(tfhm) <- tfhm$ID ; tfhm <- tfhm[-1]
tfhm <- tfhm[, colnames(grp2)]

pix1 = apply(tfhm, 1, function(x){

  ixix = which(x < 0.05)
  # ixix = setdiff(ixix, c(11,13))
  length(ixix) >= 1

  })
# pix1 = names(sort(pix1, decreasing=T))

tfhm1 <- tfhm[pix1,]

grpuphm <- Heatmap(grp2, 
        name='DEG',
        row_order = c('SS', 'NASH', 'Fibrosis', 'HCC'), 
        height = unit(2, "cm"),
        rect_gp = gpar(col = "black", lwd  = 0.2),
        border=T,
        cluster_columns = F, column_order = colnames(grp2), 
        col = structure(c('blue','white','red'), names=c('down','none','up'))
        )

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

tfuphml <- grpuphm %v% tfuphmm

draw(tfuphml, row_title_gp = gpar(fontsize = 0.2), heatmap_legend_side='left', padding = unit(c(0,20,10,50),'points'))

