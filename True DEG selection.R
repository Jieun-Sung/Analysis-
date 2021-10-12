library(dplyr)
library(parallel)
library(DESeq2)
library(ggplot2)
library(Cairo)
library(ComplexHeatmap)
library(reshape2)

library(gtools)
library(gridExtra)

X11.options(type='cairo')

#####################################################################################################################
load("/spstorage/USERS/sung/DEGanalysis/OUTPUT/GSE111842/GSE111842_truedeg.Rdata")

save(list = ls(), file = "/spstorage/USERS/sung/DEGanalysis/OUTPUT/GSE104310/GSE104310_truedeg.Rdata")
#####################################################################################################################
#import data
#####################################################################################################################
output_dir = "/spstorage/USERS/sung/DEGanalysis/OUTPUT/GSE104310/expression-estimate"
files <- list.files(output_dir)

gtf = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()

#only protein coding gene 
#nrow of gtf = 2913213 / nrow of pcg = 19988
pcg = gtf %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
pcg$gene_id = gsub("\\.[0-9]+","",pcg$gene_id)  #ID만 가져오기 

count.expr.list <- mclapply(1:length(files), mc.cores = 10, function(xx){
    file = files[xx]
    path = paste0(output_dir,"/",file,"/",file,".featureCounts.tsv")
    countmatrix = read.table(path, header = TRUE)
    countmatrix$Geneid <- gsub("\\.[0-9]+","",countmatrix$Geneid)
    countmatrix.pcg <- filter(countmatrix, countmatrix$Geneid %in% pcg$gene_id)

    return (countmatrix.pcg)}
)
names(count.expr.list) <- files
countdf0 = do.call(cbind, count.expr.list)
colnames(countdf0[1])

#####################################################################################################################
#raw count to tpm 
#####################################################################################################################

#count : raw count without filtering
#count0 : raw count filtered
#tpms : raw count into tpm 

raw_to_tpm <- function(count, genelength){
  rpk <- count/(genelength / 1000)
  coef <- sum(rpk) / 1e6
  tpm <- rpk/coef
  return(tpm)
}

#GSE111842
genelength <- count.expr.list[[1]]$Length
count <-countdf0[,grep('X.spstorage.USERS.sung.DEGanalysis.OUTPUT.GSE111842.bam.files', colnames(countdf0))]
rownames(count) <- countdf0$SRR6835872.Geneid
rownames(count) <- gsub("\\.", "", rownames(count))
colnames(count) <- files

#GSE104310
genelength <- count.expr.list[[1]]$Length
count <-countdf0[,grep('X.spstorage.USERS.sung.DEGanalysis.OUTPUT.GSE104310.bam.files', colnames(countdf0))]
rownames(count) <- countdf0$SRR6080215.Geneid
rownames(count) <- gsub("\\.", "", rownames(count))
colnames(count) <- files

#coldata 
#111842
type <- c(paste0(rep('normal',6), 1:6), paste0(rep('CTC', 16), 1:16))
coldata <- data.frame(files, type)
coldata$strain <- gsub("\\d", "", coldata$type)
coldata$patient <- c(1:6, 1:16)

#104310
type <- c('normal1', 'tumor1', 'normal2', 'tumor2', 'normal3', 'tumor3', 'normal4', 'tumor4', 'tumor9', 'tumor10', 'normal5', 'tumor5', 'normal6',
  'tumor6', 'normal7', 'tumor7', 'tumor11', 'tumor12', 'normal8', 'tumor8')
coldata <- data.frame(files, type)
coldata$strain <- gsub("\\d", "", coldata$type)
coldata$patient <- c(1,1,2,2,3,3,4,4, 9, 10, 5,5,6,6,7,7,11, 12, 8,8)

#136630
type <- c(paste0(rep('tumor', 7), 1:7), paste0(rep('normal', 5), 1:5))
coldata <- data.frame(files, type)
coldata$strain <- gsub("\\d", "", coldata$type)
coldata$patient <- c(1,2,3,4,5,6,7,2,3,4,6,7)

#control count + filtering
control_SRR <- coldata[coldata$strain == 'normal',]$files
count_control <- count[,which(colnames(count) %in% control_SRR)]
count0 <- filter(count, apply(count_control, 1, mean) > 4)
tpms <- apply(count, 2, function(x) raw_to_tpm(x, genelength))
count_control <- count0[,which(colnames(count0) %in% control_SRR)]
tpm_control <- tpms[which(rownames(tpms) %in% rownames(count_control)), which(colnames(tpms) %in% control_SRR)]

#group binning sorted by mean tpm 
control.mean.tpm <- as.data.frame(apply(tpm_control, 1, mean))
colnames(control.mean.tpm) <- "meanTPM"
control.mean.tpm$Geneid <- rownames(count_control)
control.mean.tpm <- control.mean.tpm %>% arrange(Geneid) %>% arrange(meanTPM)

gene_order <- control.mean.tpm$Geneid
group <- rep(c(1:100), each = round(length(gene_order)/100))
group <- group[1:length(gene_order)]
group[is.na(group)] <- 100

#####################################################################################################################
#DESEQ - combination
#####################################################################################################################

#combination
#111842
comb <- as.data.frame(combinations(6, 3, as.character(control_SRR)))

#104310
comb <- as.data.frame(combinations(8, 4, as.character(control_SRR)))

#DESEq for each combination

num = nrow(comb)/2
cutoff = c(0.1, 0.05, 0.01)

degl <- lapply(1:length(cutoff), function (x){

    cut <- cutoff[x]

    degll <- mclapply(1:num, mc.cores = 10, function(xx){

    comb0 <- unlist(comb[xx,])
    xy = nrow(comb) + 1 - xx

    comb1 <- unlist(comb[xy,])

    coldata0 <- data.frame(control_SRR = comb0, type = 0)
    coldata1 <- data.frame(control_SRR = comb1, type = 1)

    coldata00 <- rbind(coldata0, coldata1)
    coldata000 <- arrange(coldata00, as.character(control_SRR))

    #DESeq
    dds <- DESeqDataSetFromMatrix(countData = count_control, colData = coldata000, design = ~ type)
    ddsres <- DESeq(dds)
    res <- as.data.frame(results(ddsres))

    #upDEG <- res[res$pvalue < cut & res$log2FoldChange > 0,]
    #dnDEG <- res[res$pvalue < cut & res$log2FoldChange < 0,]
    res$DEG <- 0  
    res[is.na(res$pvalue),]$pvalue <- 1
    res[which(res$pvalue < cut),]$DEG <- 1

    res <- res[gene_order, ]
    res.order <- cbind(res, group = group)

    return (res.order)
    })

    return (degll)
})
names(degl) <- cutoff

#####################################################################################################################
#DESEQ - as routine, control + treated comparison
#####################################################################################################################

dds <- DESeqDataSetFromMatrix(countData = count0, colData = coldata, design = ~ strain)
ddsres <- DESeq(dds)
res <- results(ddsres)
res <- as.data.frame(res)

originaldegl <- mclapply(1:length(cutoff), mc.cores = 3, function(xx){

    cut = cutoff[xx]
    
    #upDEG <- res[res$pvalue < cut & res$log2FoldChange > 0,]
    #dnDEG <- res[res$pvalue < cut & res$log2FoldChange < 0,]
    res$DEG <- 0  
    #res[is.na(res$pvalue),]$pvalue <- 1
    res[which(res$pvalue < cut),]$DEG <- 1

    res <- res[gene_order, ]
    res.order <- cbind(res, group = group)

    return (res.order)
    })
names(originaldegl) <- cutoff

ratio <- mclapply(1:3, mc.cores = 3, function(xx){
    ratio <- rep(0, 100)

    originaldegcut <- originaldegl[[xx]]

    for (i in 1:nrow(originaldegcut)) {
        if (originaldegcut[i,]$DEG == 1) {
            j = originaldegcut[i,]$group
            ratio[j] = ratio[j] + 1
    }}
    
    ratio <- ratio / table(group)
    ratio <- as.data.frame(ratio)

    return (ratio)
    })

ratiodf <- do.call(cbind, ratio)
ratiodf <- ratiodf[,c(1,2,4,6)]
colnames(ratiodf) <- c("group", "deg0.1", "deg0.05", "deg0.01")

ggplot(data = ratiodf) +
    geom_line(aes(x = group, y = deg0.1, group = 1), color = 'red') +
    geom_line(aes(x = group, y = deg0.05, group = 1), color = 'blue') +
    geom_line(aes(x = group, y = deg0.01, group = 1), color = 'green') +
    #geom_point() + 
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

#####################################################################################################################
#plot
#####################################################################################################################

ratiodf0.1 <- aggregate(Freq ~ group, data = ratiodf.0.1, mean)
ratiodf0.05 <- aggregate(Freq ~ group, data = ratiodf.0.05, mean)
ratiodf0.01 <- aggregate(Freq ~ group, data = ratiodf.0.01, mean)

ggplot() +
    geom_line(data = ratiodf0.1, aes(x = group, y = Freq, group = 1), color = 'black') +
    geom_line(data = ratiodf, aes(x = group, y = deg0.1, group = 1), color = 'blue') +
    labs(x="Group", y="DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

ggplot() +
    geom_line(data = ratiodf0.05, aes(x = group, y = Freq, group = 1), color = 'black') +
    geom_line(data = ratiodf, aes(x = group, y = deg0.05, group = 1), color = 'blue') +
    labs(x="Group", y="DEG ratio") +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

ggplot() +
    geom_line(data = ratiodf0.01, aes(x = group, y = Freq, group = 1), color = 'black') +
    geom_line(data = ratiodf, aes(x = group, y = deg0.01, group = 1), color = 'blue') +
    labs(x="Group", y="DEG ratio") +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'gray', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())


#####################################################################################################################
#TRUE DEG Selection 
#####################################################################################################################

#모든 gene에 대해서 deseq pvalue, DEG 여부, group 있음 

truedegl <- lapply(1:3, function(x){

    cut <- cutoff[x]
    degcut = degl[[x]]
    originaldeg = originaldegl[[x]][,c('pvalue', 'DEG', 'group')]
    candidatedeg = as.data.frame(subset(originaldeg, DEG == 1)[,c('DEG', 'group')])

    truedegg <- mclapply(1:num, mc.cores = 10, function(xx){

        bgline <- degcut[[xx]][,c('pvalue', 'DEG', 'group')]
        #bgline[is.na(bgline$pvalue)]$pvalue <- 1
        bglinedeg <- as.data.frame(subset(bgline, DEG == 1)[,c('DEG', 'group')])
        #originaldeg에서 bgline 빼야됨. 
        truedeg <- candidatedeg[!(rownames(candidatedeg) %in% rownames(bglinedeg)),]

        return(truedeg)

    })
    return(truedegg)
})
names(truedegl) <- cutoff

trueratio <- lapply(1:length(cutoff), function (x){

    cut <- cutoff[x]
    res <- truedegl[[xx]]

    ratio_cut <- mclapply(1:num, mc.cores = 10, function(xx){

        temp <- res[[xx]]

        ratio <- rep(0, 100)

        for (i in 1:nrow(temp)){

        j = temp[i,]$group

        ratio[j] = ratio[j] + 1 }

        ratio <- ratio / table(group)
        ratio <- as.data.frame(ratio)
        
        return(ratio)
    })

    return(ratio_cut)
    })

names(trueratio) <- cutoff

trueratiodf.0.1 = do.call(rbind, trueratio[[1]])
trueratiodf.0.05 = do.call(rbind, trueratio[[2]])
trueratiodf.0.01 = do.call(rbind, trueratio[[3]])

cutoff0.1 <- ggplot(data = trueratiodf.0.1, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'cornsilk') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="True DEG ratio") +
    #geom_hline(yintercept = 0.1, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

cutoff0.05 <- ggplot(data = trueratiodf.0.05, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'cornsilk') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="True DEG ratio") +
    #geom_hline(yintercept = 0.05, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

cutoff0.01 <- ggplot(data = trueratiodf.0.01, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'cornsilk') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="True DEG ratio") +
    #geom_hline(yintercept = 0.01, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(cutoff0.1,cutoff0.05, cutoff0.01, nrow=3, ncol=1)


