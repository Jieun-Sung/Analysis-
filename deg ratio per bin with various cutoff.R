library(dplyr)
library(parallel)
library(DESeq2)
library(ggplot2)
library(Cairo)
library(ComplexHeatmap)

library(gtools)
library(gridExtra)

X11.options(type='cairo')

#directory load 
load("/spstorage/USERS/sung/DEGanalysis/OUTPUT/GSE111842/GSE111842_cutoff.Rdata")

#####################################################################################################################
#import data
#####################################################################################################################
output_dir = "/spstorage/USERS/sung/DEGanalysis/OUTPUT/GSE111842/expression-estimate"
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

genelength <- count.expr.list[[1]]$Length
count <-countdf0[,grep('X.spstorage.USERS.sung.DEGanalysis.OUTPUT.GSE111842.bam.files', colnames(countdf0))]
rownames(count) <- countdf0$SRR6835872.Geneid
rownames(count) <- gsub("\\.", "", rownames(count))
colnames(count) <- files

control_SRR <- coldata[coldata$strain == 'normal',]$files
count_control <- count[,which(colnames(count) %in% control_SRR)]
count0 <- filter(count, apply(count_control, 1, mean) > 4)
tpms <- apply(count, 2, function(x) raw_to_tpm(x, genelength))

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


#filtering
control_SRR <- as.character(coldata[coldata$strain == 'normal',]$files)
count_control <- count[,which(colnames(count) %in% control_SRR)]
count0 <- filter(count, apply(count_control, 1, mean) > 4)
count_control <- count0[,which(colnames(count0) %in% control_SRR)]

tpm_control <- tpms[which(rownames(tpms) %in% rownames(count_control)), which(colnames(tpms) %in% control_SRR)]

control.mean.tpm <- as.data.frame(apply(tpm_control, 1, mean))
colnames(control.mean.tpm) <- "meanTPM"
control.mean.tpm$Geneid <- rownames(count_control)
control.mean.tpm <- control.mean.tpm %>% arrange(Geneid) %>% arrange(meanTPM)

gene_order <- control.mean.tpm$Geneid
group <- rep(c(1:100), each = round(length(gene_order)/100))
group <- group[1:length(gene_order)]

#####################################################################################################################
#Control combination and DESEQ
#####################################################################################################################

#combination
#111842
comb <- as.data.frame(combinations(6, 3, as.character(control_SRR)))

#104310
comb <- as.data.frame(combinations(8, 4, as.character(control_SRR)))

#DESEq for each combination

num = nrow(comb)/2

deseqreslist <- mclapply(1:num, mc.cores = 10, function(xx){
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
    result <- results(ddsres)

    sizefactor <- estimateSizeFactors(dds)
    dispersion <- estimateDispersions(sizefactor)

    disp <- as.data.frame(dispersion@rowRanges@elementMetadata@listData$dispersion)
    rownames(disp) <- dispersion@rowRanges@partitioning@NAMES
    colnames(disp) <- c("dispersion")

    #check if the count files have the same gene ordering 
    gene_names <- sapply(disp, function(x) rownames(disp))
    as.data.frame(apply(gene_names, 2, function(x) all(x == rownames(dispersion))))

    totaldisp <- cbind(disp, as.data.frame(result))

    return(totaldisp)

    })


#####################################################################################################################
#ratio between controls 
#####################################################################################################################

ratio.0.1 <- rep(0, 100)
ratio.0.05 <- rep(0, 100)
ratio.0.01 <- rep(0, 100)

cutoff <- c(0.1, 0.05, 0.01)

ratio <- lapply(1:length(cutoff), function (x){

	cut <- cutoff[x]

	ratio_cut <- mclapply(1:num, mc.cores = 10, function(xx){

    res <- deseqreslist[[xx]]

    upDEG <- res[res$pvalue < cut & res$log2FoldChange > 0,]
    dnDEG <- res[res$pvalue < cut & res$log2FoldChange < 0,]

    res.order <- as.data.frame(res[gene_order, ])
    res.order <- cbind(res.order, group = group)

    upratio <- rep(0, 100)
    dnratio <- rep(0, 100)

    for (i in 1:nrow(res.order)){

    j = res.order[i,]$group

    if (rownames(res.order[i,]) %in% rownames(upDEG)){
            upratio[j] = upratio[j] + 1
         }

    else if (rownames(res.order[i,]) %in% rownames(dnDEG)){
            dnratio[j] = dnratio[j] + 1
         }
    }

    upratio <- upratio / table(group)
    upratio <- as.data.frame(upratio)
    #idx = upratio$Freq == 0
    #if (isTRUE(any(idx, 0)) == TRUE) {
    #    upratio[idx,]$Freq = 1
    #}
    #upratio$Freq = -log2(upratio$Freq)
    
    dnratio <- dnratio / table(group)
    dnratio <- as.data.frame(dnratio)
    #idx = dnratio$Freq == 0
    #if (isTRUE(any(idx, 0)) == TRUE) {
    #    dnratio[idx,]$Freq = 1
    #}
    #dnratio$Freq = log2(dnratio$Freq)
    
    upratio$type = as.factor("up")
    dnratio$type = as.factor("dn")

    ratio <- rbind(upratio, dnratio)
    
    return(ratio)

	})

	return(ratio_cut)

	})

names(ratio) <- cutoff

ratiodf.0.1 = do.call(rbind, ratio[[1]])
ratiodf.0.05 = do.call(rbind, ratio[[2]])
ratiodf.0.01 = do.call(rbind, ratio[[3]])

ratiodf.0.1_up = ratiodf.0.1[ratiodf.0.1$type == "up",]
ratiodf.0.1_dn = ratiodf.0.1[ratiodf.0.1$type == "dn",]
ratiodf.0.05_up = ratiodf.0.05[ratiodf.0.1$type == "up",]
ratiodf.0.05_dn = ratiodf.0.05[ratiodf.0.1$type == "dn",]
ratiodf.0.01_up = ratiodf.0.01[ratiodf.0.1$type == "up",]
ratiodf.0.01_dn = ratiodf.0.01[ratiodf.0.1$type == "dn",]

up0.1 <- ggplot(data = ratiodf.0.1_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.1 <- ggplot(data = ratiodf.0.1_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.1,dn0.1, nrow=2, ncol=1)

up0.05 <- ggplot(data = ratiodf.0.05_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.05 <- ggplot(data = ratiodf.0.05_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.05,dn0.05, nrow=2, ncol=1)

up0.01 <- ggplot(data = ratiodf.0.01_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.01 <- ggplot(data = ratiodf.0.01_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.01,dn0.01, nrow=2, ncol=1)


#####################################################################################################################
#ratio difference 
#"TRUE" DEG
#####################################################################################################################

#DESeq

cutoff <- c(0.1, 0.05, 0.01)

deseqratio <- mclapply(1:length(cutoff), mc.cores = 3, function(xx){

    cut = cutoff[xx]
    dds <- DESeqDataSetFromMatrix(countData = count0, colData = coldata, design = ~ strain)
    ddsres <- DESeq(dds)
    res <- results(ddsres)

    res$pvalue[is.na(res$pvalue)] = 1

    upDEG <- res[res$pvalue < cut & res$log2FoldChange > 0,]
    dnDEG <- res[res$pvalue < cut & res$log2FoldChange < 0,]

    res.order <- as.data.frame(res[gene_order, ])
    res.order <- cbind(res.order, group = group)

    upratio <- rep(0, 100)
    dnratio <- rep(0, 100)

    for (i in 1:nrow(res.order)){

    j = res.order[i,]$group

    if (rownames(res.order[i,]) %in% rownames(upDEG)){
            upratio[j] = upratio[j] + 1
         }

    else if (rownames(res.order[i,]) %in% rownames(dnDEG)){
            dnratio[j] = dnratio[j] + 1
         }
    }

    upratio <- upratio / table(group)
    upratio <- as.data.frame(upratio)
    #idx = upratio$Freq == 0
    #if (isTRUE(any(idx, 0)) == TRUE) {
    #    upratio[idx,]$Freq = 1
    #}
    #upratio$Freq = -log2(upratio$Freq)

    dnratio <- dnratio / table(group)
    dnratio <- as.data.frame(dnratio)
    #idx = dnratio$Freq == 0
    #if (isTRUE(any(idx, 0)) == TRUE) {
    #    dnratio[idx,]$Freq = 1
    #}
    #dnratio$Freq = log2(dnratio$Freq)

    upratio$type = as.factor("up")
    dnratio$type = as.factor("dn")

    deseqratio <- rbind(upratio, dnratio)

    return(deseqratio)

    })

names(deseqratio) <- cutoff

deseqratiodf.0.1 = as.data.frame(do.call(cbind, deseqratio[[1]]))
deseqratiodf.0.05 = as.data.frame(do.call(cbind, deseqratio[[2]]))
deseqratiodf.0.01 = as.data.frame(do.call(cbind, deseqratio[[3]]))

ratiodf.0.1 = do.call(rbind, ratio[[1]])
ratiodf.0.05 = do.call(rbind, ratio[[2]])
ratiodf.0.01 = do.call(rbind, ratio[[3]])

ratiodiff.0.1 <- as.data.frame(cbind(group = ratiodf.0.1_up$group, controlFreq = ratiodf.0.1_up$Freq, testFreq = deseqratio_up$Freq, 
    type = ratiodf.0.1_up$type))
attach(ratio0.1_up)
ratio0.1_up$degratio <- testFreq - controlFreq


up0.1 <- ggplot(data = ratiodf.0.1_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.1 <- ggplot(data = ratiodf.0.1_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.1, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.1,dn0.1, nrow=2, ncol=1)

up0.05 <- ggplot(data = ratiodf.0.05_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.05 <- ggplot(data = ratiodf.0.05_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.05, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.05,dn0.05, nrow=2, ncol=1)

up0.01 <- ggplot(data = ratiodf.0.01_up, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'bisque1') +
    #ggtitle("UP DEG ratio") +
    labs(x="Group", y="UP DEG ratio") +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

dn0.01 <- ggplot(data = ratiodf.0.01_dn, aes(x = group, y = Freq, group = group)) +
    geom_boxplot(fill = 'darkolivegreen2') +
    #ggtitle("DN DEG ratio") +
    labs(x="Group", y="DN DEG ratio") +
    geom_hline(yintercept = 0.01, linetype = 'solid', col = 'red', size = 1.0) +
    geom_vline(xintercept = c(16, 71, 97, 100), linetype = 'dotted', col = 'black', size = 1.0)  +
    theme(title = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.title.y = element_text(face = "bold", size = 20, color = "black")) +
    theme(axis.text.y = element_text(size = 15, color = "black")) +
    theme(axis.text.x = element_blank())

grid.arrange(up0.01,dn0.01, nrow=2, ncol=1)



