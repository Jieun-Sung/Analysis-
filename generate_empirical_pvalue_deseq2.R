###################################################################
###  Generate Empirical Pvalue - DESEQ2
###################################################################

library(yaml)
library(data.table)
library(tidyverse)
library(parallel)

setDTthreads(10)

args <- commandArgs(trailingOnly=TRUE)
cfg <- read_yaml(args[1])
ncore <- args[2]
print(cfg)

# Input directory 
count_dir <- cfg$count_dir
#count_dir <- '/spstorage/DATA/kmap_shares/KMAP-2K-EXPRESSION-PROFILES/HDFn/Batch-Effect-Removed-HTseq-Count-Matrix'

pinfo_dir <- cfg$pinfo_dir
#pinfo_dir <- "/spstorage/DATA/kmap_process/deseq_diffres/sample_HTseq_Metainfo.csv"

# Output directory
fig_dir <- cfg$fig_dir
out_dir <- cfg$out_dir
# fig_dir <- '/spstorage/USERS/sung/projects/KMAP_Signature/fig'
# out_dir <- '/spstorage/USERS/sung/projects/KMAP_Signature/output'

if(!dir.exists(fig_dir)) dir.create(fig_dir)
if(!dir.exists(out_dir)) dir.create(out_dir)

source('/spstorage/USERS/sung/projects/KMAP_Signature/script/generate_empirical_pvalue_functions.R')

##====================================================================================
# ncore <- cfg$ncore
nbin = 100                  # number of bins 
w = 4  						#w :window size  	
sampling = 3 				#s : sampling size
iter = 50 					#iter : iteration number 	

##====================================================================================


##============ Load data ============##
pinfo <- read.csv(pinfo_dir, header = TRUE)
pinfo$sampleid <- gsub('-','.', pinfo$sampleid)
platenum <- paste0(rep(unique(pinfo$plate), 3), '_', rep(unique(pinfo$dose), each = 15))

print(platenum)

rawcount_kmap_file = list.files(count_dir)	
countmatrix_list <- mclapply(rawcount_kmap_file, mc.cores = ncore, function(file){
    path = paste0(count_dir,"/",file)
    countmatrix = data.frame(fread(path, header = TRUE, stringsAsFactors = FALSE))
    rownames(countmatrix) <- countmatrix[,1]; countmatrix <- countmatrix[,-1]
    return (countmatrix)}
)
names(countmatrix_list) <- platenum

DMSO_count_filtered_list <- mclapply(countmatrix_list, mc.cores = ncore, function(countmatrix){
	DMSO_count <- countmatrix[,grep('DMSO', colnames(countmatrix))]
	#nrow(DMSO_count[rowSums(DMSO_count)==0,])
	DMSO_count_filtered <- filter(DMSO_count, apply(DMSO_count, 1, function(row) all(row!=0)))
	DMSO_count_filtered <- DMSO_count_filtered[order(apply(DMSO_count_filtered, 1, mean)),,drop = FALSE]
	#nrow(DMSO_count_filtered)
	
	return(DMSO_count_filtered)
})
names(DMSO_count_filtered_list) <- platenum

drug_count_list <- mclapply(countmatrix_list, mc.cores = ncore, function(countmatrix) { 
	temp <- strsplit(colnames(countmatrix), '[.]')
	temp <- unique(unlist(lapply(temp, function(xx) xx[4])))
	temp <- temp[!grepl('DMSO', temp)] ; temp <- temp[!grepl('Niclo', temp)]

	drug_count <- lapply(temp, function(name){
		countidx <- grep(name, colnames(countmatrix), value = TRUE)

		return(countmatrix[, countidx])
	})
	names(drug_count) <- temp
	return(drug_count)
})

deseqres_list <- mclapply(1:length(countmatrix_list), mc.cores = core, function(xx){
	if (xx %% 100 == 0) {print(paste0(xx, "th drug processing"))}

	drug_count <- drug_count_list[[xx]]
	version = strsplit(names(drug_count_list)[[xx]], '[.]')[[1]][1]
	dose = strsplit(names(drug_count_list)[[xx]], '[.]')[[1]][2]

	deseqres <- deseq_drug(drug_count, version, dose)
	return(deseqres)
	})
names(deseqres_list) <- names(drug_count_list)


##============ Calculation ============##

deseq_pval_list <- lapply(1:length(DMSO_count_filtered_list), function(xx){
	tcheck <- proc.time()
	print(paste0(xx, "th plate processing"))

	DMSO_count <- DMSO_count_filtered_list[[xx]]
	pval_list <- deseq_DMSO(DMSO_count, w)
	names(pval_list) <- rownames(DMSO_count)
	
	print((proc.time() - tcheck)/60)

	return(pval_list)
})
names(deseq_pval_list) <- platenum

emppval_from_deseq_perbin_list <- DMSO_pval_to_empiricalpval(deseq_pval_list)
names(emppval_from_deseq_perbin_list) <- platenum

loess_deseq_list <- loess_empirical(emppval_from_deseq_perbin_list)
names(loess_deseq_list) <- platenum

deseqres_list <- mclapply(1:length(drug_count_list), mc.cores = core, function(xx){
	if (xx %% 100 == 0) {print(paste0(xx, "th drug processing"))}

	drug_count <- drug_count_list[[xx]]
	version = strsplit(names(drug_count_list)[[xx]], '[.]')[[1]][1]
	dose = strsplit(names(drug_count_list)[[xx]], '[.]')[[1]][2]

	deseqres <- deseq_drug(drug_count, version, dose)
	return(deseqres)
	})
names(deseqres_list) <- names(drug_count_list)

deseq_group_quantile_list <- group_quantile(deseq_pval_list)
names(deseq_group_quantile_list) <- platenum

emp_deseq_res_list <- predict_empiricalPvalue(deseqres_list, loess_deseq_list, deseq_pval_list, deseq_group_quantile_list)
names(emp_deseq_res_list) <- names(deseqres_list)
