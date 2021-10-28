filtering <- function(genecount, filtering_num){
	print("zero count: ", nrow(genecount[rowSums(genecount)==0,]))
	genecount <- filter(genecount, apply(genecount, 1, mean) > filtering_num)
	genecount <- genecount[order(apply(genecount, 1, mean)),,drop = FALSE]

	return(genecount)
}