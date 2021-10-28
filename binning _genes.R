# binning genenames 
group <- rep(c(1:nbin), each = round(length(geneset)/nbin))
group <- group[1:length(geneset)]
names(group) <- names(geneset)
group_quantile <- data.frame(group = group, percentile = seq(from = 1, to = nbin, length = length(geneset)))

grouping_by_bin <- function(geneset, nbin){
	#geneset = dataframe 
	#if list, use len instead of nrow 
	name_per_bin <- lapply(1:nbin, function(xx){
		name_per_bin_temp = {}
		for (i in 1:length(group)){
			if (group[i] == xx) {
				name_per_bin_temp = append(name_per_bin_temp, names(group[i]))
			}
		}
		return(name_per_bin_temp)
	})
	return (name_per_bin)
}
