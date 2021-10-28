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

  enrRes = do.call(rbind, enrRes)
  enrRes$ID = names(refGMT)
  enrRes$adjpVal = p.adjust(enrRes$pVal, method = c("BH"))
  enrRes$logP = -log10(enrRes$adjpVal)
  enrRes = enrRes[,c('ID','adjpVal', 'pVal','logP','oddsRatio','tan','int','bg')]

  return(enrRes)
}