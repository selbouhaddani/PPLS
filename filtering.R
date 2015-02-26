filter_mar <- function(rna=rna){
#first, calculate the maximum of gene expression per each gene and take the quantiles, we are interested in top 25% of the gene expressions
maxGE <- apply(t(rna), 1, max)
sumGEmax <- summary(maxGE)
#next, take the IQR of each gene, to capture the variability and check the top 25% genes according to the size of IQR
IQRGE <- apply(t(rna), 1, IQR, na.rm=TRUE) 
sumGEIQR <- summary(IQRGE)
#selected genes/probes is the intersection of the two previous sets
filter2 <- (intersect(which(maxGE> sumGEmax[5]), which(IQRGE> sumGEIQR[5])))
return(filter2)
}