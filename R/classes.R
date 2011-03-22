setClass("MIM", representation(mRNA = "matrix", miRNA = "matrix")) 

setMethod("initialize", "MIM", function(.Object, mRNA, miRNA){
    .Object@mRNA <- mRNA
    .Object@miRNA <- miRNA
    .Object
})


meetThreshold <- function(correlation,miRNA, threshold=0.95){
    mt <- abs(correlation[miRNA,]) >= threshold 
    mt
}

compareChambers <- function(miRNA, left, right){
    left <- meetThreshold(left, miRNA)
    right <- meetThreshold(right, miRNA)
    inBoth <- left*right
    print(sum(inBoth, na.rm=T))
    inBoth <- ifelse(inBoth == 0, FALSE, TRUE)
}

permutateSamples <- function(mRNA, miRNA, pmut=24){
    corr <- cor(t(miRNA), t(mRNA))
    # Initialze empty matrix
    t <- mat.or.vec(dim(corr)[1], dim(corr)[2])
    colnames(t) <- rownames(mRNA)
    rownames(t) <- rownames(miRNA)
    # Hmm right now not accounting for the paired nature of the samples
    # Need to permute only the left or the right samples
    for (x in c(1:pmut)){
	    index <-seq(1,dim(miRNA)[2]) 
	    permute <- sample(index)
	    if (all(permute == index)){
	    }else{
	        permute_miRNA <- miRNA[,permute] 
	        cor_permute <- cor(t(permute_miRNA), t(mRNA))
	        temp <- (abs(corr) <= abs(cor_permute))*1
	        t<-t + temp
	        print(paste("Permutation ", x))
	}
    }
    # Returns a truth matrix
    # Must beat 1/24 to be p <= 0.05
    t <- (t <= 1)
    t
}

comparemiRNAmRNA <- function(correlation_matrix, mapping){
    require(ggplot2)
    correlations <- c()
    for (x in seq(1, dim(mapping)[1])){
	t <- 0
	tryCatch(t<-correlation_matrix[mapping[x,1], mapping[x,2]],
	error=function(...){t<<-NA})
	print(t)
	correlations <- append(correlations, t)
    }
    
    correlations
}


plotBackgroundThreshold <- function(expMat, threshold = seq(10,400,
by=10)){
    require(ggplot2)
    meetThresh <- c()
    for (x in threshold){
	meetThresh <- append(meetThresh, sum(rowSums(expMat)> x))
    }
    test <- data.frame(threshold = threshold, genes = meetThresh)
    plot(test$threshold, test$genes, type="l") 

}

removeDuplicates <- function(expressionMatrix){
# Removes rows that have the same counts for a transcripts as the counting
# method cannot resolve differences in the isoforms/overlapping transcripts. 
    t <- rep(1, times= dim(expressionMatrix)[1])
    for (x in seq(1,dim(expressionMatrix)[2])){
	temp <- duplicated(expressionMatrix[,x])
	t <- temp * t
    }
    t <- ifelse(t==0, FALSE, TRUE)
    expressionMatrix <- expressionMatrix[!t,]
    expressionMatrix
    
}

getSignficianatBoth <- function(miRNA,cor.left, cor.right, threshold=0.96,
sign=T){
    # Get the significant genes that are correlated in both left and right
    # :TODO Need to make sure the directionality is the same.
    sig.left <- which(abs(cor.left[miRNA,])>= threshold)
    sig.right <- which(abs(cor.right[miRNA,]) >= threshold)
    #sig.left<-colnames(cor.left)[sig.left]
    #sig.right<-colnames(cor.right)[sig.right]
    matches <- sig.right[sig.right %in% sig.left]
    matches
}


targetEnrich <- function(miRNA,targetscan_db, both=TRUE){
if (both == TRUE){
    predicted <- which(miRNA == targetscan_db$Representative.miRNA)
    predicted <- targetscan_db[predicted, "Gene.Symbol"]
    # These are global cor.left, cor.right are global variables
    print(head(predicted))
    empirical <- getSignficianatBoth(tolower(miRNA), cor.left, cor.right, threshold =
    0.96)
    matches <- empirical %in% predicted
    matches
} else 	{
:


}


}


