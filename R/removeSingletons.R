################################################################
#
# R script to remove singleton genes that are expressed in only 
# 1 sample.  These genes really mess up dispersion calculations.
#
################################################################

removeSingletons <- function(DGEList, threshold=10, top="None"){
    if(!("DGEList" %in% class(DGEList)))
	stop("Input should be of class DGEList.")

    nonsingletons <- apply(DGEList$counts, 1, isHighest)
    newDGEList <- DGEList
    newDGEList$counts <- newDGEList$counts[nonsingletons,]
    return(newDGEList) 

}

isHighest <- function(x, threshold=10, top="None"){
    # Top should be a percentage ie 0.6
    highest <-which(x==max(x))
    if( sum(x[-highest]) > threshold){
	return(TRUE)
    }else{
	return(FALSE)
    }

}	
