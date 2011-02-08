# Determines if there is enrichment in the differentially regulated genes 

enrichTest <- function(DGELRT, sigmicro, targetScan,threshold=0.05) {
    # Returns a vector of Fisher's Exact generated p-values of how enriched the
    # significnat genes are enriched.  
    require(edgeR)
    # Instead of passing DGEList AND sigGenes pass in an lrt.model
    allTargets <- table(targetScan[row.names(DGELRT$table) %in% targetScan$Gene.Symbol,]$miRNA.family) 
    #allTargets <- table(targetScan[targetScan$miRNA.family,]) 
    sigGenes<-topTags(DGELRT, n=20000)$table[topTags(DGELRT, n=20000)$table$FDR
    <=threshold,]
    sigTargets <- table(targetScan[targetScan$Gene.Symbol %in%
    row.names(sigGenes),]$miRNA.family)

    # Make sure sigmicro is converted to sig
    all_hits <- c()
    sig_hits <- c()
    for(t in sigmicro){
	tryCatch(z<-sigTargets[[t]], error=function(...){z<-0}) 
	tryCatch(y<-allTargets[[t]], error=function(...){y<-0})
	all_hits <- append(all_hits, y)
	sig_hits <- append(sig_hits, z)
    }

    # Do Fisher's Exact test (God I don't know if this is right
    p_values <- c()
    # Table of genes
    temp <- matrix(c(0,0, 0,0), nrow=2,ncol=2)
    for( i in seq(1,length(sigmicro))){
	temp[2,2] <- all_hits[i] - sig_hits[i] 
	temp[1,2] <- sig_hits[i]# Target different
	temp[1,1] <- length(rownames(sigGenes))-sig_hits[i] #Different/NotT
	temp[2,1] <-
	length(rownames(DGELRT$table))-length(rownames(sigGenes))-all_hits[i]+sig_hits[i]
	#Nondiff/Not T
	p <-fisher.test(temp)$p.value 
	p_values <- append(p_values,p )
    }
    p_values
}

beatRandommiRNAs <- function(DGELRT, sigmicro, targetScan,
threshold=0.05,sample=10){
    # Generates a table of p-values representing a random list of microRNAs
    # that have the same length as sigmicro
    data <- enrichTest(DGELRT, sigmicro, targetScan)

    for( i in seq(1,sample)){
	random_miRNAs <-sample(names(table(targetScan$miRNA.famil)), length(sigmicro), replace=F)
	rand <- enrichTest(DGELRT, random_miRNAs, targetScan)
	data <- cbind(data, rand)	
    }
    data
}

plotPvalue <- function(dataframe, xval = c(0,0.1,0.2,0.3,0.4,0.5)){
    plot(xval, c(0,0.1,0.2,0.3,0.9,1), type="n") 
    for(i in seq(1,dim(dataframe)[1])){
    points(xval, log10(dataframe[i,]),type="l", pch=19)
}
}
