# Functions that determine if there is enrichment in the differentially regulated genes 
#
enrichTest <- function(DGELRT, sigmicro, targetScan,threshold=0.05,
context_threshold= 0) {
    # Returns a vector of Fisher's Exact generated p-values of how enriched the
    # significnatly expressed genes are in being miRNA targets.  
    require(edgeR)
    targetScan <- targetScan[targetScan$Total.context.score <= context_threshold,]
    allTargets <- table(targetScan[targetScan$Gene.Symbol %in% row.names(DGELRT$table), ]$miRNA.family) 
    sigGenes<-topTags(DGELRT, n=20000)$table[topTags(DGELRT, n=20000)$table$FDR
    <=threshold,]
    sigTargets <- table(targetScan[targetScan$Gene.Symbol %in% row.names(sigGenes),]$miRNA.family)
    # Make sure sigmicro is converted to sig
    all_hits <- c()
    sig_hits <- c()
    for(t in sigmicro){
	    tryCatch(z<-sigTargets[[t]], error=function(...){z<-0}) 
	    tryCatch(y<-allTargets[[t]], error=function(...){y<-0})
	    all_hits <- append(all_hits, y)
	    sig_hits <- append(sig_hits, z)
    }

    p_values <- c()

    temp <- matrix(c(0,0, 0,0), nrow=2,ncol=2)
    for( i in seq(1,length(sigmicro))){
        
	    temp[2,2] <- all_hits[i] - sig_hits[i] # NotDifferent/Target  
	    temp[1,2] <- sig_hits[i]# Different/Target
	    temp[1,1] <- length(rownames(sigGenes))-sig_hits[i] #Different/NotT
	    temp[2,1] <-
        length(rownames(DGELRT$table))-length(rownames(sigGenes))-all_hits[i]+sig_hits[i]#NotDifferent/NotTarget

	#Nondiff/Not T
	    tryCatch(p <-fisher.test(temp)$p.value, error=function(...){p<-NA}) 
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
	random_miRNAs <-sample(names(table(targetScan$miRNA.family)), length(sigmicro), replace=F)
	rand <- enrichTest(DGELRT, random_miRNAs, targetScan)
	data <- cbind(data, rand)	
    }
    data
}

beatRandommRNAs <- function(DGELRT, sigmicro, targetScan) {
    # Generates a table of p-values 

}

checkCoherence <- function(miRNA,DGELRT, targetScan, threshold=0.05, context_threshold=0) {
    # For each gene returns whether there is a negative coherance.  
    sigGenes<-topTags(DGELRT, n=20000)$table[topTags(DGELRT, n=20000)$table$FDR
    <=threshold,]
    
    targetScan <- targetScan[targetScan$Total.context.score <=context_threshold,]
    spec_target<- targetScan$miRNA.family==miRNA
    targetScan <- targetScan[spec_target,]
    overlap<-targetScan$Gene.Symbol %in% rownames(sigGenes)
    overlap <- as.vector(targetScan$Gene.Symbol[overlap])
    overlap <- sigGenes[overlap,]
    neg <- sum(overlap$logFC>=0)
    total <- dim(overlap)[1]
    ratio <-((total-neg)/total) 
    ratio 
    
    
}
check corr <- function(mirFam, DGELIST, targetScan, context_threshold=0) {
    targetScan <- targetScan[targetScan$Total.context.score <=context_threshold,]
    spec_target <- targetScan$miRNA.family==mirFam



}

varyingThresholds <- function(DGELRT, sigmicro, targetScan) {
    thresholds <- seq(0,-0.4,by=-0.02)
    p.values <-seq(1, length(sigmicro))
    for(i in thresholds){
	    t <- enrichTest(DGELRT, sigmicro, targetScan, context_threshold=i)
	    p.values <- cbind(p.values, t)
	}
    p.values <- p.values[,-1]
    t(p.values)
}

fix_context_score <- function(targetScan) {


}


convertMir <- function(DGELRT, mirFamilies, threshold = 0.05){
# Function to convert MiRBase IDs to seed sequences
    mirBase.id <- rownames(topTags(DGELRT, n=20000)$table[topTags(DGELRT, n=20000)$table$FDR
    <=threshold,])
    seeds <- mirFamilies[match(mirBase.id, mirFamilies$MiRBase.ID),"Seed.m8"]
    mirBase.id <-mirBase.id[!is.na(seeds)]
    seeds <- seeds[!is.na(seeds)]
    temp <- list()
    temp[["seeds"]] <- seeds
    temp[["mirBase.ID"]] <- mirBase.id
    return(temp)

}


plotPvalue <- function(dataframe, xval = c(0,0.1,0.2,0.3,0.4,0.5)){
    require(ggplot2)
    plot(xval, c(0,0.1,0.2,0.3,0.9,1), type="n") 
    for(i in seq(1,dim(dataframe)[1])){
        points(xval, -1*log10(dataframe[i,]),type="l", pch=19)
    }
}


