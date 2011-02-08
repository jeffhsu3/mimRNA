##########################################################################
#
# Plots paired samples with line connecting paired samples.  GLMDGE is a class
# from a glmFit from EdgeR
#
##########################################################################

plotPairs <- function(gene, DGEGLM){
    if(!("DGEGLM" %in% class(DGEGLM)))
	stop("Input should have class \"DGEGLM\".")
    
    num_contrasts <- dim(DGEGLM$design)[2]-1
    classes <- DGEGLM$design[,num_contrasts+1]
    plot(classes,DGEGLM$fitted.values[gene,],  pch=19,
	ylab="Adjusted Counts",xlab="Classes", xaxt="n", xlim=c(-0.5,1.5),
	main=gene)
    # :TODO Need to change the labels
    axis(1, at=c(0,1), labels=c("Left Atria", "Right Atria"))
    lastone <- rep(0, dim(DGEGLM$design)[1])
    # Adding the lines to connect pairs
    
    for(x in seq(1,num_contrasts-1)){
	lastone <- DGEGLM$design[,x+1]+lastone
	index <- ifelse(DGEGLM$design[,x+1]==1,T,F)
	points(classes[index],DGEGLM$fitted.values[gene,index], col=x,
	type="l")
    }
    lastone <- ifelse(lastone==1, F,T)
    points(classes[lastone], DGEGLM$fitted.values[gene,lastone],
    col=num_contrasts, type="l")
    leg.txt <- colnames(DGEGLM$design)[2:num_contrasts]
    leg.txt <- append(leg.txt,"patientDonor_182")
    legend("topright",legend=leg.txt, fill=seq(1,length(leg.txt),by=1))
}

##########################################################################
#
# Plot residual
#
##########################################################################


plotPairsOrig <- function(gene, DGEList, design){
    require(ggplot2)
    if(!("DGEList" %in% class(DGEList)))
	stop("Input should have class \"DGEList\".")
    # Need to change how the normalization works
    
    num_contrasts <- dim(design)[2]-1
    classes <- design[,num_contrasts+1]
    effective.library <- (DGEList$samples$lib.size*DGEList$samples$norm.factors)
    rpm <- DGEList$counts[gene,]/effective.library*1e6 
    Classes <- factor(ifelse(classes==0, "Left Atria", "Right Atria"))
    print(Classes)
    patient <- c(0,0,1,1,2,2,3,3)
    temp <- cbind(rpm, patient)
    temp <- as.data.frame(temp)
    temp <- cbind(temp, Classes)
    print(temp)
    plot(classes,rpm,  pch=19,
	ylab="Reads per Million",xlab="Classes", xaxt="n", xlim=c(-0.5,1.5),
	main=gene)
    axis(1, at=c(0,1), labels=c("Left Atria", "Right Atria"))
    lastone <- rep(0, dim(design)[1])

    # Adding the lines to connect pairs
    for(x in seq(1,num_contrasts-1)){
	lastone <- design[,x+1]+lastone
	index <- ifelse(design[,x+1]==1,T,F)
	points(classes[index],DGEList$counts[gene,index]/effective.library[index]*1e6, col=x,
	type="l")
    }
    lastone <- ifelse(lastone==1, F,T)
    points(classes[lastone],
    DGEList$counts[gene,lastone]/effective.library[lastone]*1e6
    ,col=num_contrasts, type="l")
    leg.txt <- colnames(design)[2:num_contrasts]
    leg.txt <- append(leg.txt,"patientDonor_182")
    legend("topright",legend=leg.txt, fill=seq(1,length(leg.txt),by=1))
    p<-ggplot(temp, aes(Classes, rpm, group=patient))
    p +
    geom_line()+geom_point()+scale_x_discrete("Classes")+scale_y_continuous("Reads
    per Million")+opts(title=paste("Left-Right Differences for", gene))
}

plotPairsGG <- function(gene, DGEList, design){
    require(ggplot2)
    if(!("DGEList" %in% class(DGEList)))
	stop("Input should have class \"DGEList\".")
    # Need to change how the normalization works
    
    num_contrasts <- dim(design)[2]-1
    classes <- design[,num_contrasts+1]
    effective.library <- (DGEList$samples$lib.size*DGEList$samples$norm.factors)
    rpm <- DGEList$counts[gene,]/effective.library*1e6 
    Classes <- factor(ifelse(classes==0, "Left Atria", "Right Atria"))
    patient <- c(0,0,1,1,2,2,3,3)
    temp <- cbind(rpm, patient)
    temp <- as.data.frame(temp)
    temp <- cbind(temp, Classes)
    p<-ggplot(temp, aes(Classes, rpm, group=patient))
    p +
    geom_line()+geom_point(size=4)+scale_x_discrete("")+scale_y_continuous("Reads
    per Million")+opts(title=paste("LR Differences in",
    gene),plot.title=theme_text(size=48), axis.title.x=theme_text(size=30),
    axis.title.y=theme_text(angle=90,size=30), axis.text.y=theme_text(size=20),
    axis.text.x=theme_text(size=20))
}
##########################################################################
#
# Remove single expressors.  If one gene is expressed in only 1 sample and zero
# in the others, this makes it so that under the negative binomial it appears
# that the gene is differentially expressed between classes (since
# dispersion is shrunk towards a common value
#
##########################################################################

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
