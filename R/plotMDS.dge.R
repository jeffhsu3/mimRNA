plotMDS.dge <- function (x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), ...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Last modified 13 August 2010
{
#   Check input
    require(ggplot2)
    if(is.matrix(x)) x <- DGEList(counts=x)
    if(!is(x,"DGEList")) stop("x must be a DGEList or a matrix")
    nprobes <- nrow(x)
    nsamples <- ncol(x)
    if(is.null(labels)) labels <- 1:nsamples
    if(ndim < 2) stop("dim.plot must be at least two")
    if(nsamples < ndim) stop("Too few samples")

    x$samples$group <- factor(rep.int(1,nsamples))
    twd <- estimateTagwiseDisp(estimateCommonDisp(x), prior.n = 10, grid.length = 500)
    o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:min(nprobes,top)]
    subdata <- x$counts[o,]
    dd <- matrix(0, nrow = nsamples, ncol = nsamples)
    for (i in 2:(nsamples)) 
    	for (j in 1:(i - 1))
            dd[i, j] = sqrt(estimateCommonDisp(DGEList(counts=subdata[,c(i,j)]))$common.dispersion)
    a1 <- cmdscale(as.dist(dd), k = ndim)
    label =c(0,-1.5,0,0,0,0,0,1.5)
    text_label <- c("Right_1", "Left_1", "Right_2", "Left_2", "Right_3",
    "Left_3","Right_4", "Left_4") 
    temp <- data.frame(Dimension2=a1[,dim.plot[2]],Dimension1=a1[,dim.plot[1]],sample=x$samples$patient,tissue=x$samples$tissue) 
    #g <- ggplot(data=temp, aes(x=Dimension1, y=Dimension2, shape=sample,
    #group=tissue ))
    #print(g + geom_vline(aes(yintercept=0,color="red"))+opts(title="MDS of Left-Right miRNA"))
    #plot(rep(1, nsamples), a1[, dim.plot[2]], type = "n", ylab = paste("Dimension", dim.plot[2]), ... )
    #text(1, a1[, dim.plot[2]], labels = labels, col=col, cex=cex)
    #return(invisible(dd))
    temp
     
}

plotMDS.GG <- function(dataframe){
    g<-ggplot(dataframe, aes(Dimension1, Dimension2, shape=sample,color=tissue))
    g + geom_point(size=7)+geom_vline(xintercept=0,col="red")+opts(title="miRNA counts MDS", plot.title=theme_text(size=48), axis.title.x=theme_text(size=30),
    axis.title.y=theme_text(angle=90,size=30), axis.text.y=theme_text(size=20),
    axis.text.x=theme_text(size=20))
}
