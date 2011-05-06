# A class to compare matrices and plot them

setClass("MAM", representation(mat1 = "matrix", mat2 = "dataframe"))

setMethod("initialize", "MAM", function(.Object, mat1, mat2, annot=1){
    # input is a dataframe with the same dimension and some sort of identifier
    # column
    # Check if mat are the same dimension
    if(dim(mat1) == dim(mat2)){
    # Need to check the NA case
        match1 <- match(mat1[,annot], mat2[,annot])
        match2 <- match(mat1[,annot], mat2[,annot])
        if all(match1==match2)==FALSE {
            mat1 <- mat1[match2,]
        }
        .Object@mat1 <- mat1[,-annot] 
        .Object@mat2 <- mat2[,-annot] 
        .Object@annot <- mat1[,annot]
        .Object
    }
})

setMethod("show", "MAM", 
function(object){
    cat("Matrix and Matrix Object \n")
    cat("Matrix1:", paste(head(object@mat1)), "\n")
    cat("Matrix2:", paste(head(object@mat2)), "\n")
})

setMethod("plot", "MAM", 
function(object, design, treatment=1, ...){
    # The design is a model matrix and treatment determines which particular
    # treatment in which to lay out the 2-D array of correlation plots. 
    
    # A basic plot to plot the correlation
    require(ggplot2)
    # :TODO fix this so that it incorporates expreimental design

    features <- dim(object@mat1)[1]
    num_samples <- dim(object@mat1)[2]
    samples <- rep(seq(1, num_samples), each=features)
    melt1 <- melt(object@mat1)
    melt2 <- melt(object@mat2)
    mc <- dataframe(samples = sample ,mat1 = melt1, mat2 = melt2)
    
    g <- ggplot(mc, aes(melt1, melt2)) + facet_grid(sample~.) +
    stat_abline(col="red") + geom_point()


})

setMethod("permuteRows", "MAM",
function(object, permute=1000,){

})

setMethod("


mutualMatch <- function(data1,data2, identifier=1) {


}
