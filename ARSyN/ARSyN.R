ARSyN <- function(data=data, Covariates=Covariates, Join=TRUE, Interaction=TRUE, Variability = 0.75,beta=2)
{
####################################
### --- Compute Inputs for ASCA
####################################

 X <- t(data) #conditions x genes
 Num.factors <- nrow(Covariates)
 labels.factors <- rownames(Covariates)

 Design <- list(NULL,NULL,NULL)
 for (i in 1:Num.factors)
 {
	x <- as.character(Covariates[i,])
	Design[[i]] <- make.ASCA.design(x)
 }

####################################
### --- Execute ASCAmodel 
####################################

 my.asca <- ARSyNmodel(Factors=Num.factors,X=X,Designa=Design[[1]],Designb=Design[[2]],Designc=Design[[3]],Join=Join,Interaction=Interaction,Variability=Variability,beta=beta)

#################################### 
### --- Writing filtered matrix 
####################################

 X.filtered <- X
M<-length(my.asca)-1

for (i in 1:(M-1))
{
	X.filtered <- X.filtered-my.asca[[i]]$E 
}

X.filtered <- X.filtered-my.asca[[M]]$TP

data.filtered <- t(X.filtered)
data.filtered
}