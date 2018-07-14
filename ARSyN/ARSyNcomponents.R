ARSyNcomponents<-function(asca=asca,Variability=0.75,beta=beta)
{
# This program selects the number of components that explain more than the Variability%
# If Variability="average" the number of components will be those that explain more than 
# the average variation of the principal components
# For residuals model the number of components selected are beta*average-variability.

MODEL<-asca[-length(asca)]
M<-length(MODEL)-1
output<-NULL

if(Variability=="average")
{
	library(Matrix)
	for (i in 1:M)
	{ 
	lim<-1/rankMatrix(MODEL[[i]]$X)[1]	
	t<-table(MODEL[[i]]$var.exp[,1]>lim)
	if(length(t)==1) {t[2]=0}
	t<-t[2]
	names(t)<-names(MODEL)[i]
	output<-c(output,t)
	}
}
if(Variability!="average")
{
	lim<-Variability
	for (i in 1:M)
	{ 
	t<-which(MODEL[[i]]$var.exp[,2]>lim)[1]
	names(t)<-names(MODEL)[i]
	output<-c(output,t)
	}	
}
### Residuals model
library(Matrix)
i=M+1 
lim <- beta*1/rankMatrix(MODEL[[i]]$X)[1]	
	t<-table(MODEL[[i]]$var.exp[,1]>lim)
	if(length(t)==1) {t[2]=0}
	t<-t[2]
	names(t)<-names(MODEL)[i]
	output<-c(output,t)

	
output
}
