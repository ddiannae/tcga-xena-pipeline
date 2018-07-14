make.ASCA.design <- function(x)
{
x<-as.factor(x)
levels<-unique(x)
n<-length(x)
p<-length(levels)

Design<-matrix(0,nrow=n,ncol=p)

for (i in 1:n){
	for (j in 1:p){
	
	if (x[i]==levels[j])
	{
	Design[i,j]=1
	}
	}
}
colnames(Design)<-levels

output<-Design
output
}
