ARSyNmodel<-function(X = X,Factors=2,Designa = Designa,Designb = Designb,Designc = Designc,Variability="average",
Join =TRUE,Interaction=TRUE,beta=beta)
{

if(Factors==1)
{
  Fac0<-c(1,2)
  names(Fac0)<-c("Model.a","Model.res")
  asca0<- ASCA.1f(X=X, Designa=Designa, Fac=Fac0)
  Fac<-ARSyNcomponents(asca0,Variability=Variability,beta=beta)
  for (i in 1:length(Fac)){
	Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
	}
	asca<- ASCA.1f(X=X, Designa=Designa, Fac=Fac0)
}



if(Factors==2)
{
  Fac0<-c(1,2,2,2)
  names(Fac0)<-c("Model.a","Model.b","Model.ab","Model.res")
  if(Join){ names(Fac0)[3]<-c("Model.bab")}
  asca0<- ASCA.2f(X=X, Designa=Designa, Designb=Designb,Fac=Fac0,Join=Join,Interaction=Interaction)
  Fac<-ARSyNcomponents(asca0,Variability=Variability,beta=beta)
  for (i in 1:length(Fac)){
	Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
	}
	asca<- ASCA.2f(X=X, Designa=Designa, Designb=Designb,Fac=Fac0,Join=Join,Interaction=Interaction)
}

if(Factors==3)
{
  Fac0= c(0,2,2,2,2,2,2,2)
  names(Fac0)<-c("Model.a","Model.b","Model.c","Model.ab","Model.ac","Model.bc","Model.abc","Model.res")
  if(Join[1]){ names(Fac0)[4]<-c("Model.bab")}
  if(Join[2]){ names(Fac0)[5]<-c("Model.cac")}
  asca0<- ASCA.3f(X=X, Designa=Designa, Designb=Designb,Designc = Designc,Fac=Fac0,Join=Join,Interaction=Interaction)
  Fac<-ARSyNcomponents(asca0,Variability=Variability,beta=beta)
  for (i in 1:length(Fac)){
	Fac0[names(Fac[i])]<-Fac[names(Fac[i])]
	}
	asca<- ASCA.3f(X=X, Designa=Designa, Designb=Designb,Designc = Designc,Fac=Fac0,Join=Join,Interaction=Interaction)
}

Input<-asca$Input
asca$Input<-c(Input,Factors)
names(asca$Input)<-c(names(Input),"Factors")

output<-asca
output
}
