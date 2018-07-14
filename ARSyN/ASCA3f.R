ASCA.3f<-function(X = X,Designa = Designa,Designb = Designb,Designc = Designc,Fac = c(1,2,2,2,2,2,2,2),
Join =c(TRUE,TRUE),Interaction=c(TRUE,TRUE,TRUE,TRUE))

{
#--------------------------------------------------------------------------------------
#  Dimensions of the matrices: 
#  X (p x n) contains expression values of n genes (in columns) and p conditions (in rows)
#  Designa (p x I) contains 0's and 1's for the TIME-POINTS in the experiment
#  Designb (p x J) EXPERIMENTAL GROUP FACTOR 1  
#  Designc (p x K) ANOTHER FACTOR  

#  Join = c(TRUE,TRUE) if the analyses of the model b and ab and c and ac is studied jointly
#  Interaction = c(TRUE,TRUE,TRUE,TRUE) to consider interaction "ab", "ac", "bc" and "abc" in the separated model

n<-ncol(X)
p<-nrow(X)
I<-ncol(Designa)
J<-ncol(Designb)
K<-ncol(Designc)

Faca=Fac[1]# number components Model a (time)
Facb=Fac[2] # number components Model b  (second factor)
Facc=Fac[3] # number components Model c  (third factor)
Facab=Fac[4] # number components Model ab (interaction) 
Facac=Fac[5] # number components Model ac (interaction) 
Facbc=Fac[6] # number components Model bc (interaction) 
Facabc=Fac[7]  # number components Model abc (interaction) 
Facres=Fac[8]
#----------------------- Calculate Overall Mean -------------------------------------

offset<-apply(X,2,mean)
Xoff<-X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))

#-----------------------  PART I: Submodel a (TIME) ---------------------------------

Model.a<-ASCAfun1(Xoff,Designa,Faca)
Xres<-Xoff-Model.a$X

#-------------------------- PART II.1: Submodel b.ab-----------------------------------
if(!Join[1]) 
{
	Model.b<-ASCAfun1(Xoff,Designb,Facb)
	if (Interaction[1]) {
		Model.ab<-ASCAfun2(Xoff,Designa,Designb,Facab)
		}
}

if(Join[1])
{
	Model.bab<-ASCAfun12(Xoff,Designa,Designb,Facab)
}
#-------------------------- PART II.2: Submodel (c.ac) -------------------------------
if(!Join[2]) 
{
	Model.c<-ASCAfun1(Xoff,Designc,Facc)
	if (Interaction[2]) {
		Model.ac<-ASCAfun2(Xoff,Designa,Designc,Facac)
	}
}

if(Join[2])
{
	Model.cac<-ASCAfun12(Xoff,Designa,Designc,Facac)
}
#-------------------------- PART II.3: Submodel (bc) --------------------------------

	if (Interaction[3]) {
		Model.bc<-ASCAfun2(Xoff,Designb,Designc,Facbc)
	}
#-------------------------- PART II.4: Submodel (abc) --------------------------------

	if (Interaction[4]) {
		Model.abc<-ASCAfun.triple(Xoff,Designa,Designb,Designc,Facabc)
	}
# ------------------------Collecting models ------------------------------------------

models <- ls(pattern="Model")
output <- vector(mode="list")
Xres <- Xoff
for (i in 1: length(models)) {
	mymodel <- get(models[i], envir=environment())
	output <- c(output, list(mymodel))
	Xres <- Xres - mymodel$X
	rm(mymodel)
	gc()
}
names(output) <- models

#------------------------- PART III: Submodel res -----------------------------------

Model.res<-ASCAfun.res(Xres,Facres)

Model.res<-list(Model.res)
names(Model.res)<-c("Model.res")
output<-c(output,Model.res)



#------------------------- Add Input data to the Output ----------------------------

Input<-list(X, Designa, Designb, Designc, Fac, Join,Interaction,Xoff)
names(Input)<-c("X", "Designa", "Designb", "Designc", "Fac", "Join","Interaction","Xoff")
Input<-list(Input)
names(Input)<-"Input"
output<-c(output,Input)

output 

}

