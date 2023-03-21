rm(list=ls(all=TRUE)) #Clear memory
require(readr)
require(dplyr)
require(nnet)
require(tidyr)
require(ggplot2)
require(VGAM)
require(hnp)
require(pmultinom)
require(haven)
require(scatterplot3d)
require(gridExtra)


set.seed(1401)

ml <- read_dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")
glimpse(ml)

#Grouped Data
with(ml, table(math, prog)) 
dados<-read.csv("math.csv", head=TRUE, sep=";", dec=",") 
attach(dados)

############################################################################################################
#Histogram
Math<-rep(math,each=3)
Categ<-rep(1:3,times=length(math))
Freq<-as.vector(t(dados[,2:4]))
dados.plot<-cbind.data.frame(Math,Categ,Freq)

attach(dados.plot)

dados.plot$Categ<-as.factor(dados.plot$Categ)
levels(dados.plot$Categ)<-c("Academic","General", "Vocation")

#Histogram
ggplot(dados.plot, aes(x = Math, y = Freq))+geom_col(color="black",aes(fill = Categ),lwd=0.5)+
  theme(legend.position = "top",legend.title=element_blank())+labs(x= "Math Score", y = "Frequency")
############################################################################################################
#Plot 3d
colors <- c("red", "orange", "blue")
colors <- colors[as.numeric(dados.plot$Categ)]

scatterplot3d(dados.plot, type = "h", pch = " ", lwd = 3,
             y.ticklabs = c("Academic","","General","","Vacation"),
              xlab = "Math escore", ylab = "Category", zlab = "Frequency",
              y.margin.add = 0.2,
              color =colors, box = F,scale.y =2)

################################################################################################################
#models
mod0<- vglm(cbind(Academic,General,Vocation)~1,multinomial(refLevel=1),data=dados) 

mod1 <- vglm(cbind(Academic,General,Vocation)~math,multinomial(refLevel=1),data=dados)

#Ratio likelihood test
TRV <- 2*(logLik(mod1)-logLik(mod0))
gl <- length(coef(mod1))-length(coef(mod0))
p <- 1-pchisq(TRV,gl)
cbind(TRV, gl, p)

#AIC models
AIC(mod0)
AIC(mod1)

summary(mod1)

#Model Coefficients
coef(mod1, matrix = TRUE)
exp(coefficients(mod1)) #odds ratio

#Confidence Intervals for Parameters
confint(mod1)


#estimated probabilities
prob_pred<-fitted(mod1)


#observed probabilities - estimated probabilities
prob_final<-mod1@y-prob_pred

###########################################################################################
#hnp for proposed distances

n<-nrow(dados)

m<-vector()
for (i in 1:n) {
  m[i]<- sum(dados[i,2:4])
  
}


##Euclidean distance
d_eucl<- function(obj) {
  r <- resid(obj)
  l2_r <- apply(r, 1, function(x) dist(rbind(x,rep(0,length(x)))))
  return(as.numeric(l2_r))
}

#Mahalanobis distance
d_Mahal <- function(obj) {
  r <- resid(obj)
  k<-ncol(r)
  Sr<-solve(cov(r)+diag(rep(1e-6,k)))
  D2 <- mahalanobis(r, rep(0,nrow(r)),Sr,inverted = T)
  return(D2)
}

######################################################################################
#Implementing new class for hnp

s_fun <- function(n, obj) {
  pred<-fitted(obj)
  newresp<-matrix(NA,n,3)
  for (i in 1:n) {
    newresp[i,]<- t(rmultinom(1, size = m[i], pred[i,]))
    
  }
  newresp
}


########################################################################################################
#wrong model
f_fun.w <- function(newresp) {
  vglm(cbind(newresp[,1],newresp[,2],newresp[,3]) ~ 1, multinomial(refLevel=1),data=dados) 
}

#correct model
f_fun.r <- function(newresp) {
  vglm(cbind(newresp[,1],newresp[,2],newresp[,3])~ math, multinomial(refLevel=1),data=dados) 
}

#Plots
par(mfrow=c(1,2))

hnpE.w<-hnp(mod0, newclass = TRUE, diagfun = d_eucl, simfun = s_fun, fitfun = f_fun.w, print =T,sim = 1000,main="(a)",ylab="values of Euclidian Distance")


hnpE.r<-hnp(mod1, newclass = TRUE, diagfun = d_eucl,
            simfun = s_fun, fitfun = f_fun.r, print =T,sim = 1000,main="(b)",ylab="Values of Euclidian Distance")



hnpM.w<-hnp(mod0, newclass = TRUE, diagfun = d_Mahal,sim = 1000, simfun = s_fun, fitfun = f_fun.w,  print =T,main="(a)",ylab="values of Mahalanobis distance")


hnpM.r<-hnp(mod1, newclass = TRUE, diagfun = d_Mahal, simfun = s_fun, fitfun = f_fun.r, 
            print = T,sim = 1000,main="(b)",ylab="Values of Mahalanobis distance")

