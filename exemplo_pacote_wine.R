rm(list=ls(all=TRUE)) #Clear memory
# Loading the library
require(pmultinom)
require(hnp)
require(rattle.data)
require(dplyr)
require(nnet)
require(moments)
require(ggplot2)
require(readr)
require(tidyr)
###############################################################################################
#RQRs: PHI^-1[F(y-;prob_hat)+u*P(y;prob_hat)] with u~Uniform(0,1)

RQR.r <- function(yf,pred.r){
  res.quantile.r<-vector()
  for (i in 1:n) {
    if(yf[i]==1){
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(0,1,1),size = 1,probs=pred.r[i,],method="exact") +
                                   dmultinom(c(1,0,0),size = 1, pred.r[i,]) * runif(1))
      
    }else if(yf[i]==2){
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(1,0,1),size = 1,probs=pred.r[i,],method="exact") +
                                   dmultinom(c(0,1,0),size = 1, pred.r[i,]) * runif(1))
      
    }else{
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(1,1,0),size = 1,probs=pred.r[i,],method="exact")+
                                   dmultinom(c(0,0,1),size = 1, pred.r[i,]) * runif(1))} 
  }
  
  # res.quantile.rf<-(res.quantile.r-mean(res.quantile.r))/sd(res.quantile.r)
  return(res.quantile.r)
}
#######################################################################################################
set.seed(0417)
# Loading the wine data
data(wine)

# Checking the structure of wine dataset
str(wine)

# Using sample_frac to create a sample of wine values 
train <- sample_frac(wine, 1)
attach(train)

#Frequencies and plot
table(Type)

levels(train$Type)<-c("1", "2", "3")
ggplot(train, aes(x = Type,fill = Type)) +
  geom_bar(width=0.3,show.legend = FALSE) + ylab('Frequency')+xlab("Cultivar")

# Setting the basline 
train$Type <- relevel(train$Type, ref = "1")

# Training the multinomial model

multinom.fit0<- multinom(Type ~ 1, data = train)
multinom.fit1 <- multinom(Type ~ Phenols, data = train)
multinom.fit2 <- multinom(Type ~  Phenols+Magnesium, data = train)
multinom.fit3 <- multinom(Type ~ Phenols*Magnesium , data = train)


# Checking the model
anova(multinom.fit0,multinom.fit1,multinom.fit2,multinom.fit3)
summary(multinom.fit2)

AIC <- c(multinom.fit0$AIC,
         multinom.fit1$AIC,
         multinom.fit2$AIC,
         multinom.fit3$AIC)

confint(multinom.fit2,type = "Wald") 


#Plot obs probabilities vs estimated probabilities
pred.r<-predict(multinom.fit2, type = "prob")

categ_predict<-predict(multinom.fit2,se.fit=TRUE,interval = T) 
prob<-prop.table(table(categ_predict))

probs<-as.vector(t(prob)) #predict

p1<-as.vector(t(prop.table(table(Type)))) #observed


probfinal<-data.frame(cultivar=rep(1:3, times=2),
                      response=rep((1:3),times=2))


datafinal<-cbind(probfinal,proba=c(probs,p1),tipo=rep(c("Estimated","Observed"),each=3))

datafinal$cultivar<-as.factor(datafinal$cultivar)
datafinal$response<-as.factor(datafinal$response)

levels(datafinal$response)<-c("1","2","3")
levels(datafinal$cultivar)<-c("Cultivar 1","Cultivar 2", "Cultivar 3")

ggplot(datafinal, aes(x=response, y=proba, colour=tipo)) +
  geom_point(size=3) +
  xlim("1","2","3")+
  xlab('\n Cultivar \n')+ ylab('Proportion\n')+ylim(0,0.41)+
  scale_colour_manual(name="",breaks=c('Estimated','Observed'),
                      values=c('blue','red'))+ theme(legend.position="top")

#######################################################################################
#functions for residuals
yf<- as.numeric(train$Type)
n<-length(yf)

res.r<-RQR.r(yf,pred.r)

res.r.n <- (res.r-mean(res.r))/sd(res.r)

mean_x<-mean(res.r.n); sd_x<-sd(res.r.n)

hist(res.r.n, freq = FALSE, main = "",xlab = "Randomized quantile residuals",  col = c("lightblue"))
curve(dnorm(x, mean = mean_x, sd = sd_x), add = TRUE, col = "red",lwd = 2)

kurtosis(res.r.n);skewness(res.r.n)

#QQplot
qqnorm(res.r.n);abline(a=0,b=1)

#Shapiro-Wilk test
shapiro.test(res.r.n)$p.value
###################################################################################################
## hnp and Residual vs fitted values 

value_fit<-predict(multinom.fit2, type = "class")

ym<-matrix(NA,178,3)

for (i in 1:n) {
  if(value_fit[i]==1){
    ym[i,] <-c(1,0,0)
    
  }else if(value_fit[i]==2){
    ym[i,] <-c(0,1,0)
    
  }else{
    ym[i,] <-c(0,0,1)}}

x<-ym*pred.r
xf<-t(x[1:n,])[which(t(x[1:n,])>0)]


par(mfrow=c(1,2))

hnp(res.r.n,print=T,sim = 1000,ylab="Randomized quantile residual", main="(a)")

plot(xf,res.r.n,ylab = "Randomized quantile residual",xlab="Fitted values",main="(b)")
abline(h=0,col="red")
