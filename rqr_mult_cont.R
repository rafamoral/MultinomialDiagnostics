rm(list=ls(all=TRUE)) 
require(nnet) #Fit Multinomial Log-linear Models
require(hnp) #half-normal plot using envelope simulation
require(pmultinom) #Calculate cdf of multinomial dist
require(ggplot2)
require(dplyr)
require(tidyr)
require(moments)
require(gridExtra)
#######################################################################################################
#RQRs: PHI^-1[F(y-;prob_hat)+u*P(y;prob_hat)] with u~Uniform(0,1)
RQR.r <- function(yf,y,pred.r){
  res.quantile.r<-vector()
  for (i in 1:n) {
    if(yf[i]==1){
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(0,1,1),size = 1,probs=pred.r[i,],method="exact")+
                                   dmultinom(y[i,],size = 1, pred.r[i,]) * runif(1))
      
    }else if(yf[i]==2){
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(1,0,1),size = 1,probs=pred.r[i,],method="exact")+
                                   dmultinom(y[i,],size = 1, pred.r[i,]) * runif(1))
      
    }else{
      res.quantile.r[i] <- qnorm(pmultinom(upper=c(1,1,0),size = 1,probs=pred.r[i,],method="exact")+
                                   dmultinom(y[i,],size = 1, pred.r[i,]) * runif(1))} 
  }
  return(res.quantile.r)
}

RQR.w <- function(yf,y,pred.w){
  res.quantile.w<-vector()
  
  for (i in 1:n) {
    if(yf[i]==1){
      res.quantile.w[i] <- qnorm(pmultinom(upper=c(0,1,1),size = 1,probs=pred.w[i,],method="exact") +
                                   dmultinom(y[i,],size = 1, pred.w[i,])* runif(1))
      
      
    }else if(yf[i]==2){
      res.quantile.w[i] <- qnorm(pmultinom(upper=c(1,0,1),size = 1,probs=pred.w[i,],method="exact") +
                                   dmultinom(y[i,],size = 1, pred.w[i,])* runif(1))
      
      
    }else{
      res.quantile.w[i] <- qnorm(pmultinom(upper=c(1,1,0),size = 1,probs=pred.w[i,],method="exact") +
                                   dmultinom(y[i,],size = 1, pred.w[i,])* runif(1))} 
  }
  return(res.quantile.w)
}


f <- function(x,prob) {
  y <- t(apply(prob, 1, rmultinom,n=1, size = 1))
  yf<-factor(apply(y, 1, function(x) which(x==1)))
  
  dfM <- data.frame(x,y)
  
  #Fit of model
  fit.w <- multinom(y~ 1,data = dfM,trace="FALSE") #null model
  fit.r <- multinom(y ~ x,data = dfM,trace="FALSE") #correct model
  
  #Predicted probabilities
  pred.w<-predict(fit.w, type = "prob")
  pred.r<-predict(fit.r, type = "prob")
  
  
  p_value <- anova(fit.w, fit.r)[2,7] #p-value of the test
  
  
  #Residuals
  res.mult.w.n<-RQR.w(yf,y,pred.w)
  res.mult.r.n<-RQR.r(yf,y,pred.r)
  
  invisible(capture.output(myhnp.w<-hnp(res.mult.w.n, print = T,scale=T,plot.sim = "FALSE",sim = 1000)))
  invisible(capture.output(myhnp.r<-hnp(res.mult.r.n, print=T,scale=T, plot.sim = "FALSE",sim = 1000)))
  
  #Percentage of points outside the envelope
  npoints.w<-myhnp.w$out
  perc.w<-round((npoints.w/myhnp.w$total)*100,2)
  
  npoints.r<-myhnp.r$out
  perc.r<-round((npoints.r/myhnp.r$total)*100,2)
  
  #Difference of points outside the envelope between the incorrect and correct models
  differ<-perc.w-perc.r
  
  
  ### Shapiro Wilk Test
  test.w<-shapiro.test(res.mult.w.n)$p.value
  test.r<-shapiro.test(res.mult.r.n)$p.value
  
  mean.w <- mean(res.mult.w.n)
  sd.w<- sd(res.mult.w.n)
  kurt.w<- kurtosis(res.mult.w.n)
  skew.w<-skewness(res.mult.w.n)
  
  mean.r <- mean(res.mult.r.n)
  sd.r<- sd(res.mult.r.n)
  kurt.r<- kurtosis(res.mult.r.n)
  skew.r<-skewness(res.mult.r.n)
  
  
  #To show the results
  
  Values_Model<-c("Perc.w"=perc.w,"Perc.r"=perc.r,"Differ"=differ)
  Value_test<-c("SW_test.w"=test.w,"SW_test.r"=test.r,"p-value"=p_value)
  Stats.w <-c("media.w"=mean.w,"sd.w"=sd.w,"kurt.w"=kurt.w,"skew.w"=skew.w)
  Stats.r<-c("media.r"=mean.r,"sd.r"=sd.r,"kurt.r"=kurt.r,"skew.r"=skew.r)
  
  Result_final<-c(Values_Model,Value_test,Stats.w,Stats.r)
  
  return(Result_final)
  
  
  
}

#Residual using multinomial
#3 categories and continuous covariate
set.seed(1401)
#Results of plots
n<-100#number of samples
x<-rnorm(n) #covariate
z2<-1.38-2.7*x
z3<- 3.51-5.11*x

den<-1+exp(z2)+exp(z3)
p1<-1/den; p2<-exp(z2)/den; p3<-exp(z3)/den
prob<-cbind(p1,p2,p3)

n_replic<-1000

sim_scenario <- replicate(n_replic,f(x,prob))



Res_m1<-rep("Randomized quantile residual",each=2*n_replic)
Model<-rep(c("Null model", "Model 1"),each=n_replic)
Resp<-as.vector(t(sim_scenario[1:2,]))
Resp_test<-round(as.vector(t(sim_scenario[4:5,])),5)

Final<-tibble::tibble(Res_m1,Model,Resp,Resp_test)


## Boxplots
p1 <- ggplot(Final, aes(x = Model, y = Resp)) + labs(x = "Model", y = "Points outside of envelope (%)")+
  geom_boxplot(aes(fill=Model)) +
  facet_wrap(~Res_m1) +
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)); p1

p2<- ggplot(Final[1:n_replic,], aes(x = Resp_test)) + labs(x = "P-value (Shapiro-Wilk test)", y = "Frequency")+
  geom_histogram(aes(fill=Model),color="darkblue", fill="lightblue")+
  facet_wrap(Model~.) +
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))

p2.1<- ggplot(Final[(n_replic+1):(2*n_replic),], aes(x = Resp_test)) + labs(x = "P-value (Shapiro-Wilk test)", y = "Frequency")+
  geom_histogram(aes(fill=Model),color="black", fill="pink")+
  facet_wrap(Model~.) +
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))


grid.arrange(print(p2+labs(title = "(a)")),print(p2.1+labs(title = "(b)")),ncol=2)


#To obtain the number of outliers
out <- ggplot_build(p1)[["data"]][[1]][["outliers"]]

n_outliers<-c(length(out[[2]]),length(out[[1]]))

#mean and sd of points outside the envelope and median of values LR which p-values was less than 0.01

mean_value<-round(apply(sim_scenario[-(3:14),],1,mean),3)

sd_value<-round(apply(sim_scenario[-(3:14),],1,sd),3)

mean_SW.w<-round(mean(sim_scenario[4,]),5)
mean_SW.r<-round(mean(sim_scenario[5,]),5)

number_SW.w<-length(sim_scenario[4,][which(sim_scenario[4,]<0.05)])
number_SW.r<-length(sim_scenario[5,][which(sim_scenario[5,]> 0.05)])


n_pvalue<-length(sim_scenario[6,][which(sim_scenario[6,]<0.01)])

mean_RQR.w<-round(mean(sim_scenario[7,]),5)
sd_RQR.w<-round(mean(sim_scenario[8,]),5)
kurt_RQR.w<-round(mean(sim_scenario[9,]),5)
skew_RQR.w<-round(mean(sim_scenario[10,]),5)

mean_RQR.r<-round(mean(sim_scenario[11,]),5)
sd_RQR.r<-round(mean(sim_scenario[12,]),5)
kurt_RQR.r<-round(mean(sim_scenario[13,]),5)
skew_RQR.r<-round(mean(sim_scenario[14,]),5)


#Final results


Residuals_mult<-c("Incorrect RQR Multinomial",  "Correct RQR Multinomial")
Results<-tibble(Residuals_mult,mean_value,sd_value,n_outliers)

Values<-tibble(n_pvalue,mean_SW.w,mean_SW.r,number_SW.w,number_SW.r)

Values_estat.w<-tibble(mean_RQR.w,sd_RQR.w,kurt_RQR.w,skew_RQR.w)

Values_estat.r<-tibble(mean_RQR.r,sd_RQR.r,kurt_RQR.r,skew_RQR.r)


print(knitr::kable(Results, caption = 'Values of incorrect and correct models considering RQR Multinomial - Mean and standard deviation (sd) of percentage number of points outside, and number of outliers.'))

print(knitr::kable(Values, caption = 'Number of times that p-value < 0.01, Shapiro Wilk test (SW) mean for models, number of times that SW < 0.05 for the incorrect model, and number of times that SW > 0.05 for the correct model.'))

print(knitr::kable(Values_estat.w, caption = 'Mean, sd., kurtosis and skewness of RQR wrong.'))

print(knitr::kable(Values_estat.r, caption = 'Mean, sd., kurtosis and skewness of RQR rigth.'))

###########################################################################################
#To plot
y <- t(apply(prob, 1, rmultinom,n=1, size = 1))
yf<-factor(apply(y, 1, function(x) which(x==1)))

dfM <- data.frame(x,y)

#Fit of model
fit.w <- multinom(y~ 1,data = dfM,trace="FALSE") #null model
fit.r <- multinom(y ~ x,data = dfM,trace="FALSE") #correct model

#Predicted probabilities
pred.w<-predict(fit.w, type = "prob")
pred.r<-predict(fit.r, type = "prob")

#Residuals
res.w<-RQR.w(yf,y,pred.w)
res.r<-RQR.r(yf,y,pred.r)
res.mult.w.n<- (res.w-mean(res.w))/sd(res.w)
res.mult.r.n<-(res.r-mean(res.r))/sd(res.r)

myhnp.w<-hnp(res.mult.w.n, print = T,sim = 1000,ylab="Randomized quantile residual")

myhnp.r<-hnp(res.mult.r.n, print= T,sim = 1000,ylab="Randomized quantile residual")

Final1<-Final[1:n_replic,]

Final2<-Final[(n_replic+1):(2*n_replic),]

par(mfrow=c(1,2))


hist(Final1$Resp_test,  xlab = "P-value (Shapiro-Wilk test)", ylab = "Frequency",col="red",main = "(a)")


hist(Final2$Resp_test, xlab = "P-value (Shapiro-Wilk test)", ylab = "Frequency",col="blue",main = "(b)")


