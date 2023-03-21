#Multinomial Simulation
rm(list=ls(all=TRUE)) #Clear memory
########################################################################
#Required packages
require(MASS) #To obtain Ordered Logistic or Probit Regression
require(flexmix) #To obtain Kullback-Leibler Divergence
require(stats) #To obtain Mahalanobis and Euclidean Distance 
require(nnet) #Fit Multinomial Log-linear Models
require(hnp) # To obtain Half-Normal Plot with Simulation Envelope
require(ggplot2)
require(dplyr)
require(tidyr)
require(gridExtra)

#############################################################################3

f <- function(m,x,xf,prob) {
  y <- t(apply(prob, 1, rmultinom,n=1, size = m))
  
  # Value of Y and X together
  dfM <- cbind.data.frame(y, x,xf)
  
  #Fit of model
  fit.w <- multinom(y ~ 1, data = dfM,trace=FALSE) #null model
  fit.r <- multinom(y ~ x+xf, data = dfM,trace=FALSE) #correct model
  
  p_value <- anova(fit.w, fit.r)[2,7] #p-value of the test
  lr_stat <- anova(fit.w, fit.r)[2,6]# Likelihood ratio stat
  ###########################################################################
  #hnp for proposed residuals by the distances
  ##Euclidean distance
  d_eucl<- function(obj) {
    r <- resid(obj)
    l2_r <- apply(r, 1, function(x) dist(rbind(x,rep(0,length(x)))))
    return(as.numeric(l2_r))
  }
  ######################################################################################
  #Mahalanobis distance
  d_Mahal <- function(obj) {
    r <- resid(obj)
    k<-ncol(r)
    Sr<-solve(cov(r)+diag(rep(1e-6,k)))
    D2 <- mahalanobis(r, rep(0,nrow(r)),Sr,inverted = T)
    return(D2)
  }
  #####################################################################################
  #Implementing new class for hnp
  n<-length(x)
  s_fun <- function(n, obj) {
    pred<-predict(obj, type = "prob")
    newresp<- t(apply(pred, 1, function(x) rmultinom(1, size = m, x)))
    newresp
  }
  
  f_fun.w <- function(newresp) {
    multinom(newresp ~ 1, data = dfM) 
  }
  
  f_fun.r <- function(newresp) {
    multinom(newresp ~ x+xf, data = dfM) 
  }
  ##################################################################################################
  #hnp for null model vs correct model
  #Euclidian Dist
  #Use of this function for suppressing convergence message
  invisible(capture.output(my_hnpE.w<-hnp(fit.w, newclass = TRUE, diagfun = d_eucl, simfun = s_fun, fitfun = f_fun.w, how.many.out = T,plot.sim = "FALSE")))
  invisible(capture.output(my_hnpE.r<-hnp(fit.r, newclass = TRUE, diagfun = d_eucl, simfun = s_fun, fitfun = f_fun.r, how.many.out = T,plot.sim = "FALSE")))
  
  #Percentage of points outside the envelope
  n_pointsE.w<-my_hnpE.w$out
  percE.w<-round((n_pointsE.w/my_hnpE.w$total)*100,2)
  
  n_pointsE.r<-my_hnpE.r$out
  percE.r<-round((n_pointsE.r/my_hnpE.r$total)*100,2)
  
  #Difference of points outside the envelope between the incorrect and correct models
  dif_E<-percE.w-percE.r
  
  #############################################################################################
  #Mahalanobis Dist
  invisible(capture.output(my_hnpM.w<-hnp(fit.w, newclass = TRUE, diagfun = d_Mahal, simfun = s_fun, fitfun = f_fun.w, how.many.out =  T,plot.sim = "FALSE")))
  
  invisible(capture.output(my_hnpM.r<-hnp(fit.r, newclass = TRUE, diagfun = d_Mahal, simfun = s_fun, fitfun = f_fun.r, how.many.out =  T,plot.sim = "FALSE")))
  
  #Percentage of points outside the envelope
  n_pointsM.w<-my_hnpM.w$out
  percM.w<-round((n_pointsM.w/my_hnpM.w$total)*100,2)
  
  n_pointsM.r<-my_hnpM.r$out
  percM.r<-round((n_pointsM.r/my_hnpM.r$total)*100,2)
  
  #Difference of points outside the envelope betweenthe incorrect and correct models
  dif_M<-percM.w-percM.r
  
  #############################################################################################
  #To show the results
  
  Incorrect_Cor_Model<-c("Eucl_Inc"=percE.w,"Eucl_Cor"=percE.r,"Mahal_Inc"=percM.w,"Mahal_Cor"=percM.r)
  Difer<-c("Dif_ED"=dif_E, "Dif_M"=dif_M)
  Value_test<-c("p-value"=p_value,"LR"=lr_stat )
  
  Result_final<-c(Incorrect_Cor_Model, Difer,Value_test)
  
  return(Result_final)
}

n<-50#number of samples
x<-rnorm(n) # continuous covariate
xf<-gl(2, n / 2, labels = c("Control", "Treatment")) # categorical covariate

#Assigning the values of intercepts and betas 
#The first level is assigned the role of reference.

z2<- 1.38-2.7*x+1.35*(xf=="Treatment")
z3<- 3.51-5.11*x+2.49*(xf=="Treatment")

den<-1+exp(z2)+exp(z3)
p1<-1/den; p2<-exp(z2)/den; p3<-exp(z3)/den
prob<-cbind(p1,p2,p3)


for(m in c(5,10,15)){
  cat("The denominator is equal",m,"\n")
  
  n_replic<-1000
  sim_scenario_1 <- replicate(n_replic,f(m,x,xf,prob))
  
  Dist<-rep(c("Euclidean Distance","Mahalanobis Distance"),each=2*n_replic)
  Model<-rep(c("Model null", "Model 2"),each=n_replic,times=2)
  Resp<-as.vector(t(sim_scenario_1[1:4,]))
  
  Final<-tibble::tibble(Dist,Model,Resp)
  
  ## Boxplots
  p1 <- ggplot(Final[1:(2*n_replic),], aes(x = Model, y = Resp)) + labs(y = "Number of points outside envelope (%)")+
    geom_boxplot(aes(fill=Model)) +
    facet_wrap(~Dist) +
    theme_bw()+
    theme(legend.position = "none",
          text = element_text(size=10),
          axis.text.x = element_text(angle=0))+ scale_y_continuous(limits=c(0,100))
  
  p1.1 <- ggplot(Final[(2*n_replic+1):(4*n_replic),], aes(x = Model, y = Resp)) + labs(y = "Number of points outside envelope (%)")+
    geom_boxplot(aes(fill=Model)) +
    facet_wrap(~Dist) +
    theme_bw()+
    theme(legend.position = "none",
          text = element_text(size=10),
          axis.text.x = element_text(angle=0))+ scale_y_continuous(limits=c(0,100))
  
  
  #Graph of differece of points 
  Dist2<-rep(c("Distância Euclidiana","Distância de Mahalanobis"),each=n_replic)
  seq<-rep(1:n_replic,times=2)
  Resp2<-as.vector(t(sim_scenario_1[5:6,]))
  
  Mean_Resp2<-round(apply(sim_scenario_1[-c((1:4),7,8),],1,mean),3)
  mean_differ<-rep(Mean_Resp2,each=n_replic)
  
  Final2<-tibble::tibble(Dist2,seq,Resp2, mean_differ)
  
  p2<-ggplot(Final2[1:n_replic,], aes(x=seq, y=Resp2)) + 
    geom_line(color="blue") +
    xlab("Iterações") +ylab("Diferença entre pontos fora do envelope (%)") +
    geom_hline(yintercept = 0)+
    geom_line(aes(y =mean_differ),color = "red",linetype="longdash") +
    facet_wrap(~Dist2)
  
  p2.1<-ggplot(Final2[(n_replic+1):(2*n_replic),], aes(x=seq, y=Resp2)) + 
    geom_line(color="blue") +
    xlab("Iterações") +ylab("Diferença entre pontos fora do envelope (%)") +
    geom_hline(yintercept = 0)+
    geom_line(aes(y =mean_differ),color = "red",linetype="longdash") +
    facet_wrap(~Dist2)
  
  grid.arrange(print(p1+labs(title = "(a)")),print(p1.1+labs(title = "(b)")),ncol=2)
  
  #To obtain the number of outliers
  out <- ggplot_build(p1)[["data"]][[1]][["outliers"]]
  out1 <- ggplot_build(p1.1)[["data"]][[1]][["outliers"]]
  
  
  number_outliersMD_Cor<-length(out1[[1]]);number_outliersMD_Inc<-length(out1[[2]])
  
  number_outliersED_Cor<-length(out[[1]]);number_outliersED_Inc<-length(out[[2]])
  
  
  n_outliers<-c(number_outliersED_Inc,number_outliersED_Cor,number_outliersMD_Inc,number_outliersMD_Cor)
  
  #mean and sd of points outside the envelope and median of values LR which p-values was less than 0.01
  #number of times that p-value was less than 0.01
  mean_value<-round(apply(sim_scenario_1[-(5:8),],1,mean),3)
  
  sd_value<-round(apply(sim_scenario_1[-(5:8),],1,sd),3)
  
  mean_LR<-round(mean(sim_scenario_1[8,][which(sim_scenario_1[7,]<0.01)]),3)
  
  n_pvalue<-length(sim_scenario_1[7,][which(sim_scenario_1[7,]<0.01)])
  
  #Final results
  Distances<-c("Incorrect Euclidian Distance","Correct Euclidian Distance","Incorrect Mahalanobis Distance","Correct Mahalanobis Distance")
  Results<-tibble(Distances,mean_value,sd_value,n_outliers)
  
  Values<-tibble(mean_LR,n_pvalue)
  
  
  print(knitr::kable(Results, caption = 'Values of incorrect and correct models considering Euclidean Distance and Mahalanobis Distance - Mean of percentage number of points outside, standard deviation (sd) of percentage number of points outside, and number of outliers.'))
  
  print(knitr::kable(Values, caption = 'Mean of likelihood ratio test (LR) and number of times that p-value < 0.01 between the incorrect and correct models.'))
  
}
