#Simulating Lognormal data and fitting three models including the Bayesian quantile regression (BQR), 
# and Bayesian quantile regressions of two-part models with hurdles at 0 and 4, then computing the mean square error, prediction error, and 
#The width of credible intervals corresponding to the estimates based on the three models


library(R2jags)

sim <- function(quant = tau, n = n, M =10){

  #Providing free space for assining the results  
  widcred <- matrix(NA, nrow = M, ncol = 4)
  widcred0 <- matrix(NA, nrow = M, ncol = 4)
  widcred4 <- matrix(NA, nrow = M, ncol = 4)
  MSE<-  matrix(NA, nrow = M, ncol = 1) 
  MSE0<-  matrix(NA, nrow = M, ncol = 1) 
  MSE4<-  matrix(NA, nrow = M, ncol = 1) 
  prederrors <- matrix(NA, nrow = M, ncol = n)
  prederrors0 <- list()
  prederrors4 <- list()  
  
  for(i in 1:M){


#Simulating data from a lognormal distribution      
tau<-quant
epsilon<-rnorm(n)
x1<-rlnorm(n,2,2)
x2<-rnorm(n,0.5,.5)
sigma=0.4
mu =2-0.2*x1+0*x2 + epsilon
y=rlnorm(n,mu,sigma)
f<-function(t) ifelse(t<4, floor(t),t)
data<-data.frame(cbind(y=f(y),x1=x1,x2=x2))

#Data withoout zeros
data0<-data[which(data$y>0),]

#Data withoout zero, one, two, three      
data4<-data[which(data$y>=4),]



#Whole data with a mass point at values<4
f_jit<-function(t) ifelse(t<4, t-runif(1),t)
y_total <-f_jit(data$y)  
rm(yy)
rm(xx1)
rm(xx2)
rm(intercept)
rm(nn)
rm(tau1)
yy=y_total
xx1=data$x1
xx2=data$x2
tau1=tau
nn=n      
intercept<-rep(1,length(y_total))
      
#fitting BQR for the whole data (or no hurdle model)

      cat("  

        model{
    # Likelihood
    
    for(i in 1:nn){
 
    mu[i] <- beta[1]*intercept[i] +beta[2]*xx1[i] + beta[3]*xx2[i]+ (1 - 2*tau1)*w[i]/(tau1*(1 - tau1))
    prec[i] <- delta*tau1*(1 - tau1)/(2*w[i])
    yy[i] ~ dnorm(mu[i],prec[i])

    }
    
    # prior distributions
    
    for(i in 1:3){
    
    beta[i] ~ dnorm(0, 0.001)
  
    }
    delta ~ dgamma(0.1, 0.1)
    
    for(i in 1:nn){w[i] ~ dexp(delta)}
    
    }
    ", file="sim.txt")


data_sim <- list("yy", "nn", "intercept", "xx1", "xx2","tau1")
sim.keeps <- c("beta")      

    
rm(bqr)
bqr<-jags(model.file = "sim.txt", data = data_sim,
              n.iter = 10000, n.chains = 3, n.thin =90, n.burnin = 1000,
               parameters.to.save = sim.keeps, DIC =TRUE)
      
rm(Bestimates)  
Bestimates<-bqr$BUGSoutput$summary
      
#Extraction of the estimates of parameters         
rm(Beta) 
Beta <-round(Bestimates[,1],5)
      
#Computation of prediction errors      
prederrors[i, ] <- (Beta[1] + Beta[2]*x1 +
                              Beta[3]*x2)-y_total    

#Computaion of the width of credible intervals        
widcred[i, ] <- c(round(Bestimates[2,c(3,7)],5),round(Bestimates[3,c(3,7)],5))
      
#Comupation of the MSE      
p1<- (-.2-Beta[2])^2
p2<-(0-Beta[3])^2      
MSE[i,]<-(p1+p2)/2
      

      
      
#Finding the quantile level of interest in the data withoout zero      
x_tau=quantile(y_total,tau)
vec = data0$y
percentiles = quantile(vec,seq(0.01,0.98,0.01))
k=seq(0.01,0.98,0.01)
tau0=k[which(abs(x_tau-percentiles)==min(abs(x_tau-percentiles)))]


#Finding the quantile level of interest in the data withoout values<4    
vec = data4$y
percentiles = quantile(vec,seq(0.01,0.98,0.01))
k=seq(0.01,0.98,0.01)
tau4=k[which(abs(x_tau-percentiles)==min(abs(x_tau-percentiles))) ]


#Fitting BQR with hurdle at zero (or hurdle model at 0)
f_jit_0<-function(t) ifelse(t>0 & t<4, t-runif(1),t)
rm(yy)
rm(xx1)
rm(xx2)
rm(intercept)
rm(nn)
rm(tau1)
yy=data0$y
xx1=data0$x1
xx2=data0$x2
tau1=tau0[length(tau0)]
intercept<-rep(1,length(data0$y))
nn=length(data0$y)  
   
rm(bqr0)
bqr0<-jags(model.file = "sim.txt", data = data_sim,
              n.iter = 10000, n.chains = 3, n.thin =90, n.burnin = 1000,
               parameters.to.save = sim.keeps, DIC =TRUE)
      
  
  
rm(Bestimates0)  
Bestimates0<-bqr0$BUGSoutput$summary
      
#Extraction of the estimates of parameters         
rm(Beta0) 
Beta0 <-round(Bestimates0[,1],5)
      
#Computation of prediction errors          
prederrors0[[i]]  <-  (Beta0[1] + Beta0[2]*data0$x1 +
                              Beta0[3]*data0$x2)-f_jit_0(data0$y)   

#Computaion of the width of credible intervals        
widcred0[i, ] <- c(round(Bestimates0[2,c(3,7)],5),round(Bestimates0[3,c(3,7)],5))
      
#Comupation of the MSE        
p1<- (-.2-Beta0[2])^2
p2<-(0-Beta0[3])^2      
MSE0[i,]<-(p1+p2)/2
      

#Fitting BQR with hurdle at 4 (or hurdle model at 4)
rm(yy)
rm(xx1)
rm(xx2)
rm(intercept)
rm(nn)
rm(tau1)
yy=data4$y
xx1=data4$x1
xx2=data4$x2
tau1=tau4[length(tau4)]
intercept<-rep(1,length(data4$y))
nn=length(data4$y)     
    

 
rm(bqr4)
bqr4<-jags(model.file = "sim.txt", data = data_sim,
                n.iter = 10000, n.chains = 3, n.thin =90, n.burnin = 1000,
               parameters.to.save = sim.keeps, DIC =TRUE)
      
      
rm(Bestimates4)   
Bestimates4<-bqr4$BUGSoutput$summary
      
#Extraction of the estimates of parameters         
rm(Beta4)  
Beta4 <-round(Bestimates4[,1],5)
      
#Computation of prediction errors          
prederrors4[[i]]  <- (Beta4[1] + Beta4[2]*data4$x1 +
                              Beta4[3]*data4$x2)-data4$y

#Computaion of the width of credible intervals      
widcred4[i, ] <- c(round(Bestimates4[2,c(3,7)],5),round(Bestimates4[3,c(3,7)],5))
      
#Comupation of the MSE        
p1<- (-.2-Beta4[2])^2
p2<-(0-Beta4[3])^2      
MSE4[i,]<-(p1+p2)/2
            
        
    
    }  
      #Returning the results
    
     list(w = widcred, w0 = widcred0,w4 = widcred4,pred = prederrors, pred0 = prederrors0,pred4 = prederrors4,MSE=MSE,MSE0=MSE0,MSE4=MSE4)
}
   


#Making a data frame for the prediction errors of the three models when n=500 and tau=0.85  

#Computaion of the prediction errors related to all three models
q85n500<-sim(quant = 0.85, n = 500, M =100)
save(q85n500, file ="q85n500.RData") 

#Extraction of the prediction errors related to BQR for the whole data
vec_85_500=c(q85n500$pred)

#Converting the list of prediction errors related to BQR with hurdle at zero to a matrix 
cbind.fill0 <- function(...){
  nm <- q85n500$pred0
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

                        
rm(k)
k=c(cbind.fill0() )
vec_0_85_500=k[which(k!="NA")]

 
                        
#Converting the list of prediction errors related to BQR with hurdle at 4 to a matrix                         
cbind.fill4 <- function(...){
  nm <- q85n500$pred4
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

                        
rm(k)
k=c(cbind.fill4() )
vec_4_85_500=k[which(k!="NA")]
                                       

#Making a data frame for the prediction errors of the three models for n=500 and tau=0.85                        
T=length(vec_85_500)+length(vec_0_85_500)+length(vec_4_85_500)                         
pred_errors_85_500 <- data.frame(Prediction_Error = c(vec_85_500,vec_0_85_500,vec_4_85_500),
                       Quantile_Level = rep(c("0.85"), each = T),
                       Sample_Size = rep(c("500"),each= T),
                       Model= rep(c("No Hurdle", "Hurdle at 0", "Hurdle at 4"), c(length(vec_85_500),length(vec_0_85_500), length(vec_4_85_500)))
                             
                             )

#pred_errors_85_500                        
                        

#Making a data frame for all prediction errors related to all quantiles in {0.75, 0.85, 0.90, 0.95}
#and the sample size of 500
total_pred_500 <- rbind(pred_errors_75_500,pred_errors_85_500,pred_errors_90_500,pred_errors_95_500)
save(total_pred_500, file = "total_pred_500.RData")


#The above process should be repeated to provide a data ftame for all quantiles in {0.75, 0.85, 0.90, 0.95}
#and all sample sizes in {500,1000,3000}
load("total_pred_3000.RData")
load("total_pred_1000.RData")
load("total_pred_500.RData")

total_pred_error <- rbind(total_pred_500,total_pred_1000,total_pred_3000)
save(total_pred_error, file = "total_pred_error.RData")


#Drawing figure 1 in the article for the comparison of the prediction errors of the models

#png(file ="prediction_error.png", width =7500 , height = 3000, units = "px", res = 800)
library(ggplot2)
options(repr.plot.width =17, repr.plot.height = 6)
ggplot(total_pred_error,aes(x = Sample_Size, y = Prediction_Error, fill = Model)) + 
  geom_boxplot(lwd=.3,outlier.size=.005)+
scale_fill_manual(values = c("azure1","lightblue3", "steelblue")) +
facet_wrap(~ Quantile_Level,ncol=4) +
  ylim(c(-10,30)) +
theme(axis.text.y = element_text( size=10),
     axis.text.x = element_text(size = 7)
        ,axis.title=element_text(size=10),
     axis.title.y = element_text(size = 14),
       axis.title.x = element_text(size = 14)
      ,legend.title = element_text( size = 12),
  legend.text = element_text( size = 6)
     )+

labs(x = "Sample size",y ="Prediction errors"  ) 
#dev.off()





#Extraction of the results for the the mean square errors (MSE) of the coefficients based on three models

#Making a data frame of the MSE for n=1000 and the quantiles of 0.75 to 0.99 (of the whole data)

q75n500<-sim(quant = 0.75, n = 1000, M =100)
q76n500<-sim(quant = 0.76, n = 1000, M =100)
q77n500<-sim(quant = 0.77, n = 1000, M =100)
q78n500<-sim(quant = 0.78, n = 1000, M =100)
q79n500<-sim(quant = 0.79, n = 1000, M =100)
q80n500<-sim(quant = 0.80, n = 1000, M =100)
q81n500<-sim(quant = 0.81, n = 1000, M =100)
q82n500<-sim(quant = 0.82, n = 1000, M =100)
q83n500<-sim(quant = 0.83, n = 1000, M =100)
q84n500<-sim(quant = 0.84, n = 1000, M =100)
q85n500<-sim(quant = 0.85, n = 1000, M =100)
q86n500<-sim(quant = 0.86, n = 1000, M =100)
q87n500<-sim(quant = 0.87, n = 1000, M =100)
q88n500<-sim(quant = 0.88, n = 1000, M =100)
q89n500<-sim(quant = 0.89, n = 1000, M =100)
q90n500<-sim(quant = 0.90, n = 1000, M =100)
q91n500<-sim(quant = 0.91, n = 1000, M =100)
q92n500<-sim(quant = 0.92, n = 1000, M =100)
q93n500<-sim(quant = 0.93, n = 1000, M =100)
q94n500<-sim(quant = 0.94, n = 1000, M =100)
q95n500<-sim(quant = 0.95, n = 1000, M =100)
q96n500<-sim(quant = 0.96, n = 1000, M =100)
q97n500<-sim(quant = 0.97, n = 1000, M =100)
q98n500<-sim(quant = 0.98, n = 1000, M =100)
q99n500<-sim(quant = 0.99, n = 1000, M =100)

MSE_1000 <- data.frame(MSE = c(  
                                  q75n1000$MSE,
                                  q75n1000$MSE0,
                                  q75n1000$MSE4,
    
                                  q76n1000$MSE,
                                  q76n1000$MSE0,
                                  q76n1000$MSE4,
    
    
                                  q77n1000$MSE,
                                  q77n1000$MSE0,
                                  q77n1000$MSE4,
    
                                  q78n1000$MSE,
                                  q78n1000$MSE0,
                                  q78n1000$MSE4,
    
                                  q79n1000$MSE,
                                  q79n1000$MSE0,
                                  q79n1000$MSE4,
    
                                  q80n1000$MSE,
                                  q80n1000$MSE0,
                                  q80n1000$MSE4,
    
                                  q81n1000$MSE,
                                  q81n1000$MSE0,
                                  q81n1000$MSE4,
    
    
                                  q82n1000$MSE,
                                  q82n1000$MSE0,
                                  q82n1000$MSE4,
    
                                  q83n1000$MSE,
                                  q83n1000$MSE0,
                                  q83n1000$MSE4,
    
                                  q84n1000$MSE,
                                  q84n1000$MSE0,
                                  q84n1000$MSE4,
    
                                  q85n1000$MSE,
                                  q85n1000$MSE0,
                                  q85n1000$MSE4,
    
                                  q86n1000$MSE,
                                  q86n1000$MSE0,
                                  q86n1000$MSE4,
    
    
                                  q87n1000$MSE,
                                  q87n1000$MSE0,
                                  q87n1000$MSE4,
    
    
                                  q88n1000$MSE,
                                  q88n1000$MSE0,
                                  q88n1000$MSE4,
    
                                  q89n1000$MSE,
                                  q89n1000$MSE0,
                                  q89n1000$MSE4,
    
    
                                  q90n1000$MSE,
                                  q90n1000$MSE0,
                                  q90n1000$MSE4,
    
    
                                  q91n1000$MSE,
                                  q91n1000$MSE0,
                                  q91n1000$MSE4,
    
                                  q92n1000$MSE,
                                  q92n1000$MSE0,
                                  q92n1000$MSE4,
    
    
                                  q93n1000$MSE,
                                  q93n1000$MSE0,
                                  q93n1000$MSE4,
    
    
                                  q94n1000$MSE,
                                  q94n1000$MSE0,
                                  q94n1000$MSE4,
    
                                  q95n1000$MSE,
                                  q95n1000$MSE0,
                                  q95n1000$MSE4,
    
    
                                  q96n1000$MSE,
                                  q96n1000$MSE0,
                                  q96n1000$MSE4,
    
                                  q97n1000$MSE,
                                  q97n1000$MSE0,
                                  q97n1000$MSE4,
    
                                  q98n1000$MSE,
                                  q98n1000$MSE0,
                                  q98n1000$MSE4,
    
    
                                  q99n1000$MSE,
                                  q99n1000$MSE0,
                                  q99n1000$MSE4
    
    
                                ),
                                  
    
 Model = rep(c("No Hurdle", "Hurdle at 0", "Hurdle at 4"), each = 100),
 Quantile_Level = rep(c("0.75", "0.76", "0.77","0.78","0.79","0.80","0.81","0.82","0.83","0.84","0.85","0.86","0.87","0.88","0.89","0.90","0.91","0.92","0.93","0.94","0.95","0.96","0.97","0.98","0.99"), each = 300),
 Sample_Size=rep(c("1000"), each =7500)
                       
                       )

#Drawing figure 2 in the article for the comparison of the MSE of the models

#png(file ="MSE.png",  width =7500 , height = 4000, units = "px", res = 800)
options(repr.plot.width =14, repr.plot.height = 6)
ggplot(subset(MSE_1000, Quantile_Level %in% c("0.75", "0.76", "0.77","0.78","0.79","0.80","0.81","0.82","0.83","0.84",
                                              "0.85","0.86","0.87","0.88","0.89","0.90","0.91","0.92","0.93","0.94","0.95",
                                              "0.96","0.97","0.98","0.99")))+
  geom_boxplot(aes(x = Quantile_Level, y = MSE, fill = Model),
               position = position_dodge(0.6)) +
  scale_fill_manual(values = c("azure1","lightblue3", "steelblue")) +
  ylim(c(0,11))+
  theme(axis.text.y = element_text( size=10),
axis.text.x = element_text(size = 7)
        ,axis.title=element_text(size=10),
     axis.title.y = element_text(size = 14),
       axis.title.x = element_text(size = 14)
      ,legend.title = element_text( size = 12),
  legend.text = element_text( size = 6)
     
     )+

labs(
       x = expression(tau),
       y ="Mean squared error of the parameter estimates" ) 

#dev.off()





#Extraction of the results for the the width of credible intervals based on three models

#Making a data frame of the width of credible intervals when n=500 and quantiles are in {0.75, 0.85, 0.90, 0.95}
q75n500<-sim(quant = 0.75, n = 500, M =100)
q85n500<-sim(quant = 0.85, n = 500, M =100)
q90n500<-sim(quant = 0.90, n = 500, M =100)
q95n500<-sim(quant = 0.95, n = 500, M =100)




width_int_500 <- data.frame(Width = c(
    
                                  q75n500$w[,2]  - q75n500$w[,1],
                                  q75n500$w[,4]  - q75n500$w[,3],
                                  q75n500$w0[,2] - q75n500$w0[,1],
                                  q75n500$w0[,4] - q75n500$w0[,3],
                                  q75n500$w4[,2] - q75n500$w4[,1],
                                  q75n500$w4[,4] - q75n500$w4[,3],
    
                                  q85n500$w[,2]  - q85n500$w[,1],
                                  q85n500$w[,4]  - q85n500$w[,3],
                                  q85n500$w0[,2] - q85n500$w0[,1],
                                  q85n500$w0[,4] - q85n500$w0[,3],
                                  q85n500$w4[,2] - q85n500$w4[,1],
                                  q85n500$w4[,4] - q85n500$w4[,3],
                                 
                                  q90n500$w[,2]  - q90n500$w[,1],
                                  q90n500$w[,4]  - q90n500$w[,3],
                                  q90n500$w0[,2] - q90n500$w0[,1],
                                  q90n500$w0[,4] - q90n500$w0[,3],
                                  q90n500$w4[,2] - q90n500$w4[,1],
                                  q90n500$w4[,4] - q90n500$w4[,3],
                                 
                                  q95n500$w[,2]  - q95n500$w[,1],
                                  q95n500$w[,4]  - q95n500$w[,3],
                                  q95n500$w0[,2] - q95n500$w0[,1],
                                  q95n500$w0[,4] - q95n500$w0[,3],
                                  q95n500$w4[,2] - q95n500$w4[,1],
                                  q95n500$w4[,4] - q95n500$w4[,3]
                                 
                                
                                ),
                                  
    
                        
                        
                        Parameter = rep(c("beta1", "beta2"), each = 100),
                        Model = rep(c("No Hurdle", "Hurdle at 0", "Hurdle at 4"), each = 200),
                        Quantile_Level = rep(c("0.75", "0.85", "0.90","0.95"), each = 600),
                       Sample_Size=rep(c("500"), each =600)
                       
                       )

#width_int_500
save(width_int_500, file = "width_int_500.RData")


#The above process should be repeated to provide data frames for all sample sizes in {500,1000,3000}
load("width_int_500.RData")
load("width_int_1000.RData")
load("width_int_3000.RData")

#Providing a data frame for all quantiles and sample sizes
width_int_total <- rbind(width_int_500,width_int_1000,width_int_3000)
save(width_int_total, file = "width_int_total.RData")


#Drawing figure 3 in the article for the comparison of the width of credible intervals of the models

#png(file ="Sim_Width_beta1.png",  width =7500 , height = 4000, units = "px", res = 800)
options(repr.plot.width =14, repr.plot.height = 6)
ggplot( subset(width_int_total,Parameter=="beta1" & Quantile_Level %in% c("0.75", "0.85", "0.90","0.95")))+
  geom_boxplot(aes(x =Quantile_Level , y = Width, fill = Model),
               position = position_dodge(0.8)) +
  scale_fill_manual(values = c("azure1","lightblue3", "steelblue")) +
  facet_wrap(~ Sample_Size , nrow = 1) +
 ylim(c(0,0.8))+
   theme(axis.text.y = element_text( size=10),
     axis.text.x = element_text(size = 7)
        ,axis.title=element_text(size=10),
     axis.title.y = element_text(size = 14),
       axis.title.x = element_text(size = 14)
      ,legend.title = element_text( size = 12),
  legend.text = element_text( size = 6)
     
     )+

labs(x =expression(tau),
       y =expression("Width of credible intervals for" ~beta[1])  ) 

#dev.off()

#Drawing figure 4 in the article for the comparison of the width of credible intervals of the models

#png(file ="Sim_Width_b2_1.png",  width =7500 , height =4000, units = "px", res = 800)
options(repr.plot.width =14, repr.plot.height = 6)
ggplot( subset(width_int_total,Parameter=="beta2"& Quantile_Level %in% c("0.75", "0.85", "0.90","0.95")))+
  geom_boxplot(aes(x =Quantile_Level  , y = Width, fill = Model),
               position = position_dodge(0.8)) +
 scale_fill_manual(values = c("azure1","lightblue3", "steelblue")) +
  facet_wrap(~ Sample_Size , nrow = 1) +
 ylim(c(0,6.5))+
  theme(axis.text.y = element_text( size=10),
   axis.text.x = element_text(size = 7)
        ,axis.title=element_text(size=10),
     axis.title.y = element_text(size = 14),
       axis.title.x = element_text(size = 14)
      ,legend.title = element_text( size = 12),
  legend.text = element_text( size = 6)
     
     )+

labs(x = expression(tau),
       y =expression("Width of credible intervals for" ~beta[2])  ) 

#dev.off()




