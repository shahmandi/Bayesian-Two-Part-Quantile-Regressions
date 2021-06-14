#Required packages for MCMC computation based on Hamiltonian Markov Chain sampler 
library(rstan)
library(bayesplot)

#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
#input="Arts_Humanities.csv"
#input="Architecture.csv"
input="Emergency_Nursing.csv"
#input="Literature_Literary.csv"
#input="Media_Technology.csv"
#input="Religious_Studies.csv"
#input="Visual_Arts_Performing.csv"

df1 <-read.csv(input,sep=",",header=TRUE)

df<-df1
logic_con<-df$len_first_name >-1
df<-df[which(logic_con != FALSE),]


#Defining the independent variables
indi_collab<-as.numeric(df$len_first_name)
indi_collab<-log(indi_collab)
l_title<-as.numeric(df$len_title)
journal_inter<-as.numeric(df$gini_year_coun_list)


#Specifing the dependent variable
citation<-as.numeric(df$CitedBy)
nonzero <- which(citation != 0)

n.cit <- length(citation)
citation.star<- citation
citation.star[nonzero] <- log(citation.star[nonzero]+1 - runif(length(nonzero)))


intercept.cit <- rep(1, n.cit)



#Stan code
cat("

//Definition of the inputs imported for feeding the model
data {
  int<lower=0> N;
  int<lower=0> k;
  matrix[N, k] X;
  vector[N] y;
  real<lower=0, upper=1> tau;
}


//Definition of the parameters in the model
parameters {
 
  vector[k] beta;
  real<lower=0> delta;
  vector<lower=0>[N] w;
}


//Computation of transformed parameters
transformed parameters {

 
  
  vector[N] perc;
  vector[N] perc2;
  vector[N] mu;
  vector[N] mu2;

 
  
  perc  <- (tau * (1 - tau) * delta) ./ (2 * w);
  for(n in 1:N) perc2[n] <- inv_sqrt(perc[n]);

  mu  <- X*beta;
  mu2  <- (1 - 2 * tau) / (tau * (1 - tau)) * w + mu;
}


model {


  //Priors used for parameters

  for(i in 1:k)
  {
   beta[i]  ~ normal(0, 3);
 
  }

  delta ~ gamma(0.001,0.001);
  w ~ exponential(delta);



  //Definition of the density function of the model
  
   for (i in 1:N) {
   y[i] ~ normal(mu2[i], perc2[i]);
    }
        
      
}

" , file="stan_qr.stan")

#Prepration of data for using in stan functon
X=as.data.frame(cbind(intercept.cit,indi_collab,l_title,journal_inter))

stan_data <- list(
    N = n.cit,
    k= 4,
    X = X,
    y = citation.star,
    tau=0.50
  )


#Stan function 
fit_rstan <- stan(file = "stan_qr.stan",data = stan_data, chains = 4, warmup = 1000,iter = 2000)
fit_rstan



#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
#input="Arts_Humanities.csv"
#input="Architecture.csv"
input="Emergency_Nursing.csv"
#input="Literature_Literary.csv"
#input="Media_Technology.csv"
#input="Religious_Studies.csv"
#input="Visual_Arts_Performing.csv"

df1 <-read.csv(input,sep=",",header=TRUE)

df<-df1
logic_con<-df$len_first_name >-1
df<-df[which(logic_con != FALSE),]


#Defining the independent variables
indi_collab<-as.numeric(df$len_first_name)
indi_collab<-log(indi_collab)
l_title<-as.numeric(df$len_title)
journal_inter<-as.numeric(df$gini_year_coun_list)


#Specifing the dependent variable with the hurdle point at 0
citation<-as.numeric(df$CitedBy)
nonzero <- which(citation != 0)

n.cit <- length(citation)
citation.star<- citation
citation.star[nonzero] <- log(citation.star[nonzero] - runif(length(nonzero)))


intercept.cit <- rep(1, n.cit)



#Stan code
cat("

//Definition of the inputs imported for feeding the model
data {
  int<lower=0> N;
  int<lower=0> k;
  matrix[N, k] X;
  vector[N] y;
  real<lower=0, upper=1> tau;
}


//Definition of the parameters in the model
parameters {
 
  vector[k] beta;
  vector[k] gammaa;
  real<lower=0> delta;
  vector<lower=0>[N] w;
}


//Computation of transformed parameters
transformed parameters {

  vector[N] perc;
  vector[N] perc2;
  vector<lower=0, upper=1>[N] pi;

  vector[N] mu;
  vector[N] mu2;
  

 
 
  perc  <- (tau * (1 - tau) * delta) ./ (2 * w);
  for(n in 1:N) perc2[n] <- inv_sqrt(perc[n]);

  for (n in 1:N) pi[n]<- inv_logit( X[n]*gammaa);

  mu  <- X*beta;
  mu2  <- (1 - 2 * tau) / (tau * (1 - tau)) * w + mu;
}


model {


  //Priors used for parameters

  for(i in 1:k)
  {
   beta[i]  ~ normal(0, 3);
   gammaa[i]  ~ normal(0, 3);
  }

  delta ~ gamma(0.001,0.001);
  w ~ exponential(delta);



  //Definition of the mass probability function of the model
  
   for (i in 1:N) {
   
   (y[i] == 0) ~ bernoulli(1-pi[i]);
    if (y[i] != 0) y[i] ~ normal(mu2[i], perc2[i]);
    }
        
      
}

//Posterior predictive checks

generated quantities {

 vector[N] y_rep;
 vector[N] y_rep_zer0;
 vector[N] y_zer0;

 real<lower = 0, upper = 1> mean_gt;
 real<lower = 0, upper = 1> max_gt;
 real<lower = 0, upper = 1> sd_gt;
 real<lower = 0, upper = 1> zero_percent;

//Simulation data for doing posterior predictive checks
  for (n in 1:N) {

    if (bernoulli_rng(1-pi[n]))
   {
      y_rep[n] = 0;
    }

   else {
  
    y_rep[n] = normal_rng(mu2[n], perc2[n]);
 

     }
}


//Providing p-value for some statistics for posterior predictive checks

    mean_gt = mean(y_rep) >= mean(y);
    max_gt = max(y_rep) >= max(y);
    sd_gt = sd(y_rep) >= sd(y);



  for (n in 1:N) {

  y_rep_zer0[n] = y_rep[n] ==0 ? 1.0 : 0.0;
  y_zer0[n] = y[n] ==0 ? 1.0 : 0.0;
}
 
  zero_percent =(sum(y_rep_zer0)/N) >= (sum(y_zer0)/N);
  


}

" , file="stan_tpqr.stan")


#Prepration of data for using in stan functon
X_h0=as.data.frame(cbind(intercept.cit,indi_collab,l_title,journal_inter))

stan_data_h0 <- list(
    N = n.cit,
    k= 4,
    X = X_h0,
    y = citation.star,
    tau=0.50
  )


#Stan function 
fit_rstan_h0 <- stan(file = "stan_tpqr.stan",data = stan_data_h0, chains = 4, warmup = 1000,iter = 2000)
fit_rstan_h0



#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
#input="Arts_Humanities.csv"
#input="Architecture.csv"
input="Emergency_Nursing.csv"
#input="Literature_Literary.csv"
#input="Media_Technology.csv"
#input="Religious_Studies.csv"
#input="Visual_Arts_Performing.csv"

df1 <-read.csv(input,sep=",",header=TRUE)

df<-df1
logic_con<-df$len_first_name >-1
df<-df[which(logic_con != FALSE),]


#Defining the independent variables
indi_collab<-as.numeric(df$len_first_name)
indi_collab<-log(indi_collab)
l_title<-as.numeric(df$len_title)
journal_inter<-as.numeric(df$gini_year_coun_list)


#Specifing the dependent variable with the hurdle point at 3
citation<-as.numeric(df$CitedBy)
nonzero <- which(citation > 3)
zero <- which(citation <=3)

n.cit <- length(citation)
citation.star<- citation
citation.star[zero] <- rep(0,length(zero))
citation.star[nonzero] <- log(citation.star[nonzero]-3 - runif(length(nonzero)))


intercept.cit <- rep(1, n.cit)



#Stan code
cat("

//Definition of the inputs imported for feeding the model
data {
  int<lower=0> N;
  int<lower=0> k;
  matrix[N, k] X;
  vector[N] y;
  real<lower=0, upper=1> tau;
}


//Definition of the parameters in the model
parameters {
 
  vector[k] beta;
  vector[k] gammaa;
  real<lower=0> delta;
  vector<lower=0>[N] w;
}


//Computation of transformed parameters
transformed parameters {

  vector[N] perc;
  vector[N] perc2;
  vector<lower=0, upper=1>[N] pi;

  vector[N] mu;
  vector[N] mu2;
  

 
 
  perc  <- (tau * (1 - tau) * delta) ./ (2 * w);
  for(n in 1:N) perc2[n] <- inv_sqrt(perc[n]);

  for (n in 1:N) pi[n]<- inv_logit( X[n]*gammaa);

  mu  <- X*beta;
  mu2  <- (1 - 2 * tau) / (tau * (1 - tau)) * w + mu;
}


model {


  //Priors used for parameters

  for(i in 1:k)
  {
   beta[i]  ~ normal(0, 3);
   gammaa[i]  ~ normal(0, 3);
  }

  delta ~ gamma(0.001,0.001);
  w ~ exponential(delta);



  //Definition of the mass probability function of the model
  
   for (i in 1:N) {
   
   (y[i] == 0) ~ bernoulli(1-pi[i]);
    if (y[i] != 0) y[i] ~ normal(mu2[i], perc2[i]);
    }


}


//Posterior predictive checks

generated quantities {

 vector[N] y_rep;
 vector[N] y_rep_zer0;
 vector[N] y_zer0;

 real<lower = 0, upper = 1> mean_gt;
 real<lower = 0, upper = 1> max_gt;
 real<lower = 0, upper = 1> sd_gt;
 real<lower = 0, upper = 1> zero_percent;

//Simulation data for doing posterior predictive checks
  for (n in 1:N) {

    if (bernoulli_rng(1-pi[n]))
   {
      y_rep[n] = 0;
    }

   else {
  
    y_rep[n] = normal_rng(mu2[n], perc2[n]);
 

     }
}


//Providing p-value for some statistics for posterior predictive checks

    mean_gt = mean(y_rep) >= mean(y);
    max_gt = max(y_rep) >= max(y);
    sd_gt = sd(y_rep) >= sd(y);



  for (n in 1:N) {

  y_rep_zer0[n] = y_rep[n] ==0 ? 1.0 : 0.0;
  y_zer0[n] = y[n] ==0 ? 1.0 : 0.0;
}
 
  zero_percent =(sum(y_rep_zer0)/N) >= (sum(y_zer0)/N);
  


}

" , file="stan_tpqr.stan")



#Prepration of data for using in stan functon
X_h3=as.data.frame(cbind(intercept.cit,indi_collab,l_title,journal_inter))

stan_data_h3 <- list(
    N = n.cit,
    k= 4,
    X = X_h3,
    y = citation.star,
    tau=0.50
  )



#Stan function 
fit_rstan_h3 <- stan(file = "stan_tpqr.stan",data = stan_data_h3, chains = 4, warmup = 1000,iter = 2000)
fit_rstan_h3
