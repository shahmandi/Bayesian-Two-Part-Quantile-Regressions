#Required packages for MCMC computation based on Gibbs sampler
library(R2jags)
library(jagsUI)


#------------------------Defining the Bayesian quantile regression model-----------


#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
input="Arts_Humanities.csv"
#input="Architecture.csv"
#input="Emergency_Nursing.csv"
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



#Jags code

cat("  
 
    model{
    # Definig the Bayesian quantile regression model
    
    for(i in 1:n.cit){
    
    me[i]<-(1 - 2*tau)*w[i]/(tau*(1 - tau))
    mu[i] <- beta[1]*intercept.cit[i] +beta[2]*indi_collab[i] + beta[3]*l_title[i] +
    beta[4]*journal_inter[i] + me[i]

    pcit[i] <- delta*tau*(1 - tau)/(2*w[i])

    citation.star[i] ~ dnorm(mu[i],pcit[i])

    }
    
    # prior distributions
    
    for(i in 1:4){
    
    beta[i] ~ dnorm(0, 0.1)

    }

    delta ~ dgamma(0.001, 0.001)
    
    for(i in 1:n.cit){w[i] ~ dexp(delta)}
    
    }
    ", file="jags_qr.txt")





cit.data <- list("citation.star", "n.cit", "tau", "intercept.cit", "indi_collab", "l_title", "journal_inter")

cit.keeps <- c("beta")

#Jags function
tau=0.50
fit_jags <-jags(model.file = "jags_qr.txt", data = cit.data,
              n.iter = 100000, n.chains = 3, n.thin =160, n.burnin = 50000,
               parameters.to.save = cit.keeps, DIC =TRUE)
fit_jags


#------------------------Defining the Bayesian two-part hurdle quantile regression model with hurdle at 0-----------
#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
input="Arts_Humanities.csv"
#input="Architecture.csv"
#input="Emergency_Nursing.csv"
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



#Jags code

cat("  
    data{for(i in 1:n.cit){ones[i] <- 1}
    C= 10000}
    
    # Defining the Bayesian two-part hurdle quantile regression model
    model{
    
    
    for(i in 1:n.cit){
    
    z[i] <- equals(citation.star[i], 0)        # I(citation.star = 0)


    Xgamma[i]<-gamma[1]*intercept.cit[i] + gamma[2]*indi_collab[i] +
    gamma[3]*l_title[i] + gamma[4]*journal_inter[i]
    
    # Modeling the probability of zeros
    q[i] <- 1/( 1+exp(Xgamma[i]))
    pi[i] <- max(0.001, min(0.999, q[i]))
    

    me[i]<- (1 - 2*tau)*w[i]/(tau*(1 - tau))
    mu[i] <- beta[1]*intercept.cit[i] +beta[2]*indi_collab[i] + beta[3]*l_title[i] +
    beta[4]*journal_inter[i] + me[i]


    pcit[i] <- delta*tau*(1 - tau)/(2*w[i])

    density_normal[i]<-pow(pcit[i]/(2*3.14159),0.5)*exp(- 0.5*pow(citation.star[i] - mu[i], 2)*pcit[i])
    likelihood[i] <-z[i]*(pi[i]) + (1 - z[i])*(1-pi[i])*density_normal[i]
    
    # Ones trick
    ones[i] ~ dbern(phi[i])
    phi[i] <- likelihood[i]/C
    
    }
    
    # prior distributions
    
    for(i in 1:4){
    
    beta[i] ~ dnorm(0, 0.1)
    
    gamma[i] ~ dnorm(0, 0.1)
    }
    delta ~ dgamma(0.001, 0.001)
    
    for(i in 1:n.cit){w[i] ~ dexp(delta)}


    # Posterior predictive checks

    # Simulation of data 
    for(i in 1:n.cit){

    bern1[i] ~ dbern(pi[i])
    norm1[i] ~ dnorm(mu[i],pcit[i])

    indicator[i] <- ifelse(bern1[i], 1, 0)
    y_rep[i] <- (0 * indicator[i]) + (norm1[i] * (1 - indicator[i]))
    }



    # Computation of the errors based on original data and simulated data
    for(i in 1:n.cit){

        res[i] <- ifelse(citation.star[i]==0,0,citation.star[i]- ((1-pi[i]) * mu[i]) ) #(y-yhat)
        res.new[i] <- ifelse(y_rep[i]==0,0,y_rep[i]- ((1-pi[i]) * mu[i]) )


        }

     # Derived parameters
      fit <- sum(res[])
      fit.new <- sum(res.new[])



    
    }
    ", file="jags_tpqr_zerohurdle.txt")


cit.data <- list("citation.star", "n.cit", "tau", "intercept.cit", "indi_collab", "l_title", "journal_inter")

cit.keeps <- c("beta", "gamma","y_rep","fit","fit.new")

#Jags function
tau=0.50
library(jagsUI)
fit_jags_h0 <-jags(model.file = "jags_tpqr_zerohurdle.txt", data = cit.data,
             n.iter = 100000, n.chains = 3, n.thin =160, n.burnin = 50000,
               parameters.to.save = cit.keeps, DIC =TRUE)
fit_jags_h0



#Posterior predictive check plot for the errors
pp.check(fit_jags_h0, observed = 'fit', simulated = 'fit.new')


#------------------------Defining the Bayesian two-part hurdle quantile regression model with hurdle at 3-----------

#Reading data and providing independent and dependent variables

#Defining the field of Scopus 
input="Arts_Humanities.csv"
#input="Architecture.csv"
#input="Emergency_Nursing.csv"
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



#Jags code

cat("  
    data{for(i in 1:n.cit){ones[i] <- 1}
    C= 10000}
    

    # Defining the Bayesian two-part hurdle quantile regression model
    model{
   
    
    for(i in 1:n.cit){
    
    z[i] <- equals(citation.star[i], 0)        # I(citation.star = 0)


    Xgamma[i]<-gamma[1]*intercept.cit[i] + gamma[2]*indi_collab[i] +
    gamma[3]*l_title[i] + gamma[4]*journal_inter[i]
    
    #Modeling the probability of zeros
    q[i] <- 1/( 1+exp(Xgamma[i]))
    pi[i] <- max(0.001, min(0.999, q[i]))

    me[i]<- (1 - 2*tau)*w[i]/(tau*(1 - tau))
    mu[i] <- beta[1]*intercept.cit[i] +beta[2]*indi_collab[i] + beta[3]*l_title[i] +
    beta[4]*journal_inter[i] + me[i]


    pcit[i] <- delta*tau*(1 - tau)/(2*w[i])

    density_normal[i]<-pow(pcit[i]/(2*3.14159),0.5)*exp(- 0.5*pow(citation.star[i] - mu[i], 2)*pcit[i])
    likelihood[i] <-z[i]*(pi[i]) + (1 - z[i])*(1-pi[i])*density_normal[i]
    
    #Ones trick
    ones[i] ~ dbern(phi[i])
    phi[i] <- likelihood[i]/C
    
    }
    
    # prior distributions
    
    for(i in 1:4){
    
    beta[i] ~ dnorm(0, 0.1)
    
    gamma[i] ~ dnorm(0, 0.1)
    }
    delta ~ dgamma(0.001, 0.001)
    
    for(i in 1:n.cit){w[i] ~ dexp(delta)}


    
    # Posterior predictive checks

    # Simulation of data 
    for(i in 1:n.cit){

    bern1[i] ~ dbern(pi[i])
    norm1[i] ~ dnorm(mu[i],pcit[i])

    indicator[i] <- ifelse(bern1[i], 1, 0)
    y_rep[i] <- (0 * indicator[i]) + (norm1[i] * (1 - indicator[i]))
    }



    # Computation of the errors based on original data and simulated data
    for(i in 1:n.cit){

        res[i] <- ifelse(citation.star[i]==0,0,citation.star[i]- ((1-pi[i]) * mu[i]) ) #(y-yhat)
        res.new[i] <- ifelse(y_rep[i]==0,0,y_rep[i]- ((1-pi[i]) * mu[i]) )


        }

     # Derived parameters
      fit <- sum(res[])
      fit.new <- sum(res.new[])



    
    }
    ", file="jags_tpqr_threehurdle.txt")


cit.data <- list("citation.star", "n.cit", "tau", "intercept.cit", "indi_collab", "l_title", "journal_inter")

cit.keeps <- c("beta", "gamma","y_rep","fit","fit.new")

#Jags function
tau=0.50
library(jagsUI)
fit_jags_h3 <-jags(model.file = "jags_tpqr_threehurdle.txt", data = cit.data,
              n.iter = 500, n.chains = 3, n.thin =1, n.burnin =100,
               parameters.to.save = cit.keeps, DIC =TRUE)
fit_jags_h3
#n.iter = 100000, n.chains = 3, n.thin =160, n.burnin = 50000,

#Posterior predictive check plot for the errors
pp.check(fit_jags_h3, observed = 'fit', simulated = 'fit.new')


