# Bayesian Hurdle Quantile Regressions

The R code is for the article [**"A Bayesian Hurdle Quantile Regression Model for Citation Analysis with Mass Points at Lower Values"**](https://direct.mit.edu/qss/article/doi/10.1162/qss_a_00147/103156/A-Bayesian-Hurdle-Quantile-Regression-Model-for). This research focuses on two challenges for the analysis of citation counts by quantile regression (QR): _discontinuity_ and _substantial mass points at lower counts_, such as zero, one, two, and three. In this research, an update of a Bayesian two-part QR model was introduced to scientometrics to address these problems. The usefulness and applicability of the method were illustrated based on both simulated and real citation count data.

## Simulation stydy
For simulation study, samples with sizes of 500, 1000 and 3000 are simulated from continuous log-normal distribution (LN) with mean(2−0.2∗x1 +0∗x2 +ε) and standard deviation 0.4 where x1 ∼LN(2,2), x2 ∼ N (0.5, 0.5), and ε ∼ N (0, 1). The log-normal distribution was chosen because it approximates the typical distribution of citation counts. The floor function was used for the simulated values less than 4 to simulate substantial mass points at 0, 1, 2 and 3. The intercept value of 2 and the coefficient −0.2 of x1 were chosen so that approximately 45% of the data will be zeros and 75% of the data less than 3, similar to much citation count data. Bayesian QR and Bayesian two-part QR models with hurdle at 0, and with hurdle at 3 will be fitted. The objective is to compare Bayesian QR with the QR parts of the two-part models with hurdles at 0 and 3. In the following, for each sample size and for each quantile level

 1. The Bayesian QR model is fitted to the whole data. 
 2. Then the quantile level of the corresponding quantile value is found in the data in which the zeros are excluded and the Bayesian QR model is fitted, (i.e. the QR part of the two-part model with hurdle at 0).
 3.  Next, the quantile value corresponding to the quantile level is found in the data in which all substantial mass points (including 0,1,2, and 3) are removed and the Bayesian QR model is fitted for the corresponding quantile (this model is the QR part of the two-part model with hurdle at 3). 
 
 This process is repeated 100 times for each sample size, and for 4 specific quantiles of 0.75, 0.85, 0.90, and 0.95 (of the whole data) separately. For each combination
- The prediction error of each model
- The mean squared error of the parameters’ estimates (intercept excluded)
- The width of the credible intervals for both independent variables x1 and x2
  
 are computed. The function **jags** from the R-library **R2jags** which is based on **Gibbs sampler** was used for **MCMC** computation corresponding to the Bayesian QR models. In this function, three chains will run. For each chain, 10000 iterations with burn-in 1000 and thinning number of 90 were considered. All the R code for this simulation is available in the file of **"R_Code_Simulation.r"**.

## Study of real citation count data
The data used in this article consists of citation counts for standard journal articles (excluding reviews) published in the following seven Scopus fields: 

- Arts and Humanities (all)
- Literature and Literary Theory
-  Religious Studies
-  Visual Arts and Performing Arts
-  Media Technology, Architecture
-  Emergency Nursing

The articles were published in the year 2010 and their information was extracted at the end of the year 2019, giving the citation counts time to mature. _The number of citations to each article_ is the dependent variable. _Collaboration, length of title, and journal internationality_ were selected as independent variables. A sequence of quantiles from 0.05 to 0.95 is considered. Ordinary Bayesian QR, Bayesian two-part QR with hurdle at 0 and Bayesian two-part QR with hurdle at 3 were fitted to the datasets. **MCMC** was calculated once with **jags** function from the R-library **R2jags** and another time with **stan** function from the R-library **rstan**. The jags function is based on the Gibbs sampler and stan function is based on **Hamiltonian Monte Carlo (HMC)**. In the following, the details of both methods are explained. 

- ### Gibbs sampling algorithm

Three chains were run in the **jags** function. Based on a trial and error process for determining theproper values of iterations, burn-in and thinning rate of the MCMC chains, it was deducedthat the autocorrelation plot as one convergence indicator for the parameters correspondingto the journal internationality in both parts of the proposed model showed a slow decreasing pattern by increasing the lag number, indicating slow mixing in the chain.  To address this problem,  large  numbers  of  these  hyperparameters  are  selected.   Three  chains  separately with iteration 100000,  50000 burn-in samples,  and the thinning rate of 160 are provided to  approximate  the  posterior  distributions  of  the  parameters  and  then  their  means  areconsidered as the estimates of the parameters in the models. All the R code for this study is available in the file of **"R_code_citation_counts_R2jags.r"** and also the data of real citation counts corresponding to the 7 fields of Scopus are in the folder of **"Data_scopus"**.



- ### Hamiltonian Monte Carlo (HMC) sampling algorithm

In practice, running the Gibbs sampler with these large numbers for the hyperparameters is computationally expensive. It motivated us to move to HMC algorithm. It is a common sampler in cases that the autocorrelation plots show slow mixing and slow convergence corresponding to parameters in a model. HMC uses the derivative of the logarithm of the target posterior distribution. The benefit of the derivative is to present a clear image for the local shape of the distribution, specifying the points (areas) with higher probability density. HMC gives momentum to the Markov chain to move towards the higher parts of the density that we prefer to have more samples from. It is worth mentioning that the discontinuity of our two-part model at hurdle point $c$ does not make any problem for applying HMC because HMC is sensitive to discontinuity in the parameter space, not the data space. In the model all parameters are continuous and HMC differentiates and explores the posterior distribution of the parameters rather than the data. There are no derivatives of the data because the data are fixed. For HMC, four chains were run in the **stan** function. For each chain, 2000 iterations with warm-up 1000, and thinning rate of 2 were used. All the R code for this study is available in the file of **"R_code_citation_counts_rstan.r"**.









