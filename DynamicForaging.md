# Dynamic Foraging Patterns in Antarctic Humpbacks
Ben Weinstein  
`r Sys.time()`  







![](DynamicForaging_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

##By Month

![](DynamicForaging_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

#Correlated random walk

*Process Model*

$$ d_{t} \sim T*d_{t-1} + Normal(0,\Sigma)$$
$$ x_t = x_{t-1} + d_{t} $$

## Parameters

For each individual:

$$\theta = \text{Mean turning angle}$$
$$\gamma = \text{Move persistence} $$

For both behaviors process variance is:
$$ \sigma_{latitude} = 0.1$$
$$ \sigma_{longitude} = 0.1$$

##Behavioral States

$$ \text{For each individual i}$$
$$ Behavior_1 = \text{traveling}$$
$$ Behavior_2 = \text{foraging}$$

$$ \alpha_{i,1,1} = \text{Probability of remaining traveling when traveling}$$
$$\alpha_{i,2,1} = \text{Probability of switching from Foraging to traveling}$$

$$\begin{matrix}
  \alpha_{i,1,1} & 1-\alpha_{i,1,1} \\
  \alpha_{i,2,1} & 1-\alpha_{i,2,1} \\
\end{matrix}$$

With the probability of switching states:

$$logit(\phi_{traveling}) = \alpha_{Behavior_{t-1}}$$

$$\phi_{foraging} = 1 - \phi_{traveling} $$

##Continious tracks

The transmitter will often go dark for 10 to 12 hours, due to weather, right in the middle of an otherwise good track. The model requires regular intervals to estimate the turning angles and temporal autocorrelation. As a track hits one of these walls, call it the end of a track, and begin a new track once the weather improves. We can remove any micro-tracks that are less than three days.
Specify a duration, calculate the number of tracks and the number of removed points. Iteratively.







How did the filter change the extent of tracks?

![](DynamicForaging_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


![](DynamicForaging_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


![](DynamicForaging_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-10-2.png)<!-- -->



sink("Bayesian/Multi_RW.jags")
cat("
    model{
    
    #Constants
    pi <- 3.141592653589
    
    ##argos observation error##
    argos_prec[1:2,1:2] <- inverse(argos_sigma*argos_cov[,])
    
    #Constructing the covariance matrix
    argos_cov[1,1] <- 1
    argos_cov[1,2] <- sqrt(argos_alpha) * rho
    argos_cov[2,1] <- sqrt(argos_alpha) * rho
    argos_cov[2,2] <- argos_alpha
    
    for(i in 1:ind){
    for(g in 1:tracks[i]){
    
    ## Priors for first true location
    #for lat long
    y[i,g,1,1:2] ~ dmnorm(argos[i,g,1,1,1:2],argos_prec)
    
    #First movement - random walk.
    y[i,g,2,1:2] ~ dmnorm(y[i,g,1,1:2],iSigma)
    
    ###First Behavioral State###
    state[i,g,1] ~ dcat(lambda[]) ## assign state for first obs
    
    #Process Model for movement
    for(t in 2:(steps[i,g]-1)){
    
    #Behavioral State at time T
    logit(phi[i,g,t,1]) <- alpha_mu[state[i,g,t-1],Month[i,g,t]] 
    phi[i,g,t,2] <- 1-phi[i,g,t,1]
    state[i,g,t] ~ dcat(phi[i,g,t,])
    
    #Turning covariate
    #Transition Matrix for turning angles
    T[i,g,t,1,1] <- cos(theta[state[i,g,t]])
    T[i,g,t,1,2] <- (-sin(theta[state[i,g,t]]))
    T[i,g,t,2,1] <- sin(theta[state[i,g,t]])
    T[i,g,t,2,2] <- cos(theta[state[i,g,t]])
    
    #Correlation in movement change
    d[i,g,t,1:2] <- y[i,g,t,] + gamma[state[i,g,t],Month[i,g,t]] * T[i,g,t,,] %*% (y[i,g,t,1:2] - y[i,g,t-1,1:2])
    
    #Gaussian Displacement
    y[i,g,t+1,1:2] ~ dmnorm(d[i,g,t,1:2],iSigma)
    }
    
    #Final behavior state
    logit(phi[i,g,steps[i,g],1]) <- alpha_mu[state[i,g,steps[i,g]-1],Month[i,g,steps[i,g]-1]] 
    phi[i,g,steps[i,g],2] <- 1-phi[i,g,steps[i,g],1]
    state[i,g,steps[i,g]] ~ dcat(phi[i,g,steps[i,g],])
    
    ##	Measurement equation - irregular observations
    # loops over regular time intervals (t)    
    
    for(t in 2:steps[i,g]){
    
    # loops over observed locations within interval t
    for(u in 1:idx[i,g,t]){ 
    zhat[i,g,t,u,1:2] <- (1-j[i,g,t,u]) * y[i,g,t-1,1:2] + j[i,g,t,u] * y[i,g,t,1:2]
    
    #for each lat and long
    #argos error
    argos[i,g,t,u,1:2] ~ dmnorm(zhat[i,g,t,u,1:2],argos_prec)
    }
    }
    }
    }
    ###Priors###
    
    #Process Variance
    iSigma ~ dwish(R,2)
    Sigma <- inverse(iSigma)
    
    ##Mean Angle
    tmp[1] ~ dbeta(10, 10)
    tmp[2] ~ dbeta(10, 10)
    
    # prior for theta in 'traveling state'
    theta[1] <- (2 * tmp[1] - 1) * pi
    
    # prior for theta in 'foraging state'    
    theta[2] <- (tmp[2] * pi * 2)
    
    ##Move persistance
    # prior for gamma (autocorrelation parameter) in state 1

    #for each month
    for (m in 1:Months){

    #Intercepts
    alpha_mu[1,m] ~ dnorm(0,0.386)
    alpha_mu[2,m] ~ dnorm(0,0.386)
    
    gamma[2,m] ~ dbeta(1.5, 2)		## gamma for state 2
    dev[m] ~ dbeta(1,1)			## a random deviate to ensure that gamma[1] > gamma[2]
    gamma[1,m] <- gamma[2,m] + dev[m] 		## gamma for state 1
    }
    
    
    ##Behavioral States
    
    #Hierarchical structure across motnhs
    
    #Variance
    alpha_tau[1] ~ dt(0,1,1)I(0,)
    alpha_tau[2] ~ dt(0,1,1)I(0,)
    
    #Probability of behavior switching 
    lambda[1] ~ dbeta(1,1)
    lambda[2] <- 1 - lambda[1]
    
    ##Argos priors##
    #longitudinal argos error
    argos_sigma ~ dunif(0,10)
    
    #latitidunal argos error
    argos_alpha~dunif(0,10)
    
    #correlation in argos error
    rho ~ dunif(-1, 1)
    
    
    }"
    ,fill=TRUE)
sink()






##Chains



![](DynamicForaging_files/figure-html/unnamed-chunk-15-1.png)<!-- -->





![](DynamicForaging_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->


## Parameter Summary

```
##    parameter           par         mean        lower       upper
## 1   alpha_mu alpha_mu[1,1]  0.435094526 -0.027983034  0.86747797
## 2   alpha_mu alpha_mu[2,1] -1.531759690 -2.185339756 -0.89387142
## 3   alpha_mu alpha_mu[1,2] -0.899195730 -1.872449788  0.45653067
## 4   alpha_mu alpha_mu[2,2] -1.430390490 -2.013469353 -0.91411079
## 5   alpha_mu alpha_mu[1,3]  0.697841455 -0.507002105  1.82881099
## 6   alpha_mu alpha_mu[2,3] -2.216104974 -3.746215000 -0.27281394
## 7   alpha_mu alpha_mu[1,4] -0.459418904 -1.204279471  0.11424481
## 8   alpha_mu alpha_mu[2,4] -2.074815920 -2.712398660 -1.50041297
## 9   alpha_mu alpha_mu[1,5] -3.387112559 -4.902844836 -2.20772464
## 10  alpha_mu alpha_mu[2,5] -0.405077298 -0.850387704  0.09181157
## 11  alpha_mu alpha_mu[1,6] -0.310181591 -2.907868690  2.35334270
## 12  alpha_mu alpha_mu[2,6] -0.311501561 -2.943725284  2.41500589
## 13     gamma    gamma[1,1]  0.940776404  0.858265027  1.02593535
## 14     gamma    gamma[2,1]  0.070803844  0.008176977  0.16242277
## 15     gamma    gamma[1,2]  1.021385925  0.817772902  1.20975967
## 16     gamma    gamma[2,2]  0.265492246  0.193066046  0.33423481
## 17     gamma    gamma[1,3]  0.797019814  0.693174712  0.92580335
## 18     gamma    gamma[2,3]  0.184816648  0.072069243  0.27359696
## 19     gamma    gamma[1,4]  1.237517268  1.136851683  1.33537207
## 20     gamma    gamma[2,4]  0.278530068  0.203796415  0.35430711
## 21     gamma    gamma[1,5]  1.519547169  1.438301436  1.58066183
## 22     gamma    gamma[2,5]  0.575766210  0.534346731  0.61982973
## 23     gamma    gamma[1,6]  0.727464417  0.219082027  1.30468416
## 24     gamma    gamma[2,6]  0.311104925  0.053555864  0.67452841
## 25     theta      theta[1] -0.001467907 -0.044871431  0.04110741
## 26     theta      theta[2]  3.120116044  0.047367175  6.20141242
```

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

#Behavioral Prediction




##Autocorrelation in behavior

![](DynamicForaging_files/figure-html/unnamed-chunk-22-1.png)<!-- -->


#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-23-1.png)<!-- -->


##Behavioral description

## Predicted behavior duration




![](DynamicForaging_files/figure-html/unnamed-chunk-25-1.png)<!-- -->


## Duration by month

![](DynamicForaging_files/figure-html/unnamed-chunk-26-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-26-2.png)<!-- -->




#Time spent in grid cell

![](DynamicForaging_files/figure-html/unnamed-chunk-28-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-28-2.png)<!-- -->



