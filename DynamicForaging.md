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



### After filitering

![](DynamicForaging_files/figure-html/unnamed-chunk-7-1.png)<!-- -->





How did the filter change the extent of tracks?

![](DynamicForaging_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


![](DynamicForaging_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


![](DynamicForaging_files/figure-html/unnamed-chunk-11-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-11-2.png)<!-- -->



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
    tmp[1] ~ dbeta(20, 20)
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
    
    gamma[2,m] ~ dbeta(1, 4)		## gamma for state 2
    dev[m] ~ dunif(0.3,1)			## a random deviate to ensure that gamma[1] > gamma[2]
    gamma[1,m] <- gamma[2,m] + dev[m]	## gamma for state 1
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


```
##      user    system   elapsed 
##   611.695     3.248 55290.091
```



##Chains

```
##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   1507267   80.5    3886542  207.6   3886542  207.6
## Vcells 150279946 1146.6  239101918 1824.3 239101680 1824.3
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  1336847  71.4    3886542  207.6   3886542  207.6
## Vcells 33713406 257.3  191281534 1459.4 239101680 1824.3
```


![](DynamicForaging_files/figure-html/unnamed-chunk-16-1.png)<!-- -->




![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->


## Parameter Summary

```
##    parameter           par       mean       lower      upper
## 1   alpha_mu alpha_mu[1,1] -2.7954064 -4.51389698 -1.2244744
## 2   alpha_mu alpha_mu[2,1] -1.0916653 -1.63281590 -0.5204797
## 3   alpha_mu alpha_mu[1,2] -1.7171449 -2.68769383 -0.7718806
## 4   alpha_mu alpha_mu[2,2] -1.2011110 -1.65941252 -0.7223104
## 5   alpha_mu alpha_mu[1,3] -2.2466673 -3.62417519 -1.1121026
## 6   alpha_mu alpha_mu[2,3] -1.8798462 -2.58036813 -1.2886152
## 7   alpha_mu alpha_mu[1,4] -1.9903005 -3.17355817 -0.9051239
## 8   alpha_mu alpha_mu[2,4] -1.2616994 -1.83810428 -0.6585734
## 9   alpha_mu alpha_mu[1,5] -1.2239859 -2.57249695 -0.2128475
## 10  alpha_mu alpha_mu[2,5] -1.8603937 -2.59863744 -1.2110516
## 11  alpha_mu alpha_mu[1,6]  0.1504398 -0.86500282  1.2811342
## 12  alpha_mu alpha_mu[2,6] -2.1600366 -3.53337135 -0.8016734
## 13     gamma    gamma[1,1]  1.2517096  1.11438382  1.3483729
## 14     gamma    gamma[2,1]  0.3209919  0.25616917  0.3865315
## 15     gamma    gamma[1,2]  1.0548418  0.89627087  1.1752491
## 16     gamma    gamma[2,2]  0.2574881  0.21317204  0.3011035
## 17     gamma    gamma[1,3]  1.0850166  0.93680570  1.1881968
## 18     gamma    gamma[2,3]  0.1700221  0.10883371  0.2290166
## 19     gamma    gamma[1,4]  1.2069696  1.07947349  1.2949223
## 20     gamma    gamma[2,4]  0.2844567  0.21786109  0.3475187
## 21     gamma    gamma[1,5]  1.2606688  1.16143356  1.3633228
## 22     gamma    gamma[2,5]  0.3022215  0.22358145  0.3893407
## 23     gamma    gamma[1,6]  1.0581631  0.90697742  1.1858615
## 24     gamma    gamma[2,6]  0.1409226  0.04346489  0.2345640
## 25     theta      theta[1]  0.0784530  0.03983875  0.1158941
## 26     theta      theta[2]  6.2076677  6.17456864  6.2387553
```

![](DynamicForaging_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

#Behavioral Prediction




![](DynamicForaging_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

## Confidence
![](DynamicForaging_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

## By individual

![](DynamicForaging_files/figure-html/unnamed-chunk-25-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

##Autocorrelation in behavior

![](DynamicForaging_files/figure-html/unnamed-chunk-26-1.png)<!-- -->


#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


##Behavioral description

## Predicted behavior duration




![](DynamicForaging_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


## Duration by month

![](DynamicForaging_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-2.png)<!-- -->




#Time spent in grid cell

![](DynamicForaging_files/figure-html/unnamed-chunk-32-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-32-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-32-3.png)<!-- -->



