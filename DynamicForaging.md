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
    
    gamma[2,m] ~ dunif(0, 0.2)		## gamma for state 2
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


```
##      user    system   elapsed 
##   438.256     3.764 43020.526
```



##Chains

```
##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   1513774   80.9    3886542  207.6   3886542  207.6
## Vcells 262821566 2005.2  538564282 4109.0 459450749 3505.4
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  1336189  71.4    3886542  207.6   3886542  207.6
## Vcells 40760370 311.0  430851425 3287.2 459450749 3505.4
```


![](DynamicForaging_files/figure-html/unnamed-chunk-16-1.png)<!-- -->





![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->


## Parameter Summary

```
##    parameter           par         mean        lower       upper
## 1   alpha_mu alpha_mu[1,1] -0.480673909 -1.366390743  0.28765934
## 2   alpha_mu alpha_mu[2,1] -0.890403782 -1.592645942 -0.16333616
## 3   alpha_mu alpha_mu[1,2] -1.398069525 -2.321879979 -0.55987737
## 4   alpha_mu alpha_mu[2,2] -0.833781616 -1.378258259 -0.19830592
## 5   alpha_mu alpha_mu[1,3] -2.143812997 -4.055915694  0.17619683
## 6   alpha_mu alpha_mu[2,3] -1.698160854 -2.877896124 -0.60579066
## 7   alpha_mu alpha_mu[1,4] -0.582959093 -1.305516506  0.05108164
## 8   alpha_mu alpha_mu[2,4] -1.356080259 -1.961044375 -0.79422693
## 9   alpha_mu alpha_mu[1,5] -0.893545829 -1.677189647 -0.18919132
## 10  alpha_mu alpha_mu[2,5] -0.037775431 -0.877001299  0.82578452
## 11  alpha_mu alpha_mu[1,6]  1.530829476  0.503866298  2.77062700
## 12  alpha_mu alpha_mu[2,6] -2.396903309 -3.842947689 -1.14575636
## 13  alpha_mu alpha_mu[1,7] -0.140954037 -2.820613339  2.35749472
## 14  alpha_mu alpha_mu[2,7]  0.021237198 -2.553351588  2.68780998
## 15     gamma    gamma[1,1]  0.975804024  0.819368763  1.13423934
## 16     gamma    gamma[2,1]  0.159215080  0.086058275  0.19796555
## 17     gamma    gamma[1,2]  0.986314539  0.860990493  1.13710077
## 18     gamma    gamma[2,2]  0.188343055  0.162888534  0.19946081
## 19     gamma    gamma[1,3]  1.032419323  0.804344610  1.16747616
## 20     gamma    gamma[2,3]  0.151735949  0.092795615  0.19541140
## 21     gamma    gamma[1,4]  1.149043639  1.087534217  1.18950862
## 22     gamma    gamma[2,4]  0.181924201  0.146514728  0.19891005
## 23     gamma    gamma[1,5]  1.177099936  1.138841628  1.19687216
## 24     gamma    gamma[2,5]  0.196867425  0.190209160  0.19983664
## 25     gamma    gamma[1,6]  0.936349862  0.840419856  1.04156468
## 26     gamma    gamma[2,6]  0.083238672  0.008531231  0.17810235
## 27     gamma    gamma[1,7]  0.560177940  0.125345899  1.06058426
## 28     gamma    gamma[2,7]  0.092193465  0.008693266  0.18579578
## 29     theta      theta[1]  0.004479108 -0.023755537  0.03428136
## 30     theta      theta[2]  6.161656754  6.099857547  6.20940066
```

![](DynamicForaging_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

#Behavioral Prediction



#Behavior Prediction
![](DynamicForaging_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

##Autocorrelation in behavior

![](DynamicForaging_files/figure-html/unnamed-chunk-24-1.png)<!-- -->


#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-25-1.png)<!-- -->


##Behavioral description

## Predicted behavior duration




![](DynamicForaging_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


## Duration by month

![](DynamicForaging_files/figure-html/unnamed-chunk-28-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-28-2.png)<!-- -->




#Time spent in grid cell

![](DynamicForaging_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-3.png)<!-- -->



