# Dynamic Foraging Patterns in Antarctic Humpbacks
Ben Weinstein  
`r Sys.time()`  







![](DynamicForaging_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

##By Month

![](DynamicForaging_files/figure-html/unnamed-chunk-5-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

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
    argos_prec[1:2,1:2] <- argos_cov[,]
    
    #Constructing the covariance matrix
    argos_cov[1,1] <- argos_sigma
    argos_cov[1,2] <- 0
    argos_cov[2,1] <- 0
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
    
    gamma[1,m] ~ dunif(0.6,1)		## gamma for state 1
    dev[m] ~ dbeta(1,1)			## a random deviate to ensure that gamma[1] > gamma[2]
    gamma[2,m] <- gamma[1,m] * dev[m]	## gamma for state 2
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
    #longitudinal argos precision
    argos_sigma <- 5
    
    #latitidunal argos precision
    argos_alpha <- 5
    
    
    }"
    ,fill=TRUE)
sink()


```
##      user    system   elapsed 
##   674.409     2.705 22880.234
```



##Chains

```
##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   1568602   83.8    3886542  207.6   3886542  207.6
## Vcells 239998185 1831.1  441914955 3371.6 437088820 3334.8
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  1350285  72.2    3886542  207.6   3886542  207.6
## Vcells 35925323 274.1  353531964 2697.3 437088820 3334.8
```

![](DynamicForaging_files/figure-html/unnamed-chunk-16-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-16-2.png)<!-- -->



![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

## Parameter Summary


```
##    parameter           par        mean      lower      upper
## 1   alpha_mu alpha_mu[1,1]  2.48290522  1.8585392  3.2463030
## 2   alpha_mu alpha_mu[2,1] -2.66885559 -3.4460343 -2.0431672
## 3   alpha_mu alpha_mu[1,2]  2.70848838  2.2056528  3.1782343
## 4   alpha_mu alpha_mu[2,2] -2.87544456 -3.4780071 -2.3658373
## 5   alpha_mu alpha_mu[1,3]  2.92797063  2.4807533  3.4157086
## 6   alpha_mu alpha_mu[2,3] -2.77074067 -3.3047394 -2.2945489
## 7   alpha_mu alpha_mu[1,4]  2.61059720  1.9778024  3.1796421
## 8   alpha_mu alpha_mu[2,4] -2.04206169 -2.6087956 -1.5256137
## 9   alpha_mu alpha_mu[1,5]  3.14189186  2.3914147  3.9980522
## 10  alpha_mu alpha_mu[2,5] -2.69244776 -3.5800998 -1.8674712
## 11     gamma    gamma[1,1]  0.81705394  0.7670661  0.8665163
## 12     gamma    gamma[2,1]  0.80132037  0.7552084  0.8524941
## 13     gamma    gamma[1,2]  0.71572924  0.6809255  0.7510245
## 14     gamma    gamma[2,2]  0.70588691  0.6722358  0.7406841
## 15     gamma    gamma[1,3]  0.64159390  0.6050942  0.6863813
## 16     gamma    gamma[2,3]  0.62373720  0.5845515  0.6675901
## 17     gamma    gamma[1,4]  0.76001234  0.7054404  0.8206015
## 18     gamma    gamma[2,4]  0.74336108  0.6853136  0.8034810
## 19     gamma    gamma[1,5]  0.86067795  0.8087098  0.9106875
## 20     gamma    gamma[2,5]  0.84126630  0.7929746  0.8887574
## 21     theta      theta[1]  0.01473992 -2.9905291  3.0011089
## 22     theta      theta[2]  6.24539132  6.2300745  6.2597500
```

![](DynamicForaging_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

#Behavioral Prediction



![](DynamicForaging_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

##Spatial Prediction

![](DynamicForaging_files/figure-html/unnamed-chunk-24-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-24-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-24-3.png)<!-- -->

## By individual

![](DynamicForaging_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

##Autocorrelation in behavior

![](DynamicForaging_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

##Behavioral description

## Predicted behavior duration



![](DynamicForaging_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

## Duration by month

![](DynamicForaging_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-2.png)<!-- -->



## Duration by month

![](DynamicForaging_files/figure-html/unnamed-chunk-32-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-32-2.png)<!-- -->

#Proportion of states by month

![](DynamicForaging_files/figure-html/unnamed-chunk-33-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-33-2.png)<!-- -->

#Time between foraging bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-34-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-3.png)<!-- -->

#Distance between bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

#Time spent in grid cell

![](DynamicForaging_files/figure-html/unnamed-chunk-36-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-36-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-36-3.png)<!-- -->


