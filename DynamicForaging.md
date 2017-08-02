# Dynamic Foraging Patterns in Antarctic Humpbacks
Ben Weinstein  
`r Sys.time()`  






```
## Source: local data frame [12 x 2]
## 
##    Animal max(timestamp, na.rm = T)
##     (int)                    (time)
## 1  112699       2012-06-17 03:57:31
## 2  121207       2013-05-09 18:49:00
## 3  121208       2013-02-18 07:52:00
## 4  121210       2013-05-05 07:44:00
## 5  123224       2013-05-24 12:13:00
## 6  123232       2013-09-28 07:28:00
## 7  123236       2013-03-18 11:26:00
## 8  131127       2016-07-15 07:58:36
## 9  131130       2016-04-30 00:30:06
## 10 131132       2016-05-10 19:44:39
## 11 131133       2016-07-05 20:26:44
## 12 131136       2016-06-30 18:49:35
```

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
    
    #for each if 6 argos class observation error
    
    for(x in 1:6){
    
    ##argos observation error##
    argos_prec[x,1:2,1:2] <- argos_cov[x,,]
    
    #Constructing the covariance matrix
    argos_cov[x,1,1] <- argos_sigma[x]
    argos_cov[x,1,2] <- 0
    argos_cov[x,2,1] <- 0
    argos_cov[x,2,2] <- argos_alpha[x]
    }
    
    for(i in 1:ind){
    for(g in 1:tracks[i]){
    
    ## Priors for first true location
    #for lat long
    y[i,g,1,1:2] ~ dmnorm(argos[i,g,1,1,1:2],argos_prec[1,1:2,1:2])
    
    #First movement - random walk.
    y[i,g,2,1:2] ~ dmnorm(y[i,g,1,1:2],iSigma)
    
    ###First Behavioral State###
    state[i,g,1] ~ dcat(lambda[]) ## assign state for first obs
    
    #Process Model for movement
    for(t in 2:(steps[i,g]-1)){
    
    #Behavioral State at time T
    phi[i,g,t,1] <- alpha[state[i,g,t-1],Month[i,g,t]] 
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
    phi[i,g,steps[i,g],1] <- alpha[state[i,g,steps[i,g]-1],Month[i,g,steps[i,g]-1]] 
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
    argos[i,g,t,u,1:2] ~ dmnorm(zhat[i,g,t,u,1:2],argos_prec[argos_class[i,g,t,u],1:2,1:2])
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
      alpha[1,m] ~ dbeta(1,1)

      alpha[2,m] ~ dbeta(1,1)
      
      
      gamma[1,m] ~ dnorm(gamma_mu,gamma_tau)		## gamma for state 1
      dev[m] ~ dbeta(1,1)			## a random deviate to ensure that gamma[1] > gamma[2]
      gamma[2,m] <- gamma[1,m] * dev[m]
    }
    
    ##Behavioral States
    
    #Hierarchical structure across months
    
    #Switching among states in inv.logit space
    alpha_mu[1] ~ dnorm(0,0.386)
    alpha_mu[2] ~ dnorm(0,0.386)

    #Variance in state change per month
    alpha_tau[1] ~ dunif(0,5)
    alpha_tau[2] ~ dunif(0,5)
    
    #Probability of behavior switching 
    lambda[1] ~ dbeta(1,1)
    lambda[2] <- 1 - lambda[1]

    #spatial autocorrelation priors, we know the true state is high autocorrelation
    gamma_mu~dnorm(0.8,100)

    #reasonable variance keeps the chains from wandering.
    gamma_tau ~ dunif(100,500)
    
    ##Argos priors##
    #longitudinal argos precision, from Jonsen 2005, 2016, represented as precision not sd
    
    #by argos class
    argos_sigma[1] <- 11.9016
    argos_sigma[2] <- 10.2775
    argos_sigma[3] <- 1.228984
    argos_sigma[4] <- 2.162593
    argos_sigma[5] <- 3.885832
    argos_sigma[6] <- 0.0565539
    
    #latitidunal argos precision, from Jonsen 2005, 2016
    argos_alpha[1] <- 67.12537
    argos_alpha[2] <- 14.73474
    argos_alpha[3] <- 4.718973
    argos_alpha[4] <- 0.3872023
    argos_alpha[5] <- 3.836444
    argos_alpha[6] <- 0.1081156
    
    
    }"
    ,fill=TRUE)
sink()


```
##       user     system    elapsed 
##    334.943     13.404 115238.677
```



##Check Gelman Rubin



##Chains

```
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells    1491782    79.7    3205452   171.2    3205452   171.2
## Vcells 1714636332 13081.7 4032004780 30761.8 3457182926 26376.3
```

```
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells    1491975    79.7    3205452   171.2    3205452   171.2
## Vcells 1835596070 14004.5 4032004780 30761.8 3457182926 26376.3
```

![](DynamicForaging_files/figure-html/unnamed-chunk-17-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-17-2.png)<!-- -->



![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-21-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-21-2.png)<!-- -->

## Parameter Summary


```
##    parameter          par         mean         lower        upper
## 1      alpha   alpha[1,1]   0.80272609   0.709781302   0.88275731
## 2      alpha   alpha[2,1]   0.18417316   0.109886082   0.28596504
## 3      alpha   alpha[1,2]   0.85779609   0.769663126   0.92651316
## 4      alpha   alpha[2,2]   0.14277837   0.074662328   0.23242453
## 5      alpha   alpha[1,3]   0.84809249   0.756291365   0.92646800
## 6      alpha   alpha[2,3]   0.08231248   0.041061798   0.13314298
## 7      alpha   alpha[1,4]   0.71974044   0.552429628   0.91616349
## 8      alpha   alpha[2,4]   0.18694099   0.062148052   0.33259023
## 9      alpha   alpha[1,5]   0.72416208   0.553740143   0.86467130
## 10     alpha   alpha[2,5]   0.07342636   0.033334537   0.13045811
## 11     alpha   alpha[1,6]   0.64232640   0.444304376   0.82654295
## 12     alpha   alpha[2,6]   0.17339837   0.060332138   0.34555909
## 13  alpha_mu  alpha_mu[1]   0.02486589  -2.725346674   2.69872887
## 14  alpha_mu  alpha_mu[2]  -0.01029971  -2.662084854   2.60340546
## 15 alpha_tau alpha_tau[1]   2.51146610   0.245123592   4.77690829
## 16 alpha_tau alpha_tau[2]   2.52500326   0.258015652   4.76923839
## 17     gamma   gamma[1,1]   0.92994797   0.874182895   0.99815805
## 18     gamma   gamma[2,1]   0.15544185   0.019085260   0.31478130
## 19     gamma   gamma[1,2]   0.84353009   0.787270981   0.90370236
## 20     gamma   gamma[2,2]   0.07979238   0.005575408   0.20071307
## 21     gamma   gamma[1,3]   0.82789181   0.749898718   0.90894382
## 22     gamma   gamma[2,3]   0.06339788   0.004223703   0.16736089
## 23     gamma   gamma[1,4]   0.92005535   0.786832814   1.04616744
## 24     gamma   gamma[2,4]   0.07579794   0.004179654   0.23179100
## 25     gamma   gamma[1,5]   0.84989878   0.741810880   0.96045772
## 26     gamma   gamma[2,5]   0.36474752   0.187056218   0.53299185
## 27     gamma   gamma[1,6]   0.90467546   0.793598913   1.03364088
## 28     gamma   gamma[2,6]   0.36325160   0.091319031   0.65163358
## 29  gamma_mu     gamma_mu   0.87402906   0.806961545   0.94251461
## 30 gamma_tau    gamma_tau 287.28682070 116.488343696 476.90033397
## 31     theta     theta[1]   0.02360181   0.005872277   0.04146352
## 32     theta     theta[2]   3.04564301   2.780647705   3.28326050
```

![](DynamicForaging_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

#Behavioral Prediction



![](DynamicForaging_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

##Spatial Prediction

![](DynamicForaging_files/figure-html/unnamed-chunk-25-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-25-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

## By individual

![](DynamicForaging_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

##Autocorrelation in behavior

![](DynamicForaging_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

### As single timeline

![](DynamicForaging_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

##Behavioral description

###Proportion of states by month

![](DynamicForaging_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-2.png)<!-- -->

###Distance between bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-31-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-31-2.png)<!-- -->

#Behavior duration


```
## Source: local data frame [177 x 7]
## Groups: Animal, Track, Bout, phistate [177]
## 
##    Animal Track  Bout  phistate MonthF      Days    Month
##     (dbl) (dbl) (int)    (fctr)  (dbl)     (dbl)   (fctr)
## 1       1     1     1 Traveling      1 0.6059028  January
## 2       1     2     1 Traveling      1 0.4953819  January
## 3       1     3     1 Traveling      1 0.6483912  January
## 4       1     4     1 Traveling      1 1.5586921  January
## 5       1     5     1 Traveling      1 0.6500116  January
## 6       1     7     1 Traveling      2 0.9866435 February
## 7       1     8     1 Traveling      2 0.5460301 February
## 8       1    10     1 Traveling      2 0.9215972 February
## 9       1    11     1 Traveling      2 2.6218171 February
## 10      1    12     1 Traveling      2 0.6120602 February
## ..    ...   ...   ...       ...    ...       ...      ...
```

```
## Source: local data frame [193 x 7]
## Groups: Animal, Track, Bout, phistate [193]
## 
##    Animal Track  Bout               phistate MonthF       Days    Month
##     (dbl) (dbl) (int)                 (fctr)  (dbl)      (dbl)   (fctr)
## 1       1     6     1 Area-restricted Search      1  2.6031250  January
## 2       1     9     1 Area-restricted Search      2  3.3496065 February
## 3       2     1     2 Area-restricted Search      2  0.1513079 February
## 4       2     1     4 Area-restricted Search      2  0.1564352 February
## 5       2     1     6 Area-restricted Search      2  0.1554977 February
## 6       2     1     8 Area-restricted Search      2  0.1646644 February
## 7       2     2     2 Area-restricted Search      3  0.1447917    March
## 8       2     2     4 Area-restricted Search      3 11.9918056    March
## 9       2     2     6 Area-restricted Search      3  1.9623843    March
## 10      2     2     8 Area-restricted Search      3  1.9876620    March
## ..    ...   ...   ...                    ...    ...        ...      ...
```

```
## Source: local data frame [370 x 7]
## Groups: Animal, Track, Bout, phistate [370]
## 
##    Animal Track  Bout               phistate MonthF      Days    Month
##     (dbl) (dbl) (int)                 (fctr)  (dbl)     (dbl)   (fctr)
## 1       1     1     1              Traveling      1 0.6059028  January
## 2       1     2     1              Traveling      1 0.4953819  January
## 3       1     3     1              Traveling      1 0.6483912  January
## 4       1     4     1              Traveling      1 1.5586921  January
## 5       1     5     1              Traveling      1 0.6500116  January
## 6       1     6     1 Area-restricted Search      1 2.6031250  January
## 7       1     7     1              Traveling      2 0.9866435 February
## 8       1     8     1              Traveling      2 0.5460301 February
## 9       1     9     1 Area-restricted Search      2 3.3496065 February
## 10      1    10     1              Traveling      2 0.9215972 February
## ..    ...   ...   ...                    ...    ...       ...      ...
```

![](DynamicForaging_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

##Proportion of time allocation
![](DynamicForaging_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```
##      Month Traveling Area-restricted Search     PropF TotalTime
## 1  January 116.31559              139.96898 0.5461467 256.28457
## 2 February 222.67517              359.08926 0.6172417 581.76443
## 3    March  83.89773              216.83172 0.7210192 300.72946
## 4    April  50.95119              166.57176 0.7657664 217.52295
## 5      May  22.33618               83.91296 0.7897754 106.24914
## 6     June  13.01757               49.52854 0.7918724  62.54611
```

## Number of bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

#Time spent in grid cell
## All years
![](DynamicForaging_files/figure-html/unnamed-chunk-35-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-35-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-35-3.png)<!-- -->

##Add in sea ice

![](DynamicForaging_files/figure-html/unnamed-chunk-36-1.png)<!-- -->


