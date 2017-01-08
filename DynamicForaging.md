# Dynamic Foraging Patterns in Antarctic Humpbacks
Ben Weinstein  
`r Sys.time()`  






```
## Source: local data frame [12 x 2]
## 
##    Animal max(timestamp, na.rm = T)
##     (int)                    (time)
## 1  112699       2012-06-17 03:57:31
## 2  121207       2013-05-09 18:49:37
## 3  121208       2013-02-18 07:52:02
## 4  121210       2013-05-05 07:44:41
## 5  123224       2013-05-24 12:13:18
## 6  123232       2013-09-28 07:28:23
## 7  123236       2013-03-18 11:26:09
## 8  131127       2016-07-15 00:58:36
## 9  131130       2016-04-29 17:30:06
## 10 131132       2016-05-10 12:44:39
## 11 131133       2016-07-05 10:55:56
## 12 131136       2016-06-30 10:57:06
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
    phi[i,g,t,1] <- alpha_mu[state[i,g,t-1],Month[i,g,t]] 
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
    phi[i,g,steps[i,g],1] <- alpha_mu[state[i,g,steps[i,g]-1],Month[i,g,steps[i,g]-1]] 
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
      alpha_mu[1,m] ~ dbeta(1,1)
      alpha_mu[2,m] ~ dbeta(1,1)
      
      gamma[1,m] ~ dbeta(3,2)		## gamma for state 1
      dev[m] ~ dbeta(1,1)			## a random deviate to ensure that gamma[1] > gamma[2]
      gamma[2,m] <- gamma[1,m] * dev[m]
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
##      user    system   elapsed 
##   282.510     3.305 42876.936
```



##Chains

```
##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   1486124   79.4    3886542  207.6   3886542  207.6
## Vcells 285576184 2178.8  522485063 3986.3 504553714 3849.5
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  1340845  71.7    3886542  207.6   3886542  207.6
## Vcells 45524386 347.4  417988050 3189.0 504553714 3849.5
```

![](DynamicForaging_files/figure-html/unnamed-chunk-16-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-16-2.png)<!-- -->



![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

## Parameter Summary


```
##    parameter           par       mean       lower      upper
## 1   alpha_mu alpha_mu[1,1] 0.80553460 0.715431330 0.87818705
## 2   alpha_mu alpha_mu[2,1] 0.18235604 0.109409823 0.28160951
## 3   alpha_mu alpha_mu[1,2] 0.87872364 0.813858283 0.93367984
## 4   alpha_mu alpha_mu[2,2] 0.12731278 0.073274089 0.19868436
## 5   alpha_mu alpha_mu[1,3] 0.88007005 0.787044480 0.95250817
## 6   alpha_mu alpha_mu[2,3] 0.08185960 0.042158792 0.13550712
## 7   alpha_mu alpha_mu[1,4] 0.76690594 0.590762258 0.94740034
## 8   alpha_mu alpha_mu[2,4] 0.16500575 0.036409764 0.32261311
## 9   alpha_mu alpha_mu[1,5] 0.73757399 0.554322677 0.89560749
## 10  alpha_mu alpha_mu[2,5] 0.08263174 0.035565821 0.14608234
## 11     gamma    gamma[1,1] 0.92250893 0.869250266 0.97087141
## 12     gamma    gamma[2,1] 0.16287621 0.034265188 0.31577771
## 13     gamma    gamma[1,2] 0.81961440 0.766039089 0.87324607
## 14     gamma    gamma[2,2] 0.09809393 0.008796683 0.22109778
## 15     gamma    gamma[1,3] 0.76710440 0.679690136 0.85766302
## 16     gamma    gamma[2,3] 0.09415885 0.006748417 0.23430990
## 17     gamma    gamma[1,4] 0.87473536 0.738201689 0.97461891
## 18     gamma    gamma[2,4] 0.09003360 0.004889299 0.23695565
## 19     gamma    gamma[1,5] 0.74700645 0.577310171 0.91112600
## 20     gamma    gamma[2,5] 0.36712146 0.165950817 0.54295450
## 21     theta      theta[1] 0.02537279 0.006994574 0.04293103
## 22     theta      theta[2] 3.10665509 2.807198675 3.39612180
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

### As single timeline

![](DynamicForaging_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

#Simulated tracks

![](DynamicForaging_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

##Behavioral description

###Proportion of states by month

![](DynamicForaging_files/figure-html/unnamed-chunk-29-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-29-2.png)<!-- -->

###Distance between bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-30-2.png)<!-- -->

#Behavior duration


```
## Source: local data frame [156 x 7]
## Groups: Animal, Track, Bout, phistate [156]
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
## Source: local data frame [165 x 7]
## Groups: Animal, Track, Bout, phistate [165]
## 
##    Animal Track  Bout               phistate MonthF       Days    Month
##     (dbl) (dbl) (int)                 (fctr)  (dbl)      (dbl)   (fctr)
## 1       1     6     1 Area-restricted Search      1  2.6031250  January
## 2       1     9     1 Area-restricted Search      2  3.3496065 February
## 3       2     1     2 Area-restricted Search      2  0.1554977 February
## 4       2     2     2 Area-restricted Search      3 10.6119097    March
## 5       2     2     4 Area-restricted Search      4  2.4937847    April
## 6       2     2     6 Area-restricted Search      4  3.9760880    April
## 7       2     2     8 Area-restricted Search      4  0.9698958    April
## 8       2     2    10 Area-restricted Search      4  0.4637847    April
## 9       2     2    12 Area-restricted Search      4  0.9947569    April
## 10      2     2    14 Area-restricted Search      4  0.9829745    April
## ..    ...   ...   ...                    ...    ...        ...      ...
```

```
## Source: local data frame [321 x 7]
## Groups: Animal, Track, Bout, phistate [321]
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

![](DynamicForaging_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

##Proportion of time allocation
![](DynamicForaging_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

```
##      Month Traveling Area-restricted Search     PropF TotalTime
## 1  January 136.27752              147.59362 0.5199318 283.87115
## 2 February 223.28979              341.49608 0.6046470 564.78587
## 3    March 109.96427              192.57084 0.6365239 302.53512
## 4    April  50.09941              162.72390 0.7645962 212.82331
## 5      May  21.40273               60.49271 0.7386578  81.89544
```

## Number of bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

#Time spent in grid cell
## All years
![](DynamicForaging_files/figure-html/unnamed-chunk-34-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-3.png)<!-- -->

##Add in sea ice

![](DynamicForaging_files/figure-html/unnamed-chunk-35-1.png)<!-- -->


