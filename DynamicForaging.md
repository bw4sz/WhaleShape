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
      logit(alpha[1,m]) <- alpha_mu_invlogit[1,m]
      alpha_mu_invlogit[1,m] ~ dnorm(alpha_mu[1],alpha_tau[1])
      
      logit(alpha[2,m]) <- alpha_mu_invlogit[2,m]
      alpha_mu_invlogit[2,m] ~ dnorm(alpha_mu[2],alpha_tau[2])
      
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
##      user    system   elapsed 
##   337.338    14.828 97254.086
```



##Chains

```
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells    1490299    79.6    3205452   171.2    3205452   171.2
## Vcells 1714633394 13081.7 4032004763 30761.8 3457122303 26375.8
```

```
##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells   1330370   71.1    3205452   171.2    3205452   171.2
## Vcells 146016614 1114.1 3225603810 24609.5 3457122303 26375.8
```

![](DynamicForaging_files/figure-html/unnamed-chunk-16-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-16-2.png)<!-- -->



![](DynamicForaging_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

## Change in autocorrelation over time

![](DynamicForaging_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# Change in transition probabilities over time

![](DynamicForaging_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

## Parameter Summary


```
##    parameter          par         mean         lower        upper
## 1      alpha   alpha[1,1]   0.83802019   0.761248126   0.90612053
## 2      alpha   alpha[2,1]   0.14317890   0.084073800   0.21495309
## 3      alpha   alpha[1,2]   0.87296859   0.807546737   0.92467307
## 4      alpha   alpha[2,2]   0.12249091   0.071840545   0.18916900
## 5      alpha   alpha[1,3]   0.86305844   0.786676613   0.92539809
## 6      alpha   alpha[2,3]   0.08023776   0.045901482   0.12386840
## 7      alpha   alpha[1,4]   0.81645308   0.668781346   0.93813263
## 8      alpha   alpha[2,4]   0.11736292   0.047220365   0.20650919
## 9      alpha   alpha[1,5]   0.78277762   0.656528137   0.88449402
## 10     alpha   alpha[2,5]   0.07049474   0.037936518   0.11173555
## 11     alpha   alpha[1,6]   0.76821363   0.619891322   0.88894542
## 12     alpha   alpha[2,6]   0.10854330   0.045065706   0.19556304
## 13  alpha_mu  alpha_mu[1]   1.58447391   1.007321285   2.15959609
## 14  alpha_mu  alpha_mu[2]  -2.16539180  -2.694175247  -1.64174412
## 15 alpha_tau alpha_tau[1]   3.34698452   1.283325269   4.85325825
## 16 alpha_tau alpha_tau[2]   3.33273827   1.274820878   4.86538508
## 17     gamma   gamma[1,1]   0.91296259   0.862549491   0.96496316
## 18     gamma   gamma[2,1]   0.17617872   0.037863949   0.32390327
## 19     gamma   gamma[1,2]   0.83134324   0.780860273   0.88737858
## 20     gamma   gamma[2,2]   0.08535385   0.006840522   0.20630532
## 21     gamma   gamma[1,3]   0.80865222   0.735332620   0.88794686
## 22     gamma   gamma[2,3]   0.07734333   0.005983412   0.19668083
## 23     gamma   gamma[1,4]   0.86514312   0.754764822   0.99033609
## 24     gamma   gamma[2,4]   0.10780413   0.006218684   0.29597768
## 25     gamma   gamma[1,5]   0.81909457   0.716411572   0.92365736
## 26     gamma   gamma[2,5]   0.38637429   0.210061579   0.55143390
## 27     gamma   gamma[1,6]   0.86634345   0.755787930   0.98500252
## 28     gamma   gamma[2,6]   0.38447246   0.110143467   0.65417152
## 29  gamma_mu     gamma_mu   0.84739652   0.783185851   0.91265592
## 30 gamma_tau    gamma_tau 298.64687846 122.076514596 477.39130746
## 31     theta     theta[1]   0.02440948   0.007619716   0.04175364
## 32     theta     theta[2]   3.07400912   2.855138705   3.29013549
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
## Source: local data frame [166 x 7]
## Groups: Animal, Track, Bout, phistate [166]
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
## Source: local data frame [177 x 7]
## Groups: Animal, Track, Bout, phistate [177]
## 
##    Animal Track  Bout               phistate MonthF       Days    Month
##     (dbl) (dbl) (int)                 (fctr)  (dbl)      (dbl)   (fctr)
## 1       1     6     1 Area-restricted Search      1  2.6031250  January
## 2       1     9     1 Area-restricted Search      2  3.3496065 February
## 3       2     1     2 Area-restricted Search      2  0.1564352 February
## 4       2     1     4 Area-restricted Search      2  0.1554977 February
## 5       2     2     2 Area-restricted Search      3 11.1066551    March
## 6       2     2     4 Area-restricted Search      3  0.9433218    March
## 7       2     2     6 Area-restricted Search      4  1.9837269    April
## 8       2     2     8 Area-restricted Search      4  3.9760880    April
## 9       2     2    10 Area-restricted Search      4  0.1050231    April
## 10      2     2    12 Area-restricted Search      4  0.9947569    April
## ..    ...   ...   ...                    ...    ...        ...      ...
```

```
## Source: local data frame [343 x 7]
## Groups: Animal, Track, Bout, phistate [343]
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
## 1  January 121.93544              142.51822 0.5389157 264.45366
## 2 February 222.33936              331.34093 0.5984337 553.68029
## 3    March 110.77988              236.50453 0.6810111 347.28441
## 4    April  44.79596              150.25373 0.7703356 195.04969
## 5      May  24.35963               82.92464 0.7729431 107.28427
## 6     June  17.38661               44.18263 0.7176088  61.56924
```

## Number of bouts

![](DynamicForaging_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

#Time spent in grid cell
## All years
![](DynamicForaging_files/figure-html/unnamed-chunk-34-1.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-2.png)<!-- -->![](DynamicForaging_files/figure-html/unnamed-chunk-34-3.png)<!-- -->

##Add in sea ice

![](DynamicForaging_files/figure-html/unnamed-chunk-35-1.png)<!-- -->


