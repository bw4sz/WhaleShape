traj<-function(gamma=gamma,theta=theta,a1,total_time=total_time,step_length){
  
  ##Constants
  #time interval in days, divided into 4 hours intervals
  steps<-length(seq(0,total_time*24,by=step_length))
  
  #generate observations along time interval
  timestep<-sort(runif(300,0,total_time) * 24)  
  
  #Position Vector
  xy<-matrix(nrow=steps,ncol=2)
  
  #Displacement Vector
  d<-matrix(nrow=steps,ncol=2)
  
  #Behavioral States (at time t)
  state<-c()
  
  #Probability of staying in behavior
  phi<-matrix(nrow=steps,ncol=2)
  
  #Mean turning angle
  theta=theta
  
  #Degree of autocorrelation
  gamma=gamma
  
  #Process variance in latitude
  sigma_lat=0.2
  
  #Process variance in longitude
  sigma_lon=0.2
  
  #Correlation in process variance
  rho=0
  
  #Multivariate Normal Variance in random walk
  Sigma<-matrix(nrow=2,ncol=2)
  Sigma[1,1] <- sigma_lon^2
  Sigma[1,2] <- rho * sigma_lon * sigma_lat
  Sigma[2,1] <- rho * sigma_lon * sigma_lat
  Sigma[2,2] <- sigma_lat^2
  
  #Transition matrix for correlated movement
  T<-array(dim=c(steps,2,2))
  
  ##Process Model
  
  ##Initial position##
  xy[1,]<-mvrnorm(mu=c(0,0),Sigma=Sigma)
  
  #First behavioral state is a random draw from two probabilities
  lambda<-c()
  lambda[1]<-0.5
  lambda[2]<-1-lambda[1]
  state[1] <- sample(c(1,2),size=1,prob=lambda)
  
  #First step is random walk
  xy[2,]<-mvrnorm(mu=xy[1,],Sigma = Sigma)
  
  #Draw random walk locations
  for (x in 2:(steps-1)){
    
    #Extract env

    #Behavior
    phi[x,1] <- a1[state[x-1]]
    phi[x,2]<- 1 - phi[x,1]
    state[x] <- sample(c(1,2),size=1,prob=phi[x,])

    #Movement correlation matrix
    T[x,1,1] <- cos(theta[state[x]])
    T[x,1,2] <- -sin(theta[state[x]])
    T[x,2,1] <- sin(theta[state[x]])
    T[x,2,2] <- cos(theta[state[x]])
    
    # Add Correlated Displacement
    d[x,]<-xy[x,] + gamma[state[x]] * T[x,,] %*% (xy[x,]  - xy[x-1,])
    
    #next position
    ## Random walk
    xy[x+1,]<-mvrnorm(n=1,mu=d[x,],Sigma=Sigma)
  }
  
  #Format to data frame
  dxy<-data.frame(xy)
  colnames(dxy)<-c("x","y")
  dxy$Step<-1:nrow(dxy)
  
  #Behavior
  dxy$State<-as.factor(c(state,NA))
  
  levels(dxy$State)<-c("Traveling","Feeding")
  
  #add time label
  dxy$Hour<-seq(0,total_time*24,by=step_length)
  
  ##Measurement model
  argosdf<-list()
  
  #first location has no interpolation error
  argos_x<-dxy[dxy$Step==1,c("x")]
  argos_y<-dxy[dxy$Step==1,c("y")]
  argosdf[[1]]<-data.frame(Step=1,time=0.5,argos_x,argos_y)
  
  for (x in 2:steps){
    
    #which time slices sit in interval
    argos_time<-timestep[step_length * (x-1) < timestep & timestep < step_length * x]
    
    #step locations
    trueloc<-dxy[dxy$Step==x,c("x","y")]
    pastloc<-dxy[dxy$Step==(x-1),c("x","y")]
    
    #for each of those observations
    stepdf<-list()
    for(i in 1:length(argos_time)){
      
      #interpolation distance
      j<-(argos_time[i]-step_length*(x-1))/step_length
      
      if(length(argos_time)==0){j=0}
      #observed locations
      #Add argos noise.
      argos_x<-(1-j)*pastloc$x+j*trueloc$x
      argos_y<-(1-j)*pastloc$y+j*trueloc$y
      
      #order dataframe, if there are no observations, fill with NA
      if(!length(argos_time)==0){
        stepdf[[i]]<-data.frame(Step=x,time=argos_time[i],argos_x,argos_y)
      } else{
        stepdf[[i+1]]<-data.frame(Step=x,time=NA,argos_x,argos_y)
      }
    }
    argosdf[[x]]<-rbind_all(stepdf)
  }
  #bind df together
  argosdf<-rbind_all(argosdf)
  
  #merge with true positions
  dxy<-merge(dxy,argosdf,by="Step")
  
  return(dxy)
}