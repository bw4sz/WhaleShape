a<-traj(gamma=c(1,0.2),theta=c(0,6),a1=inv.logit(c(-2,-1.5)),total_time = 500,step_length=12)
ggplot(a,aes(x=x,y=y)) + theme_bw()  + geom_point(aes(col=State),size=1)  + theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) + labs(x="",y="") + geom_path(aes(col=State,group=1),size=0.5) 
