a<-mxy %>% filter(Animal==2,Track %in% 2,Bout==25) %>% as.data.frame()

d<-SpatialPointsDataFrame(cbind(a$x,a$y),data=a,proj4string=CRS("+proj=longlat +datum=WGS84"))
ggplot(a,aes(x=x,y=y,col=phistate)) + borders(fill="black") + coord_cartesian(xlim=as.numeric(bbox(d)[1,]),ylim=as.numeric(bbox(d)[2,])) + theme_bw() + facet_wrap(~Track,ncol=2) + geom_path(aes(group=Track),size=1)

#as bouts
ggplot(a,aes(x=x,y=y,col=as.factor(Bout))) + borders(fill="black") + coord_cartesian(xlim=as.numeric(bbox(d)[1,]),ylim=as.numeric(bbox(d)[2,])) + theme_bw() + facet_wrap(~Track,ncol=2) + geom_path(aes(group=Track),size=1)


#without tracks
a<-mxy %>% filter(Month=="June") %>% as.data.frame()

d<-SpatialPointsDataFrame(cbind(a$x,a$y),data=a,proj4string=CRS("+proj=longlat +datum=WGS84"))
ggplot(a,aes(x=x,y=y)) + borders(fill="grey") + coord_cartesian(xlim=as.numeric(bbox(d)[1,]),ylim=as.numeric(bbox(d)[2,])) + theme_bw()  + geom_path(size=1,aes(group=Animal,col=phistate)) + facet_wrap(~Animal) + mytheme()

