f<-mxy %>% filter(Animal==9,Track %in% 12)

d<-SpatialPointsDataFrame(cbind(f$x,f$y),data=f,proj4string=CRS("+proj=longlat +datum=WGS84"))
ggplot(f,aes(x=x,y=y,col=phistate)) + borders(fill="black") + coord_cartesian(xlim=as.numeric(bbox(d)[1,]),ylim=as.numeric(bbox(d)[2,])) + theme_bw() + facet_wrap(~Track) + geom_path(aes(group=Track),size=1)

