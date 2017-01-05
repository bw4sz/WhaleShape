#Download sea ice data

#Readme: https://nsidc.org/data/seaice_index/more-about-monthly.html
library(stringr)
library(rgdal)
library(maptools)
library(raster)
library(ggmap)
library(scales)
library(broom)
library(sp)
library(rgeos)
library(tidyr)

#Build the url
pols<-list()

Months<-month.abb[1:7]
MonthNums=c("01","02","03","04","05","06","07")
Years=c(2012,2013,2015,2016)

#construct urls
urls<-list()
for(y in 1:length(Years)){
  yearurls<-list()
  for(x in 1:length(Months)){
    yearurls[[x]]<-paste('ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/shapefiles/',Months[x],'/shp_extent/extent_S_',paste(Years[y],MonthNums[x],sep=""),'_polygon_v2.zip',sep="")
  }
  urls[[y]]<-yearurls
}
urls<-unlist(urls)


for(x in 1:length(urls)){

fl<-paste("InputData/SeaIce/",str_match(urls[[x]],"extent_S_(\\d+)")[,2],sep="")
download.file(urls[[x]],destfile=paste(fl,".zip",sep=""))

#unzip
unzip(zipfile=paste(fl,".zip",sep=""),exdir="./InputData/SeaIce")
}

#cut polygon
#cutp<-readShapePoly("InputData/CutPolygon.shp")
#reproject into stereographic

# #read
# shp<-list.files("InputData/SeaIce",pattern=".shp",full.names=T)
# pols<-list()
# for(a in 1:length(shp)){
# 
#   rp<-readShapeLines(shp[[a]])
#   #define projection
#   stere <- "+proj=stere +lat_0=-90 +lat_ts=70 +lon_0=0 +datum=WGS84 +units=m"
#   proj4string(rp)<-stere
# 
#   #0 width buffer
#   #rp <- gBuffer(rp, byid=TRUE, width=0)
#   #crop
#   e<-extent(-3455516,-757334.6,66739.71,2200756)*1.1
#   rmask<-raster(rp,ext=e)
#   rmask<-disaggregate(rmask,10)
#   rasP<-rasterize(rp,rmask)
#   #presence of ice
#   rasP <- projectRaster(rasP,crs=CRS("+proj=longlat +ellps=WGS84"))
#   rasP<-rasP>0
#   pols[[a]]<-rasP
#   
#   #name the layers
#   #get naming structure
#   s<-str_match(shp[[a]],"S_(\\d+)_")[,2]
#   yr<-as.numeric(substring(s,0,4))
#   mn<-as.numeric(substring(s,5,6))
#   names(pols[[a]])<-paste(month.abb[mn],yr,sep="_")
# }
# spols<-stack(pols)
# #writeindividaul raster
# writeRaster(spols,"InputData/MonthlyIceRaster",overwrite=T)

#as individual polygons
#read
shp<-list.files("InputData/SeaIce",pattern="polygon_v2.shp",full.names=T)
polygonlist<-list()
polygons<-list()
for(a in 1:length(shp)){
  
  rp<-readShapePoly(shp[[a]])
  #define projection
  stere <- "+proj=stere +lat_0=-90 +lat_ts=70 +lon_0=0 +datum=WGS84 +units=m"
  proj4string(rp)<-stere
  
  #crop
  rp<-gBuffer(rp, byid=TRUE, width=0)
  e<-extent(-3455516,-757334.6,66739.71,2200756)*1.1
  crp<-crop(rp,e)
  plot(tcrp<-spTransform(crp,CRS("+proj=longlat +ellps=WGS84")))

  #name the layers
  #get naming structure
  s<-str_match(shp[[a]],"S_(\\d+)_")[,2]
  yr<-as.numeric(substring(s,0,4))
  mn<-as.numeric(substring(s,5,6))
  plot(tcrp)
  polygons[[a]]<-tcrp
  polygonlist[[a]]<-fortify(tcrp,region="INDEX")
  names(polygonlist)[[a]]<-paste(month.name[mn],yr,sep="_")
  names(polygons)[[a]]<-paste(month.name[mn],yr,sep="_")
}

#write each
for(x in 1:length(polygonlist)){
  write.csv(polygonlist[[x]],paste("C:/Users/Ben/Documents/Whales/InputData/SeaIce/",names(polygonlist)[[x]],"PolygonLatLong.csv",sep=""))
}

#as one list
mlist<-melt(polygonlist,id.vars=colnames(polygonlist[[1]]))

#split months and years
mlist<-mlist %>% separate(L1,c("Month","Year"),"_")
ggplot(ice[ice$Month=="April",],aes(x=long,y=lat,col=Year,group=paste(Year,group))) + facet_wrap(~Month) + geom_polygon() + theme_bw() + borders(fill="black")  + coord_cartesian(ylim = c(-69,-62),xlim=c(-55,-72)) + theme_bw()
write.csv(mlist,"InputData/SeaIce_AllYears_polygons.csv")
