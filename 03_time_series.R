# This code is for site-level time series extraction
library(sp)
library(raster)
library(rgdal)
library(zyp)

ss <- 35

### 1. Load coordinations & make shpefile
crd <- read.csv('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/Selected_Sites_Lat_Lon.csv',
                head=F)

head(crd)
str(crd)
summary(crd)

# Phenology sites on Korea shapfile 
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))

shp_wr <- shapefile('/projectnb/modislc/users/mkmoon/KoreaPhenology/shpfile/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp')
nest <- c("Korea, Republic of")
shp_kp <- shp_wr[as.character(shp_wr$NAME) %in% nest,]
plot(shp_kp)

lat <- crd$V1
lon <- crd$V2
sites <- data.frame(lon,lat)
xy <- sites[,c(1,2)]
shp_po <- SpatialPointsDataFrame(coords=xy,data=sites,proj4string=crs(shp_kp))
plot(shp_po,cex=1.2,pch=19,col=2,lwd=1.5,add=T)

# setwd('/projectnb/modislc/users/mkmoon/KoreaPhenology/shpfile/')
# writeOGR(shp_po, ".", "phenology_sites", driver="ESRI Shapefile",overwrite=T)


### 2. Extract time series
lat <- crd$V1[ss]
lon <- crd$V2[ss]
sites <- data.frame(lon,lat)
xy <- sites[,c(1,2)]
shp_on <- SpatialPointsDataFrame(coords=xy,data=sites,proj4string=crs(shp_kp))
plot(shp_on,cex=3,col=4,lwd=5,add=T)

ltm <- raster('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/ko_ltm_Greenup_h28v05.tif')
plot(ltm)
shp_on_si <- spTransform(shp_on,crs(ltm))
plot(shp_on_si,add=T,cex=3,col=4,lwd=5)

# Load data
load('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/matrix_Greenup_h28v05')
p1 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
p2 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
p3 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
p4 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
p5 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
p6 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
p7 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
p8 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
p9 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
pt <- phe_value[c(p1,p2,p3,p4,p5,p6,p7,p8,p9),]

sub <- setValues(ltm,phe_value[,1])
sub <- crop(sub,extent(shp_on_si@bbox[1]-30000,shp_on_si@bbox[1]+30000,shp_on_si@bbox[2]-30000,shp_on_si@bbox[2]+30000))
plot(sub)
plot(shp_on_si,add=T,cex=5,col=4,lwd=5)

sub <- setValues(ltm,phe_value[,1])
sub <- crop(sub,extent(shp_on_si@bbox[1]-3000,shp_on_si@bbox[1]+3000,shp_on_si@bbox[2]-3000,shp_on_si@bbox[2]+3000))
plot(sub)
plot(shp_on_si,add=T,cex=5,col=4,lwd=5)


### 3. Plot
par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
pt_mean <- apply(pt,2,mean,na.rm=T)
pt_sd <- apply(pt,2,sd,na.rm=T)
plot(pt_mean,ylim=c(90,125),axe=F,ann=F,type='l',col=rgb(1,0,0),lwd=5)
polygon(c(seq(1,15),seq(15,1)),c(pt_mean-pt_sd,rev(pt_mean+pt_sd)),
        border=NA,col=rgb(1,0,0,0.1))
axis(1,at=seq(1,15,2),seq(2001,2015,2),cex.axis=1.3)
axis(2,at=seq(0,200,10),cex.axis=1.3)
box(lty=1)
mtext('Year',1,2.7,cex=1.5)
mtext('Day of year',2,2.5,cex=1.5)

mk <- MannKendall(pt_mean)
pix_mk <- as.numeric(mk$sl)

x <- as.integer(seq(1,15))
y <- pt_mean
z <- zyp.sen(y~x)
abline(z$coefficients[1],z$coefficients[2],lty=5,lwd=1.5)

text(2,98,paste('Mann-Kendall p-value: ',round(pix_mk,3),sep=''),cex=1.5,pos=4)
text(2,95,paste('Theil-sen slope: ',round(z$coefficients[2],3),sep=''),cex=1.5,pos=4)


### 4. Load & comapre with ground data
phe_gd <- read.csv('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/test_pheno.csv',
                   head=F)
pt_gd <- phe_gd[which(phe_gd[,1]==crd[ss,1]),5:19]

par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))
pt_mean <- apply(pt,2,mean,na.rm=T)
pt_sd <- apply(pt,2,sd,na.rm=T)
plot(pt_mean,ylim=c(min(pt_mean)-10,max(pt_mean)+10),axe=F,ann=F,type='l',col=rgb(1,0,0),lwd=5)
polygon(c(seq(1,15),seq(15,1)),c(pt_mean-pt_sd,rev(pt_mean+pt_sd)),
        border=NA,col=rgb(1,0,0,0.1))
axis(1,at=seq(1,15,2),seq(2001,2015,2),cex.axis=1.3)
axis(2,at=seq(0,200,10),cex.axis=1.3)
box(lty=1)
mtext('Year',1,2.7,cex=1.5)
mtext('Day of year',2,2.5,cex=1.5)
points(1:15,pt_gd,type='o')

plot(pt_mean,pt_gd,xlim=c(90,130),ylim=c(90,130))
abline(0,1)
arrows(pt_mean+pt_sd,as.numeric(pt_gd),pt_mean-pt_sd,as.numeric(pt_gd),angle=90,code=3,length=0)

# Color
mycol <- colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# set.seed(12344)
mycol <- sample(mycol,93)
plot(1:93,col=mycol,pch=19)

# 
par(mfrow=c(3,3),oma=c(1,1,1,1),mar=c(1,1,1,1),mgp=c(2,0.5,0))
for(i in 1:90){
  # Ground data
  pt_gd <- phe_gd[i,5:19]
  
  # Extract MODIS
  lat <- phe_gd[i,1]
  lon <- phe_gd[i,2]
  sites <- data.frame(lon,lat)
  xy <- sites[,c(1,2)]
  shp_on <- SpatialPointsDataFrame(coords=xy,data=sites,proj4string=crs(shp_kp))
  shp_on_si <- spTransform(shp_on,crs(ltm))
  
  p1 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
  p2 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
  p3 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)-1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
  p4 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
  p5 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
  p6 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+0)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
  p7 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)-1
  p8 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+0
  p9 <- (ceiling((ltm@extent@ymax - shp_on_si@bbox[2])/463.3127)+1)*2400 + ceiling((shp_on_si@bbox[1] - ltm@extent@xmin)/463.3127)+1
  pt <- phe_value[c(p1,p2,p3,p4,p5,p6,p7,p8,p9),]
  
  pt_mean <- apply(pt,2,mean,na.rm=T)
  pt_sd <- apply(pt,2,sd,na.rm=T)
    
  plot(pt_mean,pt_gd,xlim=c(0,180),ylim=c(0,180),pch=19,cex=2)
  abline(0,1)
  arrows(pt_mean+pt_sd,as.numeric(pt_gd),pt_mean-pt_sd,as.numeric(pt_gd),angle=90,code=3,length=0)  

#   if(i==1){
#     plot(pt_mean,pt_gd,xlim=c(0,180),ylim=c(0,180),pch=19,col=mycol[i])
#     abline(0,1)
# #     arrows(pt_mean+pt_sd,as.numeric(pt_gd),pt_mean-pt_sd,as.numeric(pt_gd),angle=90,code=3,length=0)  
#   }else{
#     points(pt_mean,pt_gd,pch=19,col=mycol[i])
# #     arrows(pt_mean+pt_sd,as.numeric(pt_gd),pt_mean-pt_sd,as.numeric(pt_gd),angle=90,code=3,length=0)
#   }

#   for(i in 1:15){
#     if(is.na(pt_mean[i])) pt_gd[i] <- NA
#   }
#   for(i in 1:15){
#     if(is.na(pt_gd[i])) pt_mean[i] <- NA
#   }
# 
#   pt_mm <- mean(pt_mean,na.rm=T)
#   pt_gg <- mean(as.numeric(pt_gd),na.rm=T)
# 
# 
#   if(i==1){
#     plot(pt_mm,pt_gg,xlim=c(0,180),ylim=c(0,180))
#   }else{
#     points(pt_mm,pt_gg)
#   }
}


### Test github
## test 1


