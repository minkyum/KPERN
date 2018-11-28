## This code is about raster manipulation.
# Loading required libraries 
library(sp)
library(raster)
library(RColorBrewer)

### 1. Load shapefile for Korea
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(4,4,4,4),mgp=c(2.5,1,0))

shp_wr <- shapefile('/projectnb/modislc/users/mkmoon/KoreaPhenology/shpfile/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp')
plot(shp_wr)

shp_wr$NAME

nest <- c("Korea, Democratic People's Republic of","Korea, Republic of")
shp_kp <- shp_wr[as.character(shp_wr$NAME) %in% nest,]
plot(shp_kp)

### 2. Load phenology data
path <- paste('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/')
file <- list.files(path=path,pattern=glob2rx('*tif'),full.names=T)

ltm <- vector('list',3)
tre <- vector('list',3)
for(i in 1:3){
  ltm[[i]] <- raster(file[i])
  tre[[i]] <- raster(file[i+3])  
}
ltm <- merge(ltm[[1]],ltm[[2]],ltm[[3]])
tre <- merge(tre[[1]],tre[[2]],tre[[3]])

plot(ltm,main='Long-term mean')
plot(tre,main='15-years trend')
hist(ltm,main='Long-term mean')
hist(tre,xlim=c(-5,5),breaks=seq(-500,500,0.2),main='15-years trend')

### 3. Reprojection & Masking
plot(shp_kp)
shp_kp_sinu <- spTransform(shp_kp,crs(ltm))
plot(shp_kp_sinu)

plot(ltm,main='Long-term mean')
plot(shp_kp_sinu,add=T)
plot(tre,main='Trend')
plot(shp_kp_sinu,add=T)

# reprojection
pr3 <- projectExtent(ltm,crs(shp_kp))
res(pr3) <- 0.05
ltm_ll <- projectRaster(ltm,pr3)

pr3 <- projectExtent(tre,crs(shp_kp))
res(pr3) <- 0.05
tre_ll <- projectRaster(tre,pr3)

plot(ltm_ll)
plot(shp_kp,add=T)
plot(tre_ll)
plot(shp_kp,add=T)

# crop
ltm_ll_c <- crop(ltm_ll,extent(shp_kp))
tre_ll_c <- crop(tre_ll,extent(shp_kp))

plot(ltm_ll_c)
plot(shp_kp,add=T)
plot(tre_ll_c)
plot(shp_kp,add=T)

# mask
ltm_ll_m <- mask(ltm_ll_c,shp_kp)
tre_ll_m <- mask(tre_ll_c,shp_kp)

plot(ltm_ll_m)
plot(shp_kp,add=T)
plot(tre_ll_m)
plot(shp_kp,add=T)

### 4. Plot
setwd('/projectnb/modislc/users/mkmoon/KoreaPhenology/figure/')
pdf(file='ko_ltm_tre',width=12,height=7)  

par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(1,1,1,1),mgp=c(1,1,0),bty='n')

lb <- quantile(values(ltm_ll_m),0.02,na.rm=T)
ub <- quantile(values(ltm_ll_m),0.98,na.rm=T)
values(ltm_ll_m)[values(ltm_ll_m)<lb] <- lb
values(ltm_ll_m)[values(ltm_ll_m)>ub] <- ub  

mycol <- brewer.pal(11,'Spectral')
mycol <- colorRampPalette(mycol)(100)
plot(ltm_ll_m,zlim=c(lb,ub),legend=F,axes=F,col=mycol)
plot(shp_kp,add=T)
plot(ltm_ll_m,
     legend.only=T,
     col=mycol,
     zlim=c(lb,ub),
     legend.width=1,
     legend.shrink=1,
     horiz=F,
     smallplot=c(0.1,0.13,0.1,0.4),
     axis.args=list(at=c(100,110,120,130,140,150,160),
                    cex.axis=1),
     legend.args=list(text='Day of year',side=4,font=1,line=2.7,cex=1.3)) 
mtext('Long-term mean Green-up',3,-2,outer=T,adj=0.03,cex=1.3)

values(tre_ll_m)[values(tre_ll_m)< -1.2] <- -1.2
values(tre_ll_m)[values(tre_ll_m)> 1.2] <- 1.2

mycol <- brewer.pal(11,'RdBu')
mycol <- colorRampPalette(mycol)(100)
plot(tre_ll_m,zlim=c(-1.5,1.5),legend=F,axes=F,col=mycol)
plot(shp_kp,add=T)
plot(tre_ll_m,
     legend.only=T,
     col=mycol,
     zlim=c(-1.2,1.2),
     legend.width=1,
     legend.shrink=1,
     horiz=F,
     main='dd',
     smallplot=c(0.1,0.13,0.1,0.4),
     axis.args=list(at=c(-1.2,-0.8,-0.4,0,0.4,0.8,1.2),
                    cex.axis=1),
     legend.args=list(text='Trend (days/year)',side=4,font=1,line=2.7,cex=1.3)) 
mtext('15-years Green-up trend',3,-2,outer=T,adj=0.68,cex=1.3)

dev.off()