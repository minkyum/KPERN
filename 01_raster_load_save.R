## This code is for loading and saving raster files.
# Loading required libraries 
library(sp)
library(raster)
library(Kendall)
library(zyp)

# From bash code
args <- commandArgs()
print(args)
tt <- as.numeric(args[3])

if(tt==1){
  tile <- 'h27v04'
}else if(tt==2){
  tile <- 'h27v05'
}else{
  tile <- 'h28v05'
}
print(tile)


### 1. Set variables
phe <- c('Greenup','Maturity','Senescence','Dormancy','MidGreenup','MidGreendown')  
vari <- 1

### 2. Load land cover type for extent
if(tt==1){
  ext <- raster('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp/igbp.h27v04.2015.bip')
}else if(tt==2){
  ext <- raster('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp/igbp.h27v05.2015.bip')
}else{
  ext <- raster('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp/igbp.h28v05.2015.bip')
}

### 3. Load land surface phenology as matrix
phe_value <- matrix(NA,(2400*2400),15)
for(year in 2001:2015){
  ncol <- 2400
  nrow <- 2400
  nbands <- 2
  cnt <- ncol*nrow*nbands
      
  path <- paste('/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/',phe[vari],sep='')
  sstr <- paste(phe[vari],'*',tile,'*',year,sep='')
  file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  
  data <- readBin(file,what="integer",n=cnt,size=2,endian="little")
  data[data>30000] <- NA
  data1 <- data - as.numeric(as.Date(paste(year-1,'-12-31',sep='')))
  data2 <- array(data1,c(nbands, ncol, nrow))
  data2 <- aperm(data2, c(3,2,1)) #for transposing
  aa <- brick(data2,
              xmn=ext@extent@xmin,
              xmx=ext@extent@xmax,
              ymn=ext@extent@ymin,
              ymx=ext@extent@ymax,
              crs=ext@crs)
  phe_value[,(year-2000)] <- getValues(aa[[1]])
  
  print(year)
}

### 4. Long-term mean & Trend
ltm <- apply(phe_value,1,mean,na.rm=T)
ltm <- setValues(ext,ltm)

# pix_mk <- matrix(NA,nrow(phe_value),1)
pix_th <- matrix(NA,nrow(phe_value),1)
for(i in 1:nrow(phe_value)){
#   mk <- MannKendall(phe_value[i,])
#   pix_mk[i,1] <- as.numeric(mk$sl)
  
  x <- as.integer(seq(1,15))
  y <- phe_value[i,]
  z <- zyp.sen(y~x)
  pix_th[i,1] <- z$coefficients[2]   

  if(i%%100000==0) print(i)
}
tre <- setValues(ext,pix_th)

### 5. Save as tif file
setwd('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/')
writeRaster(ltm,paste('ko_ltm_',phe[vari],'_',tile,sep=''),format='GTiff',overwrite=T)
writeRaster(tre,paste('ko_tre_',phe[vari],'_',tile,sep=''),format='GTiff',overwrite=T)
save(phe_value,file=paste('matrix_',phe[vari],'_',tile,sep=''),overwrite=T)
