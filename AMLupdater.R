# Name: AMLupdater
# Version: 1.0.200330
# Description: AMLUpdater, a program used for update wrf land use data with high resolution land use data (FROM-GLC), by Airmonster team from Chengdu, the code was written in R. geo_em.d* data was processed off-line which means you should copy your data from WPS dir and then process with this program. Tests were carried out on Window with R 3.5.2.
# Author: lcw@cdaes.cn
# URI: https://github.com/airmonster/AMLupdater

library(raster)
library(RNetCDF)
library(doParallel)

# FROM-GLC version2 (2015_v1) download web URL, shared by Tsinghua University.
GLCURL="http://data.ess.tsinghua.edu.cn/data/temp/Fromglc2015tif/"
# FROM-GLC local data storage path.
DATDIR="\\path\\to\\GLC2015\\"
# Local WRF geo_em path.
GEODIR="\\path\\to\\GEOEM\\"
case="CDTEST01"

NSL=seq(0,80,10) #lat
EWL=seq(20,180,10) #lon
MAXDM=2
DLD=F

if(DLD==T){
  try(
    {
      for(xi in 1:length(NSL))
        for(yi in 1:length(EWL)){
          try({
            download.file(url=paste0(GLCURL,EWL[yi],"E","_",NSL[xi],"N",".tif"),
                          mode = "wb",
                          destfile=paste0(DATDIR,EWL[yi],"E","_",NSL[xi],"N",".tif"),
                          quiet=F)
          })
          
        } 
    }
  )
}

getclu=function(xt,vis=F){
  crp=NA
  for(xi in 1:length(NSL))
    for(yi in 1:length(EWL)){
      rasterfile=paste0(DATDIR,EWL[yi],"E","_",NSL[xi],"N",".tif")
      if(file.exists(rasterfile)){
        rasterdf=raster(rasterfile,as.is=T)
        #print(paste0("Processing:",rasterfile))
        tr=try(crop(rasterdf,xt),silent = T)
        if('try-error' %in% class(tr)){
          #print("Not overlap.")
        }else{
          if(length(crp)==1){
            crp=crop(rasterdf,xt)
          }else{
            crp=merge(crp,crop(rasterdf,xt))
          }
        }
      }
    }
  #if(vis==T) plot(crp,zlim=c(0,120))
  return(as.vector(crp))
}

get_points=function(i,j,lat_v,lon_u){
  lat1=lat_v[i,j]
  lat2=lat_v[i,j+1]
  lon1=lon_u[i,j]
  lon2=lon_u[i+1,j]
  lims=c(lat1,lat2,lon1,lon2)
  tp=extent(lims[3],lims[4],lims[1],lims[2])
  ludf=getclu(tp)
  return(ludf)
}

getpoly=function(i,j,lat_v,lon_u){
  gridpoly=data.frame(x=c(lon_u[i,j],lon_u[i+1,j],lon_u[i+1,j],lon_u[i,j],lon_u[i,j]),
                      y=c(lat_v[i,j],lat_v[i,j],lat_v[i,j+1],lat_v[i,j+1],lat_v[i,j])
  )
  return(gridpoly)
}

process.lu.par=function(i,jmx,lutype){
  nludf=data.frame()
  ubonly=F
  
  ubid=13
  gsid=14
  wtid=17
  ftid=2
  for(j in 1:jmx){
    tlu=LANDUSEF[i,j,]
    maxu=which(tlu==max(tlu))
    minu=min(tlu)
    clu=get_points(i,j,lat_v,lon_u)
    #clu=subset(clu,clu!=0)
    if(length(clu)>0){
      ntlu=tlu
      
      if(ubonly==T){
        ttlu=tlu
        ntlu[ubid]=round(length(subset(clu,clu==80))/length(clu),6) #
        ntlu[wtid]=length(subset(clu,clu==60))/length(clu) #
        
        
        if(LUINDEX[i,j]==12 | LUINDEX[i,j]==14 | LUINDEX[i,j]==10){ #MODIS ONLY
          ntlu[1]=length(subset(clu,clu==23))/length(clu) #
          ntlu[2]=length(subset(clu,clu==21))/length(clu) #
          #ntlu[3]=length(subset(clu,clu==3))/length(clu) #
          ntlu[4]=length(subset(clu,clu==22))/length(clu) #
          ntlu[5]=length(subset(clu,clu==25 | clu==26))/length(clu) #
        }
        
        
        rnlu=1-ntlu[ubid]-ntlu[wtid]
        rolu=1-tlu[ubid]-tlu[wtid]
        if(rolu==0){ 
          ntlu=ttlu*0
          ntlu[ubid]=length(subset(clu,clu==80))/length(clu) #
          ntlu[wtid]=length(subset(clu,clu==60))/length(clu) #
          if(rnlu>0) ntlu[gsid]=rnlu
        }else{
          ntlu=ttlu*rnlu/rolu
          ntlu[ubid]=length(subset(clu,clu==80))/length(clu) #
          ntlu[wtid]=length(subset(clu,clu==60))/length(clu) #
        }
      }else{ #All category
        ntlu=ntlu*0
        
        ntlu[1]=length(subset(clu,clu==23))/length(clu) #
        ntlu[2]=length(subset(clu,clu==21))/length(clu) #
        #ntlu[3]=length(subset(clu,clu==3))/length(clu) #
        ntlu[4]=length(subset(clu,clu==22))/length(clu) #
        ntlu[5]=length(subset(clu,clu==25 | clu==26))/length(clu) #
        ntlu[6]=length(subset(clu,clu==41))/length(clu) #
        ntlu[7]=length(subset(clu,clu==42))/length(clu) #
        #ntlu[8]=length(subset(clu,clu==5))/length(clu) #
        #ntlu[9]=length(subset(clu,clu==2))/length(clu) #
        ntlu[10]=length(subset(clu,clu==31 | clu==32 | clu==33))/length(clu) #
        ntlu[11]=length(subset(clu,clu==51 | clu==52 | clu==53))/length(clu) #
        ntlu[12]=length(subset(clu,clu==11 | clu==12 | clu==13| clu==15))/length(clu) #
        ntlu[13]=length(subset(clu,clu==80))/length(clu) #
        #ntlu[14]=length(subset(clu,clu==9))/length(clu) #
        ntlu[15]=length(subset(clu,clu==101 | clu==102))/length(clu) #
        ntlu[16]=length(subset(clu,clu==90 | clu==92))/length(clu) #
        ntlu[17]=length(subset(clu,clu==60))/length(clu) #
        ntlu[18]=length(subset(clu,clu==71))/length(clu) #
        ntlu[19]=length(subset(clu,clu==72))/length(clu) #
        #ntlu[20]=length(subset(clu,clu==9))/length(clu) #
        
      }
    }
    nludf=rbind(nludf,data.frame(I=i,J=j,LUI=which(ntlu==max(ntlu))[1],t(ntlu)))
  }
  
  return(nludf)
}

process.par=function(dm,case){
  geofile=paste0(GEODIR,case,"/geo_em.d0",dm,".nc")
  print(paste0("GEOFILE:",geofile))
  dmnc=open.nc(geofile,write=T)
  #NOTICE, the grid file should be modis
  lutype=att.get.nc(dmnc,"NC_GLOBAL","MMINLU")
  wrfver=att.get.nc(dmnc,"NC_GLOBAL","TITLE")
  print(paste0("Processing ",lutype," data on domain ",dm,"."))
  
  update=F
  FIX_ONLY_LARGER=F
  LANDUSEF=var.get.nc(dmnc,"LANDUSEF")
  LUINDEX=var.get.nc(dmnc,"LU_INDEX")
  LANDMASK=var.get.nc(dmnc,"LANDMASK")
  
  lat_m=var.get.nc(dmnc,"XLAT_M")
  lat_u=var.get.nc(dmnc,"XLAT_U")
  lat_v=var.get.nc(dmnc,"XLAT_V")
  
  lon_m=var.get.nc(dmnc,"XLONG_M")
  lon_u=var.get.nc(dmnc,"XLONG_U")
  lon_v=var.get.nc(dmnc,"XLONG_V")
  
  if(wrfver!="OUTPUT FROM GEOGRID V3.5"){
    lai12m=var.get.nc(dmnc,"LAI12M")
  }
  
  alb12m=var.get.nc(dmnc,"ALBEDO12M")
  gre12m=var.get.nc(dmnc,"GREENFRAC")
  
  hgt=var.get.nc(dmnc,"HGT_M")
  #for the moment only fix the urban , not focused on all
  xm=range(lon_m)
  ym=range(lat_m)
  
  tcore=detectCores(logical=T)
  cl=makeCluster(tcore)
  registerDoParallel(cl) #
  print(paste0("Start Parallel Process on a ",tcore," cores machine."))
  jmx=dim(LANDUSEF)[2]
  newludf = foreach(i=1:dim(LANDUSEF)[1],.packages="raster",.verbose=T,.combine = rbind,.export=c("DATDIR","GEODIR","EWL","get_points","getclu","jmx","LANDUSEF","lat_v","lon_u","LUINDEX","lutype","NSL","process.lu.par")) %dopar% {
    process.lu.par(i,jmx,lutype)
  }
  stopCluster(cl)
  
  pb=txtProgressBar(max=dim(LANDUSEF)[1],style=3)
  for(i in 1:dim(LANDUSEF)[1]){
    for(j in 1:dim(LANDUSEF)[2]){
      templu=subset(newludf,I==i & J==j)
      LANDUSEF[i,j,]=as.vector(t(templu[,4:length(templu[1,])]))
      LUINDEX[i,j]=templu$LUI
    }
    image(LUINDEX,col=rainbow(50))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  LANDUSEF_NORM=array(0,dim=c(dim(LANDUSEF),1))
  LANDUSEF_NORM[,,,1]=LANDUSEF
  LI_NORM=array(0,dim=c(dim(LANDUSEF)[1:2],1))
  LI_NORM[,,1]=LUINDEX
  LM_NORM=array(0,dim=c(dim(LANDUSEF)[1:2],1))
  
  LANDMASK[LUINDEX==17]=0
  LM_NORM[,,1]=LANDMASK
  
  GREENFRAC_NORM=array(0,dim=c(dim(gre12m),1))
  for(mi in 1:12){
    FGREEN=(1-LANDUSEF[,,13])*max(gre12m[,,mi][LUINDEX==13])
    gre12m[,,mi][LUINDEX==13]=FGREEN[LUINDEX==13]
  }
  GREENFRAC_NORM[,,,1]=gre12m
  
  var.put.nc(dmnc,variable = "LANDUSEF",data=LANDUSEF_NORM,c(1,1,1,1))
  var.put.nc(dmnc,variable = "LU_INDEX",data=LI_NORM,c(1,1,1))
  var.put.nc(dmnc,variable = "LANDMASK",data=LM_NORM,c(1,1,1))
  var.put.nc(dmnc,variable = "GREENFRAC",data=GREENFRAC_NORM,c(1,1,1,1))
  close.nc(dmnc)
}

###########################################

for(dm in 1:MAXDM) {
  system.time({
    process.par(dm,case)
  })
}
print("Done.")
