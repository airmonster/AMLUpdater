# AMLUpdater  

## A program used for update wrf land use data with high resolution land use data, by Airmonster team from Chengdu Academy of Environmental Sciences.  

The code was written in R. geo_em.d* data was processed off-line which means you should copy your data from WPS dir and then process with this program. Tests were carried out on Window with R 3.5.2.

This program was inspired by this link: http://bbs.06climate.com/forum.php?mod=viewthread&tid=2435  

The default land use data in this program was FROM-GLC shared by Tsinghua University: http://data.ess.tsinghua.edu.cn/fromglc2015_v1.html  

## HOW
As you know, land use was an important part of WRF static data, and the default land use data from the model were outdated for most cities in China because of the rapid urbanization. It is quite difficult for beginners to update the land use data which acquires technics including programming and gis, so we provide an easy way to start with.
FROM-GLC 2015 v1 data with remapping according to modis land use scheme in WRF was used in this program to calculate the ratio of different land use types, and updates were made to geo_em data.
* Dependency of AMLUpdater are a list of R packages including:  
*Raster* : raster process, including clip and merge.  
*RNetCDF* : RW for netcdf files.  
*doParallel* ：Parallel process for grids to improve effiency.  
* Before running the program, some parameters should be modified:  
*DATDIR* :Local path for FROM-GLC data.  
*GEODIR* :Upper folder for geo_em files.  
*case* :Under GEODIR, and inside which are geo_em.d0?.nc  
*NSL* :LAT range in format of (start-lat, end-lat, step), notice step should be 10 for FROM-GLC 2015v1.  
*EWL* : LON range in format of (start-lon, end-lon, step), notice step should be 10 for FROM-GLC 2015v1.  
*DLD* : T/F, whether data needs to be downloaded, T for the first time.  
*MAXDM* : max number of the nested domain in WRF.  

## TIPS
Inside the code, a special parameter named ubonly should be set to T if you want to update urban and waterbody only, sometimes it is a better choice for you.
In our tests, after updating the land use data, correlation of temperature were slightly improved for all 13 observation points in Chengdu, and correlation of the wind speed of 5 observation points were improved. 

## Author
lcw@cdaes.cn
