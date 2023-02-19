nc_convert <-
function(nc_data_location){
  options(digits=10)
  requireNamespace("ncdf4")

  kelp_biomass_file=ncdf4::nc_open(nc_data_location)
  
  lat<-ncdf4::ncvar_get(kelp_biomass_file,attributes(kelp_biomass_file$var)$names[1])
  long<-ncdf4::ncvar_get(kelp_biomass_file,attributes(kelp_biomass_file$var)$names[2])
  year<-ncdf4::ncvar_get(kelp_biomass_file,attributes(kelp_biomass_file$var)$names[3]) 
  quarter<-ncdf4::ncvar_get(kelp_biomass_file,attributes(kelp_biomass_file$var)$names[4])
  biomass<-ncdf4::ncvar_get(kelp_biomass_file,attributes(kelp_biomass_file$var)$names[5])
  
  #binding it together 
  biomass_data<-data.frame(Long=long,Lat=lat,Biomass=biomass)
  colnames(biomass_data)[-c(1,2)]=paste("Biomass",year,quarter,sep=".")
  
  return(biomass_data)
  
}
