decline_finder <-
function(data,baseline_threshold=0.1,scarce_cutoff=0.4,present_window=16,hist_period=100,window_lag=0,
                           lat_min=27.01,lat_max=37.5,lon_min=-123.5,lon_max=-114,table_name=NULL) {

  scarce_cutoff=(1-scarce_cutoff)
  #loading in the dataset
  if(window_lag>0){kelp_biomass_data=data[,-c((ncol(data)+1-window_lag):ncol(data))]}else{
    kelp_biomass_data=data}

  kelp_biomass_data=kelp_biomass_data[(kelp_biomass_data$Long>=lon_min  & kelp_biomass_data$Long<=lon_max
                                       &  kelp_biomass_data$Lat>=lat_min  &  kelp_biomass_data$Lat<=lat_max), ]

  #kelp_biomass_data=data[,-c((ncol(data)-window_lag):ncol(data))]


  ####MEAN CALCULATIONS#####
  mean_data=kelp_biomass_data[,-c(5:ncol(kelp_biomass_data))]#we're just using the structure of kelp_biomass_data to create a new dataframe
  colnames(mean_data)[c(3,4)]=c("mean","zero count")
  mean_data[,c(3,4)]=NA

  #Calculating long-term means and counting the number of zeroes for each pixel
  start_column=ncol(kelp_biomass_data)-(hist_period-1)#find the data corresponding to the number of years back from present day to mean period
  mean_data$mean=apply(kelp_biomass_data[,start_column:ncol(kelp_biomass_data)], 1,function(x) mean(x,na.rm = T))
  mean_data$`zero count`=apply(kelp_biomass_data[,start_column:ncol(kelp_biomass_data)], 1, function(x) length(which(x==0)))###I can probably move this to after
  #the means are calculated, after all if the means are zero then the whole vector is zeroes or NA values

  #removing rows with zero mean, finding their latitude and longitude coordinates
  zero_finder=which(mean_data$mean!=0)
  filtered_kelp_biomass_data=kelp_biomass_data[zero_finder,]
  kelp_biomass_locations=filtered_kelp_biomass_data[,1:2]

  ####getting a dataset that has only the most recent years in the data
  number_of_windows=(ncol(filtered_kelp_biomass_data)-present_window)+seq(1:present_window)#setting up how many quarters we want to look back for declines

  ####THRESHOLD CALCULATION AND FILTERTING####
  #setting a portion of the mean as the threshold for
  mean_data$thresholds=baseline_threshold*mean_data$mean

  #creating an array of T or F values indicating if the values in a given cell were above or below the mean
  biomass_decline_logical=apply(filtered_kelp_biomass_data[3:ncol(filtered_kelp_biomass_data)],2,function(x) x>mean_data$thresholds[zero_finder])
  biomass_decline_logical=cbind(filtered_kelp_biomass_data[,1:2],biomass_decline_logical)

  #finding rows in the declined window quarters where there have been decreases the whole time
  number_of_declined_windows=ncol(biomass_decline_logical)-present_window+1:present_window
  declined_rows=apply(biomass_decline_logical[,number_of_declined_windows], 1, function(x) all(x==F,na.rm = T))

  #developing a filter to choose rows that only have a certain number of zeroes, with the default being at least 40% being nonzero
  mean_data$zero_filter=mean_data$`zero count`/(ncol(kelp_biomass_data)-start_column)<scarce_cutoff
  combo_logical=declined_rows & mean_data$zero_filter[zero_finder]
  filtered_decline_locations=kelp_biomass_locations[combo_logical,]

  if(nrow(filtered_decline_locations)==0){

    print.noquote(paste("No biomass pixels have less than",scarce_cutoff,"proportion of zero values for this time period" ))
    stop(nrow(filtered_decline_locations)==0)

  }

  ####RASTER WORK AND COUNTING DECLINES REGIONALLY#####
  requireNamespace("sp")
  requireNamespace("raster")###I really should look into if I will need to replace this....
  #####IMPORTANT NOTE: The raster package may need be replaced with Terra in the future, though I will not be the one bothering to update all of the code
  #####READ THIS ARTICLE FOR MORE INFORMATION: https://www.r-bloggers.com/2022/04/r-spatial-evolution-retirement-of-rgdal-rgeos-and-maptools/

  #creating our points
  decline_spatialpoints=sp::SpatialPoints(filtered_decline_locations)
  crs(decline_spatialpoints)=sp::CRS("+proj=longlat +datum=WGS84")

  #creating our raster to cover this area for counts and finding regional patterns.The resolution and cell locations of this raster line up exactly
  #with the raster data for the long term sea surface temperature measurements
  empty_raster=raster::raster(nrow=59,ncol=45,xmn=-125,xmx=-113.75,ymn=27.5,ymx=42.25,crs="+proj=longlat +datum=WGS84")


  #choosing an extent to crop our raster
  crop_extent=raster::extent(lon_min,lon_max,lat_min,lat_max)
  empty_raster=raster::crop(empty_raster,crop_extent)

  #counting the number of non-scarce and declined locations in each cell of empty_raster
  count_raster=empty_raster;count_raster[]=0
  count=table(cellFromXY(empty_raster,decline_spatialpoints))
  count_raster[as.numeric(names(count))]=count

  #now counting the number of solely non-scarce locations in each cell
  non_scarce_locations=kelp_biomass_locations[mean_data$zero_filter[zero_finder],1:2]
  non_scarce_locations=sp::SpatialPoints(non_scarce_locations)

  count_raster_non_scarce=empty_raster;count_raster_non_scarce[]=0
  count_non_scarce=table(cellFromXY(count_raster_non_scarce,non_scarce_locations))
  count_raster_non_scarce[as.numeric(names(count_non_scarce))]=count_non_scarce

  ratio_raster=count_raster/count_raster_non_scarce

  if(!is.null(table_name)){

    raster_table=na.omit(values(ratio_raster))
    cell_filter=!(as.vector(seq(1:length(values(ratio_raster))) %in% attr(raster_table,"na.action")))
    raster_table=as.data.frame(as.vector(raster_table));colnames(raster_table)="PPD"

    Pixels=values(count_raster_non_scarce)[cell_filter]
    raster_table=cbind(raster_table,Pixels)

    cell_centers=sp::coordinates(empty_raster)[cell_filter,];colnames(cell_centers)=c("x","y")
    cell_corners=matrix(nrow=nrow(cell_centers),ncol=4);cell_corners[,c(1,3)]=-0.125;cell_corners[,c(2,4)]=0.125
    cell_corners[,c(1,2)]=cell_corners[,c(1,2)] + cell_centers[,1]; cell_corners[,c(3,4)]=cell_corners[,c(3,4)] + cell_centers[,2]
    colnames(cell_corners)=c("long_west","long_east","lat_south","lat_north")

    #finding means of the non_scarce_locations
    non_scarce_coords=as.data.frame(non_scarce_locations)
    non_scarce_filter= mean_data$Long %in% mean_data$Long & mean_data$Lat %in% non_scarce_coords$Lat
    non_scarce_mean_data=mean_data[non_scarce_filter,c(1:3)]

    #going cell by cell and taking the aggregate mean for each cell
    AVG.biomass=vector(length = nrow(cell_centers))

    for(i in 1:length(AVG.biomass)){

      temp_filter= (cell_corners[i,1] < non_scarce_mean_data$Long) & (non_scarce_mean_data$Long < cell_corners[i,2]) & (cell_corners[i,3]< non_scarce_mean_data$Lat) & (non_scarce_mean_data$Lat < cell_corners[i,4])

      temp_biomasses= non_scarce_mean_data$mean[temp_filter]
      cell_mean=mean(temp_biomasses)

      AVG.biomass[i]=cell_mean


    }

    raster_table=cbind(raster_table,AVG.biomass)

    kelp.area_km.sq=0.03*0.03*Pixels; raster_table=cbind(raster_table, kelp.area_km.sq)

    raster_table=cbind(raster_table,cell_corners,rep(NA,nrow(raster_table)),rep(NA,nrow(raster_table)))

    colnames(raster_table)[c((ncol(raster_table)-1),ncol(raster_table))]=c("present_quarter","mean_start_quarter")
    raster_table[1,ncol(raster_table)-1]=sub("Biomass.","",colnames(kelp_biomass_data)[ncol(kelp_biomass_data)])
    raster_table[1,ncol(raster_table)]=sub("Biomass.","",colnames(kelp_biomass_data)[start_column])
    raster_table[,1]<-round(raster_table[,1],3)
    raster_table[,3]<-round(raster_table[,3],1)
    raster_table[,4]<-round(raster_table[,4],3)

    assign(table_name,raster_table)
    write.table(raster_table,file=table_name,row.names=F, quote=F, sep="\t")

  }

  return(ratio_raster)

}
