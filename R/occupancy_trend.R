occupancy_trend<-function(data, present_year=2022, outFile="Out.DF.txt",test=FALSE,npermuts=1000){
requireNamespace("sp")
requireNamespace("raster")
	
#Number of quarters in the data set
nquarters<-(ncol(data)-2)
Years<-matrix(unlist(strsplit(names(data[3:(nquarters+2)]),"\\.")),ncol=3,byrow=T)[,2]

	if(!present_year%in%(2000:max(Years))stop("Error: Check that present_year is numeric and not older than 2000")


YearFilter<-Years<=present_year
Years<-Years[YearFilter]

kelp_biomass_data<-data[,c(T,T,YearFilter)]
nquarters<-(ncol(kelp_biomass_data)-2)



#Finding the indexes per cell
latV<-seq(27,37,.25)
longV<-seq(-123.5,-113.5,.25)

print("allocating pixels to 0.25 x 0.25 degree regions")
CellIndexes<-list()
CellNames<-NULL

for(long in longV){
	for(lat in latV){
	temp<-which(kelp_biomass_data$Long>long & kelp_biomass_data$Long<=(long+.25) &
			 kelp_biomass_data$Lat>lat & kelp_biomass_data$Lat<=(lat+0.25))
        if(length(temp)!=0) {
       		 	    CellIndexes<-append(CellIndexes,list(temp))
       		 	     CellNames<-c(CellNames,paste(lat,long,sep=";"))
        		    }
 }
}
names(CellIndexes)<-CellNames

#Total number of cells
print("Estimating present occupancy")
PresentOccup<-apply(kelp_biomass_data[ ,  (nquarters -1):(nquarters+2)  ],1,function(x) as.integer(any(x!=0,na.rm=T)) )

#Estimate the probability of pixel occurrence. This is simply how many years the pixel was present in at least one quarter
#First transform the data to yearly and binomial. First a vector of years
Years<-matrix(unlist(strsplit(names(kelp_biomass_data[3:(nquarters+2)]),"\\.")),ncol=3,byrow=T)[,2]

print("Estimating the Long term annual probability of pixel occupancy")
KelpPresentY.Freq<-apply(kelp_biomass_data[, 3:(nquarters+2)], 1 ,function(x)sum( tapply(x,Years,function(y)any(y!=0)),na.rm=T))

KelpYearsWithData<-apply(kelp_biomass_data[, 3:(nquarters+2)], 1 ,function(x)sum( tapply(x,Years,function(y){!all(is.na(y))}) ))

#This is the proportion of years where kelp was observed in at least one quarter
HistPropYearOc<-KelpPresentY.Freq/ KelpYearsWithData

#Deficit, this is the difference between occupancy in the present and proportion (probability of pixel being occupied each year)
PresObs.Minus.HistYearOcc<-PresentOccup-HistPropYearOc

#Now lets sum the PresObs.Minus.HistYearOcc by Region
#First a data.frame with a row per region (cells)
if(test){  #If condition for using a significance test or not
ProbDeclineDF<-data.frame(matrix(nrow=length(CellIndexes),ncol=10))
names(ProbDeclineDF)<-c("Lat","Long","N.Pixel","AVG.Biomass","RYPO","LTPAPO","POT","LowerQ","UpperQ","Sig")

for(cell in 1: length(CellIndexes)){
   coords<-unlist(strsplit(names(CellIndexes)[cell],";"))
   ProbDeclineDF$Lat[cell]<-as.numeric(coords[1])
   ProbDeclineDF$Long[cell]<-as.numeric(coords[2])
   ProbDeclineDF$N.Pixel[cell]<-length(CellIndexes[[cell]])
   ProbDeclineDF$AVG.Biomass[cell]<-round(sum(apply(kelp_biomass_data[ CellIndexes[[cell]] , (nquarters -1):(nquarters+2) ],1,max),na.rm=T),1)
   ProbDeclineDF$POT[cell]<-round(mean(PresObs.Minus.HistYearOcc[CellIndexes[[cell]]],na.rm=T),4)
   ProbDeclineDF$RYPO[cell]<-round(mean( PresentOccup[CellIndexes[[cell]]]),4)
   ProbDeclineDF$LTPAPO[cell]<-round(mean(HistPropYearOc[CellIndexes[[cell]]]),4)
   #SignTest
     	print(paste0("Testing significance for cell ", cell))
       	testSample<-sapply(HistPropYearOc[CellIndexes[[cell]]],function(p) sample(x=c(1,0), replace=T,size=npermuts,  prob=c(p,1-p)))
        RandBalance<-apply(testSample, 1, function(r) r-HistPropYearOc[CellIndexes[[cell]]])
       	avg.balance.rand<-apply((t(RandBalance)),1,mean)
   ProbDeclineDF$LowerQ[cell]<-round(quantile(avg.balance.rand,0.001),4)
   ProbDeclineDF$UpperQ[cell]<-round(quantile(avg.balance.rand,0.999),4)
   ProbDeclineDF$Sig[cell]<-ifelse(ProbDeclineDF$POT[cell]<=ProbDeclineDF$LowerQ[cell] | ProbDeclineDF$POT[cell]>=ProbDeclineDF$UpperQ[cell],
        				"Sign.","NS")
   print(paste0("Completed test for cell ",cell, " out of ",length(CellIndexes)))
   }
   ProbDeclineDF$Direction<-ifelse(ProbDeclineDF$Sig=="Sign.",  ifelse(ProbDeclineDF$POT>0,"UP","DOWN"),"NS")}else
   {
   ProbDeclineDF<-data.frame(matrix(nrow=length(CellIndexes),ncol=7))
   names(ProbDeclineDF)<-c("Lat","Long","N.Pixel","AVG.Biomass","RYPO","LTPAPO","POT")
   npermuts=5000

for(cell in 1: length(CellIndexes)){
  coords<-unlist(strsplit(names(CellIndexes)[cell],";"))
  ProbDeclineDF$Lat[cell]<-as.numeric(coords[1])
  ProbDeclineDF$Long[cell]<-as.numeric(coords[2])
  ProbDeclineDF$N.Pixel[cell]<-length(CellIndexes[[cell]])
  ProbDeclineDF$AVG.Biomass[cell]<-round(sum(apply(kelp_biomass_data[ CellIndexes[[cell]] , (nquarters -1):(nquarters+2) ],1,max),na.rm=T),1)
  ProbDeclineDF$POT[cell]<-round(mean(PresObs.Minus.HistYearOcc[CellIndexes[[cell]]],na.rm=T),4)
  ProbDeclineDF$RYPO[cell]<-round(mean( PresentOccup[CellIndexes[[cell]]]),4)
  ProbDeclineDF$LTPAPO[cell]<-round(mean(HistPropYearOc[CellIndexes[[cell]]]),4)
  #print(paste0("Completed cell ",cell, " out of ",length(CellIndexes)))
   }

  ProbDeclineDF$Direction<-ifelse(ProbDeclineDF$POT>0,"UP","DOWN")
  }
  o1<-order( ProbDeclineDF$Lat,decreasing=T)
  ProbDeclineDF<-ProbDeclineDF[o1,]

  #writing to file
  write.table(ProbDeclineDF,outFile,quote=F,row.names=F,sep="\t")

  RasterBalance<-raster::rasterFromXYZ(data.frame(x=ProbDeclineDF$Long+0.125,y=ProbDeclineDF$Lat+0.125,z=ProbDeclineDF$POT))
  crs(RasterBalance)<-"+proj=longlat +datum=WGS84"

 return(RasterBalance)
}



