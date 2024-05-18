\name{occupancy_trend}
\alias{occupancy_trend}
\title{
Estimate trend of kelp canopy by comparing present with historical pixel occupancy.
}
\description{
The function estimates a trend of kelp canopy by comparing present (yearly only) and historical occupancies. Results are summarized at the 0.25 x 0.25 degree (regions) scale.
}
\usage{
occupancy_trend(data,present_year,outFile="Out.DF.txt",test=FALSE,npermuts=1000)
}
\arguments{
     \item{data}{
          A data frame of temporal variation of kelp biomass per Landsat pixel, as produced by \code{\link{nc_convert}}. 
          }
     \item{present_year}{
          The year used to characterize the period used to estimate a trend.
          }
     \item{outFile}{
          The name for the tab-delimited output file written to the working directory. Defaults to Out.DF.txt.
          }
     \item{test}{
          A logical to decide if significance tests should be run. Defaults to FALSE.
          }
\item{npermuts}{
         Number of randomizations used in the significance test.
         }
}
\details{
   Pixel occupancy in a given reference year (present_year) is compared to the long-run proportion of years occupied across the entire time series to determine a trend
    (recovery, decline, no significant change). Landsat pixel data (30 x 30 m) kelp yearly occupancy is first converted to an annual, binomial series of pixel occupancy. 
    If at least one quarter in a given year has non-zero biomass, the pixel is considered occupied during the year. Next, long-term probability of yearly pixel occupancy 
    (LTPYPO) is subtracted from the reference year pixel occupancy RYPO (1 or 0), RYPO-LTPYPO. This simple statistic, which ranges from -1 to 1, is a good weight of the
    importance of the change relative to the long-term occupancy history. For example, an occupied pixel (RYPO of 1) with a LTPYPO of 0.8 produces a 0.2 deficit, whereas
    if LTPYPO is only 0.1, the deficit is 0.9. Pixel deficits are averaged for all pixels in a region of 0.25 x 0.25 degrees to produce a single trend value at this scale.
    A significance test is available, where randomized RYPO values are sampled with LTPYPO chances. The procedure is repeated npermuts times to generate a randomized 
    distribution for the region deficit trend under the null hypothesis that pixel occupancy follows LTPYPO. The test does not integrate the autoregressive nature of kelp
    occupancy (the effect of the previous year's occupancy) and should be improved in the future.
    
     
   }
\value{
     A raster file with the trend deficit per region 0.25 x 0.25 degree region.  A summary statistic file is also written to the working directory.
   table_name.
}
\references{
 Tennies, N., Alberto, F. (in prep) A tool for detecting giant kelp canopy biomass decline in the Californias.
 
 Bell, T, K. Cavanaugh, D. Siegel. (2023) SBC LTER: Time series of quarterly NetCDF files of kelp biomass in the canopy from Landsat 5, 7 and 8, since 1984 (ongoing) ver 19. Environmental Data Initiative.
https://doi.org/10.6073/pasta/630565d6a8bf54c7cbce6802284dd431. Accessed 2023-02-22.
}
\author{
Filipe Alberto and Nathaniel Tennies
}
\seealso{
\code{\link{decline_finder}}
\code{\link{nc_convert}}
}
\keyword{ IO }
