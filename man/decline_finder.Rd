\name{decline_finder}
\alias{decline_finder}
\title{
The main function estimating proportion of Landsat pixels with kelp biomass in decline.
}
\description{
The function first estimates if individual Landsat pixels (30x30 m) are declining or not when compared to a historical
baseline (also estimated internally). Then calculates the proportion of pixels in decline inside larger areas of 0.25 x0.25 degree (lat and long).
 The function produces an array that can be useful to map out kelp in decline. The function can also produce an output table with summary
statistics for each region (0.25 x 0.25 degree scale). 
}
\usage{
decline_finder(data, baseline_threshold = 0.1, scarce_cutoff = 0.6, present_window = 16, 
hist_period = 100, window_lag = 0, lat_min = 27.01, lat_max = 37.5, lon_min = -123.5,
       lon_max = -114, table_name = NULL)
}
\arguments{
  \item{data}{
       A data frame of temporal variation of kelp biomass per Landsat pixel, as produced by \code{\link{nc_convert}}. 
  }
  \item{baseline_threshold}{
       A proportion of the average historical biomass (per pixel) below which all quarters in present_window period,  for the same pixel,
       are required to be for the pixel to be classified as in decline. Defaults conservatively to 0.1, i.e., 10 percent of the historical
        average.
        }
  \item{scarce_cutoff}{
       If the proportion of quarters with no kelp in the time series of length hist_period is larger than this cutoff value, 
       the pixel is considered kelp-scarce and removed from analysis. Defaults to 0.6, i.e., pixel needs to have non-zero kelp biomass in at
        least 40 percent of the quarters to be retained in the analysis. 
       }
  \item{present_window}{
       The duration of the current period, in number of quarters, used to characterize present biomass status. 
       }
  \item{hist_period}{ 
       Number of quarters, counting back from the present, used to calculate hist_biomass. Defaults to 100 quarters.
        }
  \item{window_lag}{ 
       How many quarters from the most recent quarter in the time series should 
       present_window slide back. This allows to estimate how kelp in decline was classified for
        periods in the past. It defaults to zero.
       }
  \item{lat_min}{ 
       Used to limit the extent of the output files and maps by providing the minimum latitude in decimal degrees. Defaults to 27.01, the edge of the data in the SBC-LTER file in 2022. 
       }
  \item{lat_max}{
       Used to limit the extent of the output files and maps by providing the maximum latitude in decimal degrees. Defaults to 37.05, the edge of the data in the SBC-LTER file in 2022. 
       }
  \item{lon_min}{ 
       Used to limit the extent of the output files and maps by providing the minimum longitude in decimal degrees. Defaults to -123.5, the edge of the data in the SBC-LTER file in 2022. 
       }
  \item{lon_max}{
       Used to limit the extent of the output files and maps by providing the maximum longitude in decimal degrees. Defaults to -114, the edge of the data in the SBC-LTER file in 2022. 
       }
  \item{table_name}{
       A string with the name of the output file recording summary statistics for each 0.25 x 0.25 degree region in the study extent. This defaults to NULL, so be sure to
       replace with a quoted file name (string) here if you want a detailed output written to the working directory. 
       }
}
\details{
   The function offers several arguments that allow for flexible definition of the periods used for charactherizing historical and present times and to control
   overall functionality. The flexibility  allows for an easy implementation of sensitivity analysis of the effect of these arguments to the characterization of kelp decline.
   Examples are given in a package vignette.}
\value{
   A raster file with the proportion of Landsat pixels in decline for each 0.25 x 0.25 degree region in the study extent. This raster can be easily plotted to map out 
   regions of kelp in decline (see examples in package vignette). To get an optional summary statistic file written to the working directory be sure to supply a name to argument
   table_name.
}
\references{
    Tennies, N., Alberto, F. (in prep) A tool for detecting giant kelp canopy biomass decline in the Californias.
    
    Bell, T, K. Cavanaugh, D. Siegel. (2023) SBC LTER: Time series of quarterly NetCDF files of kelp
biomass in the canopy from Landsat 5, 7 and 8, since 1984 (ongoing) ver 19. Environmental Data Initiative.
https://doi.org/10.6073/pasta/630565d6a8bf54c7cbce6802284dd431. Accessed 2023-02-22.
}
\author{
   Nathaniel Tennies and Filipe Alberto
}
\seealso{
    \code{\link{nc_convert}}
}
\keyword{ IO }
