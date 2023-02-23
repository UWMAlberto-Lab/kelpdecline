\name{nc_convert}
\alias{nc_convert}
\title{
Convert nc input to a data.frame with kelp biomass variation
}
\description{
The function converts the external file in nc format to a data.frame with the time series of kelp biomass.
}
\usage{
nc_convert(nc_data_location)
}
\arguments{
  \item{nc_data_location}{ A string with the name of the nc file file containing kelp biomass temporal variation, curated by the SBC-LTER team.
}
}

\value{
    A data frame with kelp biomass recorded for each 30 x 30m Landsat pixel areas in the study extent (Central to Baja California).
    Each pixel with kelp is a row in the data frame, whereas columns contain temporal variation recorded by quarter.
}
\references{
 Tennies, N., Alberto, F. (in prep) A tool for detecting giant kelp canopy biomass decline in the Californias.
 
 Bell, T, K. Cavanaugh, D. Siegel. (2023) SBC LTER: Time series of quarterly NetCDF files of kelp biomass in the canopy from Landsat 5, 7 and 8, since 1984 (ongoing) ver 19. Environmental Data Initiative.
https://doi.org/10.6073/pasta/630565d6a8bf54c7cbce6802284dd431. Accessed 2023-02-22.
}
\author{
Nathaniel Tennies and Filipe Alberto
}
\seealso{
\code{\link{decline_finder}}
}
\keyword{ IO }
