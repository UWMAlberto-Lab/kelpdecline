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
Tennies N, Alberto F. (in prep) A tool for detecting kelp canopy biomass decline in Central, Southern and Baja California
}
\author{
Nathaniel Tennies and Filipe Alberto
}
\seealso{
\code{\link{decline_finder}}
}
\keyword{ IO }
