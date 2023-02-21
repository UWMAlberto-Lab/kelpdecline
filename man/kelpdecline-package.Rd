\name{kelpdecline-package}
\alias{kelpdecline-package}
\alias{kelpdecline}
\docType{package}
\title{
	Estmates kelp biomass decline relative to a historical baseline
}
\description{Estimates proportion of Landsat pixels with kelp biomass in decline when compared to a historical baseline. 
	The proportions are given for 0.25 x 0.25 degree, lat x long, cells in Central, Southern and Baja California.
	}
\details{The nc_convert function first converts the nc file curated by the SBC-LTER project to a data.frame with kelp biomass by Landsat pixel and quarter.
	Then the decline_finder function estimates if individual Landsat pixels (30x30 m) are declining or not when compared to a historical
	 baseline (also estimated internally). The function then calculates the proportion of pixels in decline inside larger areas of 0.25 by
	 0.25 degree (lat and long). This function produces an array that can be useful to map out kelp in decline. The function can also produce an output
	 table with summary statistics for each region (0.25 by 0.25 degree scale).
}
\author{
    Nathaniel Tennies and Fiipe Alberto

Maintainer: Filipe Alberto albertof@uwm.edu
}
\references{
    Tennies N, Alberto F. (in prep) A tool for detecting kelp canopy biomass decline in Central, Southern and Baja California.

    Bell, T. W., Cavanaugh, K. C., and Siegel, D. A. (2017). SBC LTER: Time series of quarterly netcdf files of kelp 
    biomass in the canopy from Landsat 5, 7 and 8, 1984-2016 (ongoing). Santa Barbara Coastal LTER.
}

\keyword{ package }
