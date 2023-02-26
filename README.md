# kelpdecline
R code to estimate regions of kelp biomass decline in the Californias when compared to historical baselines from Landsat. 

The kelpdecline package aims to detect regions of kelp biomass in decline. The detection is based on a time series of Landsat estimated biomass (Cavanaugh et al. 2011, Bell et al. 2017). The figures and tables produced
quantify the proportion of Landsat pixels (30 x 30 m) with biomass in decline within 0.25 x 0.25 degrees, Lat x Long, regions. The method proposed provides a framework to understand kelp decline in the context
of a historical baseline and is valuable for all kelp forest conservation purposes. A second approach based on Landsat pixel kelp occupancy is also available. This last method provides a trend estimate (-1 to 1) for a reference year. Our specific motivation for building this tool was to design an efficient sampling strategy for genetic monitoring studies (Klingbeil et
al. 2022). Further considerations are found in the companion publication describing this package (Tennies &
Alberto, in prep).

Please read the [Example_kelpdecline_run.pdf](https://github.com/falberto73/kelpdecline/blob/main/Example_kelpdecline_run.pdf) for more a detailed example on the package.
