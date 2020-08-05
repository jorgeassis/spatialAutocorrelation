## Spatial Autocorrelation

R Pipelines to reduce the spatial autocorrelation in Species Distribution Models.

Spatial autocorrelation (SA) is a common challenge in Species Distribution Models, which may result in inappropriate spatial inference and prediction. 

Here I propose a straightforward approach to reduce the effect of SA in SDM. I use a simple example focused on a brown algae species capable of producing marine forests and a set of environmental predictors known to largely explain its distribution.

### Dependences

library(raster) <br>
library(ggplot2) <br>
library(sdmpredictors) <br>
library(ecodist) <br>
library(sp) <br>
library(spThin)

## Authors

**Jorge Assis** @ [theMarineDataScientist](github.com/jorgeassis)

## License

Except where otherwise noted, the content on this repository is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).

### Appropriate credits

Assis, J. (2020) R Pipelines to reduce the spatial autocorrelation in Species Distribution Models. theMarineDataScientist
