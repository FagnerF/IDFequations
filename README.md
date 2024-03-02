# 👋 Hello! Welcome the function IDF equations for Brazil  <p align="right">

## My name is Fagner Costa

### Description
This package calculates IDF equations for Brazil. This package is the result
of the thesis research of the student Fagner Costa

#### What is the IDF equation used for?

The IDF curve, or their respective equations, relates Intensity, Duration, and Frequency and is widely used in engineering to determine the maximum rainfall in a given location.
Various hydraulic structures dimensioned are from this rain, such as dams, channels, retention basins, gutters, etc.
The determination of each curve depends on local rainfall; it is impossible to use an IDF curve developed for a city in Minas Gerais in a municipality in São Paulo.

#### But why?

As the rainfall regime varies from region to region, and the IDF curves are strictly related to rainfall, they also vary.
Usually, IDFs built are from observed historical data, assuming that, in the future, the same trend will be in situ (stationarity).

### Pre-requisites
Version R >= 4.3.2 (https://cran.r-project.org/)

Software RStudio (https://posit.co/downloads/)

### Installation package
install.packages("devtools")

library(devtools)

devtools::install_github("FagnerF/IDFequations")

library(IDFequations)

### Examples
IDFequations::ScriptIDF(StatesANA=c("ACRE","AMAPÁ"),Directory="C:/Users/Fagner/Desktop/Results/") #ANA data-base

IDFequations::ScriptIDF(StatesANA=c(),Directory="C:/Users/Fagner/Desktop/Results/") #Local data-base

IDFequations::ScriptIDFSR(Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/") #Remote sensing (MSWEP, CHIRPS, CPC, for examples) data-base. Simply provide, in a single file, the historical series in “.nc” or “.tif” formats.

IDFequations::ScriptIDFCLIMBra(CLIMBra=c("ssp245"),Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/") # CLIMBra data-base (CLIMBra= ssp245 or ssp585)

IDFequations::ScriptGraphic(Directory=c("C:/Users/Fagner/Desktop/Results/")) #Graphics

Note: To use this package, historical rainfall records of more than 20 years are required.

### Quote
Shortly
