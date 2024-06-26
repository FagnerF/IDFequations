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
IDFequations::ScriptIDF(StatesANA=c("ACRE","AMAPÁ"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozona=c()) #ANA database

IDFequations::ScriptIDF(StatesANA=c(),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozona=c()) #Local database 
Example file: https://www.mediafire.com/file/snic3aep3d25eke/Example_database.xlsx/file

IDFequations::ScriptIDFSR(Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozona=c()) #Remote-sensing (“.nc” or “.tif”) 
Example file:https://www.mediafire.com/file/518cx9nykmeo1hr/Example_database.txt/file 

IDFequations::ScriptIDFCLIMBra(CLIMBra=c("ssp245"),Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozona=c()) # CLIMBra database (ssp245 or ssp585)
Example file:https://www.mediafire.com/file/518cx9nykmeo1hr/Example_database.txt/file 

IDFequations::ScriptGraph(Directory="C:/Users/Fagner/Desktop/Results/")

Note: To use this package, are required time series over 20 years.

### Quote
Shortly
