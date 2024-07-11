# 👋 Hello! Welcome the tool GEOT-IDF Equations <p align="right">

## My name is Fagner Costa

### Prerequisites
R version = 4.4.1 (https://cran.r-project.org/)

Rtools version = 4.4 (https://cran.r-project.org/bin/windows/Rtools/)

Software RStudio (https://posit.co/downloads/)

System operational: Windows

### How to install
install.packages("devtools")

library(devtools)

devtools::install_github("FagnerF/IDFequations")

library(IDFequations)

### How to run
IDFequations::ScriptIDF(StatesANA=c("ACRE","AMAPÁ"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozone=c()) #ANA database

IDFequations::ScriptIDF(StatesANA=c(),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozone=c()) #Local database 
Example file: https://www.mediafire.com/file/snic3aep3d25eke/Example_database.xlsx/file

IDFequations::ScriptIDFSR(Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozone=c()) #Remote-sensing (“.nc” or “.tif”) 

IDFequations::ScriptIDFCLIMBra(CLIMBra=c("ssp245"),Station=c("Agua Branca-PB"),LatSR=c("-7.512"),LonSR=c("-37.637"),Directory="C:/Users/Fagner/Desktop/Results/",Method=c("CETESB"),Isozone=c()) # CLIMBra database (ssp245 or ssp585)

Note: To use this package, are required time series over 20 years.

### How to cite
Costa F. F., Rufino I. A. A., Aragão R., Castro M. A. H., Filho R. S. R. (2024). GEOT-IDF Equations: An R-based tool for intense rainfall studies in Brazil. Earth Science Informatics. https://doi.org/10.21203/rs.3.rs-4172597/v1
