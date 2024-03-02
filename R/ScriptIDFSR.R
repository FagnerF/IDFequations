#' This script calculates the IDF equations - remote sensing
#'
#' @param LatSR LonSR Input Latitude
#' @param LonSR Input Longitude
#' @param Station Input Name station
#' @param Directory Location where the IDF equations will be saved
#'
#' @importFrom stats setNames qnorm qlnorm qgamma qexp qweibull pweibull sd nlminb na.omit
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom smwrBase ppearsonIII qpearsonIII plpearsonIII qlpearsonIII
#' @importFrom PearsonDS ppearson0 qpearson0
#' @import remotes
#' @import fitdistrplus
#' @import e1071
#' @import evd
#' @import hydroGOF
#' @import optimx
#' @import lubridate
#' @import utils
#' @import ncdf4
#' @importFrom terra rast extract
#' @importFrom raster stack
#' @importFrom sp SpatialPoints
#' @importFrom tidyr pivot_longer
#'
#' @export
ScriptIDFSR <- function(Station, LatSR, LonSR, Directory) {
  tmp <- NULL

  ArquivTr <- c(2, 5, 10, 15, 20, 25, 50, 75, 100) #(years)
  TamanhoArquivTr <- length(ArquivTr)

  Probab <- 1 / ArquivTr

  ArquivDuracoes <-
    c(5, 10, 15, 20, 30, 60, 360, 480, 720, 1440) #(minutes)
  TamanhoArquivDuracoes <- length(ArquivDuracoes)

  valcrit <- c(0.254)

  stationType <- 2

  path <- file.choose()

  nc <- terra::rast(path)

  NamesDatanc <- raster::stack(path)
  TamSerie <- dim(NamesDatanc)[3]
  NamesData <- matrix(0, TamSerie, 1)
  NamesData <- gsub("X", " ", gsub("[.]", "/", names(NamesDatanc)))

  iteration <- length(Station)

  for (t in 1:iteration) {
    cat("== Working on the data ==", "\n")
    cat("Wait...", "\n")

    #Create list to receive results
    Dados <- list()
    precipitation_1 <- list()
    Data <- list()
    Ano <- list()
    PrecDiaria <- list()

    coord <-
      sp::SpatialPoints(coords = data.frame(x = as.numeric(LonSR[t]), y = as.numeric(LatSR[t])))

    poly <- terra::vect(coord)

    #Extract rainfall
    Dados <- terra::extract(nc, poly)

    colnames(Dados)[-c(1)] <-
      paste0("precipitation_", seq(1, TamSerie, 1))

    Dados <- Dados %>%
      tidyr::pivot_longer(cols = precipitation_1:paste0("precipitation_", TamSerie),
                          values_to = "PrecDiaria")

    Dados <- Dados[, -c(1, 2)]

    TabPrecDiaria <- data.frame(NamesData, Dados) %>%
      setNames(c("Data", "PrecDiaria")) %>%
      data.frame()

    MaxPrec <- TabPrecDiaria %>%
      dplyr::mutate(Ano = lubridate::year(lubridate::ymd(Data))) %>%
      dplyr::group_by(Ano) %>%
      dplyr::slice(which.max(PrecDiaria)) %>%
      dplyr::ungroup() %>%
      data.frame()

    DadosdfNumeric <- as.numeric(MaxPrec$PrecDiaria)

    ArquivPrec <- DadosdfNumeric[order(DadosdfNumeric)]

    TamanhoArquivPrec <- length(ArquivPrec)

    dfTAR <- ScriptTestAd(ArquivPrec)$dfTA$Resultados
    dfTAD <- ScriptTestAd(ArquivPrec)$dfTA$Distribuicoes

    if (dfTAR > valcrit) {
      stop("ATTENTION, no distribution fits the historical rainfall series!")

    } else{
      modnormal <- ScriptTestAd(ArquivPrec)$modnormal
      modln <- ScriptTestAd(ArquivPrec)$modln
      modgamma <- ScriptTestAd(ArquivPrec)$modgamma
      modexpo <- ScriptTestAd(ArquivPrec)$modexpo
      modweibull <- ScriptTestAd(ArquivPrec)$modweibull
      modgumbel <- ScriptTestAd(ArquivPrec)$modgumbel
      modGEV <- ScriptTestAd(ArquivPrec)$modGEV

      X <- ScriptDistProb(
        dfTAD,
        Probab,
        modnormal,
        modln,
        modgamma,
        modexpo,
        modweibull,
        modgumbel,
        modGEV,
        ArquivPrec,
        TamanhoArquivDuracoes,
        TamanhoArquivTr,
        ArquivDuracoes
      )$X

      IMaxObs <- ScriptPMax24h(X,
                               TamanhoArquivDuracoes,
                               TamanhoArquivTr,
                               ArquivDuracoes)$IMaxObs
      kotim <- ScriptOtimiza(TamanhoArquivDuracoes,
                             TamanhoArquivTr,
                             ArquivTr,
                             ArquivDuracoes,
                             IMaxObs)$kotim
      motim <- ScriptOtimiza(TamanhoArquivDuracoes,
                             TamanhoArquivTr,
                             ArquivTr,
                             ArquivDuracoes,
                             IMaxObs)$motim
      t0otim <- ScriptOtimiza(TamanhoArquivDuracoes,
                              TamanhoArquivTr,
                              ArquivTr,
                              ArquivDuracoes,
                              IMaxObs)$t0otim
      notim <- ScriptOtimiza(TamanhoArquivDuracoes,
                             TamanhoArquivTr,
                             ArquivTr,
                             ArquivDuracoes,
                             IMaxObs)$notim

      TabelaResultados <- ScriptExit(
        TamanhoArquivDuracoes,
        TamanhoArquivTr,
        kotim,
        ArquivTr,
        motim,
        ArquivDuracoes,
        t0otim,
        notim,
        ArquivPrec,
        IMaxObs,
        dfTAD
      )$TabFinal

      tmp <- rbind(tmp, TabelaResultados)

    }

    cat("\014")

  }

  TableFinal <- matrix(0, length(Station), 12)

  TableFinal <-
    data.frame(Station, LatSR, LonSR, tmp) %>%
    stats::setNames(
      c(
        "Station",
        "Latitude",
        "Longitude",
        "K",
        "m",
        "t0",
        "n",
        "Distribution",
        "NS",
        "R2",
        "RMSE",
        "MAE"
      )
    )

  View(TableFinal)

  TabFinal <- TableFinal %>%
    write.table(
      paste0(Directory, "/", "TabCoef.txt"),
      append = FALSE,
      sep = ";",
      dec = ".",
      row.names = FALSE,
      col.names = TRUE
    )

  mapa <- leaflet(data = TableFinal[,-c(8:12)]) %>%
    addTiles() %>%
    setView(lng = -49.265,
            lat = -10.939,
            zoom = 4)

  mapa <- mapa %>%
    addMarkers(
      lng = ~ as.numeric(LonSR),
      lat = ~ as.numeric(LatSR),
      label = ~ paste(
        "Station",
        "-" ,
        Station,
        ":" ,
        "*",
        "K:",
        K,
        "*",
        "m:",
        m,
        "*",
        "t0:",
        t0,
        "*",
        "n:",
        n
      )
    )

  return(mapa)

  cat("\014")

  cat("Process finished!", "\n")

}
