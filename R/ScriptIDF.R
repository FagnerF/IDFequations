#' This script calculates the IDF equations for any place in Brazil
#'
#' @param StatesANA Input data server ANA or local
#' @param Directory Location where the IDF equations will be saved
#'
#' @importFrom stats setNames qnorm qlnorm qgamma qexp qweibull pweibull sd nlminb na.omit
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom smwrBase ppearsonIII qpearsonIII plpearsonIII qlpearsonIII
#' @importFrom PearsonDS ppearson0 qpearson0
#' @import readxl
#' @import remotes
#' @import fitdistrplus
#' @import e1071
#' @import evd
#' @import hydroGOF
#' @import xml2
#' @import optimx
#' @import lubridate
#' @import utils
#'
#' @export
ScriptIDF <- function(StatesANA, Directory) {
  tmp <- NULL

  ArquivTr <- c(2, 5, 10, 15, 20, 25, 50, 75, 100) #(years)
  TamanhoArquivTr <- length(ArquivTr)

  Probab <- 1 / ArquivTr

  ArquivDuracoes <-
    c(5, 10, 15, 20, 30, 60, 360, 480, 720, 1440) #(minutes)
  TamanhoArquivDuracoes <- length(ArquivDuracoes)

  valcrit <- c(0.254)

  stationType <- 2

  if (length(StatesANA) >= 1) {
    CodStat <- ScriptCodStat(StatesANA, stationType)$CodStat

    cat("\n")

    for (i in 1:nrow(CodStat)) {
      cat("== Working on the data ==", "\n")
      cat("Wait...", "\n")
      cat("Station", "", CodStat$codstation[i], "\n")

      #Create list to receive results
      serie <- list()
      Data <- list()
      Ano <- list()
      Prec <- list()

      #Download the historical series from the ANA data-base
      station_number <- gsub(" ", "%20", CodStat$codstation[i])

      html_raw <- xml2::read_html(
        paste(
          "http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroSerieHistorica?codEstacao=",
          station_number,
          "&dataInicio=&dataFim=&tipoDados=",
          stationType,
          "&nivelConsistencia=",
          sep = ""
        )
      )

      station_df <- html_raw %>%
        xml2::xml_find_all(".//documentelement") %>%
        xml2::xml_children() %>%
        xml2::as_list()

      station_df <- lapply(station_df, function(row) {
        row[sapply(row, function(x) {
          length(x) == 0
        })] <- NA
        row %>%
          unlist() %>%
          t() %>%
          dplyr::as_tibble()
      }) %>%
        do.call(what = dplyr::bind_rows)

      #Data-base station
      station_df <- station_df %>%
        # Convert precipitation columns to numeric and datahora to date format
        dplyr::mutate(
          dplyr::across(dplyr::matches("chuva..$"), as.numeric),
          dplyr::across(dplyr::matches("data"), as.Date)
        ) %>%
        dplyr::rename(data = dplyr::any_of("datahora")) %>%
        dplyr::arrange(dplyr::across("data"))

      serie[[i]] <- station_df

      tab <-
        data.frame(serie[[i]][["data"]], serie[[i]][["maxima"]])

      dados <- tab %>%
        stats::setNames(c("Data", "Prec")) %>%
        stats::na.omit()

      MaxPrec <- dados %>%
        dplyr::mutate(Ano = lubridate::year(lubridate::ydm(Data))) %>%
        dplyr::group_by(Ano) %>%
        dplyr::slice(which.max(Prec)) %>%
        dplyr::ungroup() %>%
        data.frame()

      row_sub <- apply(MaxPrec, 1, function(row)
        all(row != 0))

      ArquivPrec <- MaxPrec[row_sub,]

      ArquivPrecNumeric <- as.numeric(MaxPrec[row_sub, 2])

      ArquivPrecCrescente <- sort(c(ArquivPrecNumeric))

      ArquivPrec <- ArquivPrecCrescente

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

    TableFinal <- matrix(0, length(CodStat$codstation), 12)

    TableFinal <-
      data.frame(CodStat$codstation, CodStat$lat, CodStat$long, tmp) %>%
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
        lng = ~ CodStat$long,
        lat = ~ CodStat$lat,
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

  } else{
    path <- file.choose()

    NomesAbas <- readxl::excel_sheets(path)

    iteration <- (length(NomesAbas) - 1)

    for (t in 1:iteration) {
      cat("== Working on the data ==", "\n")
      cat("Wait...", "\n")
      cat("Station", "-", NomesAbas[t + 1], "\n")

      #Create list to receive results
      Dados <- list()
      LatLon <- list()

      Dados <- readxl::read_xlsx(paste0(path), NomesAbas[t + 1]) %>%
        stats::na.omit()

      LatLon <- readxl::read_xlsx(paste0(path), "LatLon") %>%
        stats::na.omit()

      Dados <- Dados[Dados[, 2] != "0", ]
      Dados <- Dados[Dados[, 2] != "0.1", ]

      Dadosdf <- data.frame(Dados[, ])

      DadosdfNumeric <- as.numeric(Dadosdf[, 2])

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

    TableFinal <- matrix(0, length(NomesAbas[-c(1)]), 12)

    TableFinal <-
      data.frame(NomesAbas[-c(1)],  LatLon$Latitude, LatLon$Longitude, tmp) %>%
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
        lng = ~ LatLon$Longitude,
        lat = ~ LatLon$Latitude,
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
}
