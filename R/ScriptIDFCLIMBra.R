#' This script calculates the IDF Equations - CLIMBra
#'
#' @param CLIMBra Input Type GCMs (ssp245 or ssp585)
#' @param Station Input Name station
#' @param LatSR Input Latitude
#' @param LonSR Input Longitude
#' @param Directory Input Location where the IDF equations will be saved#'
#' @param Method Input Disaggregation method
#' @param Isozona Input Isozona
#'
#' @importFrom stats setNames qnorm qlnorm qgamma qexp qweibull pweibull sd nlminb na.omit
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom smwrBase ppearsonIII qpearsonIII plpearsonIII qlpearsonIII
#' @importFrom PearsonDS ppearson0 qpearson0
#' @import fitdistrplus
#' @import e1071
#' @import evd
#' @import hydroGOF
#' @import optimx
#' @import lubridate
#' @import utils
#' @importFrom terra rast extract
#' @importFrom raster stack
#' @importFrom sp SpatialPoints
#' @importFrom tidyr pivot_longer
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import tibble
#'
#' @export
ScriptIDFCLIMBra <-
  function(CLIMBra, Station, LatSR, LonSR, Directory, Method, Isozona) {
    tmp <- NULL
    Station <- vector("character", length = length(Station))
    Latitude <- vector("numeric", length = length(Station))
    Longitude <- vector("numeric", length = length(Station))

    ArquivTr <- c(2, 5, 10, 15, 20, 25, 50, 75, 100) # (years)
    TamanhoArquivTr <- length(ArquivTr)

    Probab <- 1 / ArquivTr

    ArquivDuracoes <-
      c(5, 10, 15, 20, 30, 60, 360, 480, 720, 1440) # (minutes)
    TamanhoArquivDuracoes <- length(ArquivDuracoes)

    valcrit <- c(0.254)

    stationType <- 2

    # URL dependendo do cenário
    url <- switch(CLIMBra,
      "ssp245" = paste0(base_url, "?fileId=4a9e943f253195a7b10eaaaeec55c83e&path=/V5/Gridded%20data/pr/ssp245/EC-EARTH3-pr-ssp245.nc&fileName=EC-EARTH3-pr-ssp245.nc"),
      "ssp585" = paste0(base_url, "?fileId=dc8a7220ea27f834e1c2c4fd98853b5a&path=/V5/Gridded%20data/pr/ssp585/EC-EARTH3-pr-ssp585.nc&fileName=EC-EARTH3-pr-ssp585.nc")
    )

    # Extrai o nome do arquivo da URL usando expressão regular
    filename <- sub(".*fileName=([^&]+).*", "\\1", url)

    # Download do arquivo
    utils::download.file(url, filename, mode = "wb", method = "curl")

    # Carrega o arquivo NetCDF
    nc <- tryCatch({
      terra::rast(filename)
    })

    # Extrai nomes das variáveis do arquivo NetCDF
    NamesDatanc <- tryCatch({
      raster::stack(filename)
    })
    TamSerie <- dim(NamesDatanc)[3]
    NamesData <- matrix(0, TamSerie, 1)
    NamesData <- gsub("X", " ", gsub("[.]", "/", names(NamesDatanc)))

    iteration <- length(Station)

    for (t in 1:iteration) {
      cat("== Working on the data ==", "\n")
      cat("Wait...", "\n")

      # Create list to receive results
      Dados <- list()
      precipitation_1 <- list()
      Data <- list()
      Ano <- list()
      PrecDiaria <- list()

      coord <-
        sp::SpatialPoints(coords = data.frame(x = as.numeric(LonSR[t]), y = as.numeric(LatSR[t])))

      poly <- terra::vect(coord, crs = "+proj=longlat +datum=WGS84")

      # Extract rainfall
      Dados <- terra::extract(nc, poly)

      colnames(Dados)[-c(1)] <-
        paste0("precipitation_", seq(1, TamSerie, 1))

      Dados <- Dados %>%
        tidyr::pivot_longer(
          cols = precipitation_1:paste0("precipitation_", TamSerie),
          values_to = "PrecDiaria"
        )

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

      modnormal <- ScriptTestAd(ArquivPrec)$modnormal
      modln <- ScriptTestAd(ArquivPrec)$modln
      modgamma <- ScriptTestAd(ArquivPrec)$modgamma
      modexpo <- ScriptTestAd(ArquivPrec)$modexpo
      modweibull <- ScriptTestAd(ArquivPrec)$modweibull
      modgumbel <- ScriptTestAd(ArquivPrec)$modgumbel
      modGEV <- ScriptTestAd(ArquivPrec)$modGEV

      if (is.null(modnormal) && is.null(modln) && is.null(modgamma) &&
          is.null(modexpo) && is.null(modweibull) && is.null(modgumbel) && is.null(modGEV)) {
        cat(paste0("ATTENTION, no distribution fits the station ", Station[t], "!"))
        # Pass to the next station
        next
      } else {
        # Calcular distribuições acumuladas
        dfTAD <- ScriptTestAd(ArquivPrec)$dfTAD$Distribuicoes
        X <- ScriptDistProb(
          dfTAD,
          Probab,
          modnormal,
          modln,
          modgamma,
          modexpo,
          modweibull,
          modgumbel,
          modGEV
        )$X

        IMaxObs <- ScriptPMax24h(
          X,
          TamanhoArquivDuracoes,
          TamanhoArquivTr,
          ArquivDuracoes,
          ArquivTr,
          Method,
          Isozona
        )$IMaxObs
        kotim <- ScriptOtimiza(
          TamanhoArquivDuracoes,
          TamanhoArquivTr,
          ArquivTr,
          ArquivDuracoes,
          IMaxObs
        )$kotim
        motim <- ScriptOtimiza(
          TamanhoArquivDuracoes,
          TamanhoArquivTr,
          ArquivTr,
          ArquivDuracoes,
          IMaxObs
        )$motim
        t0otim <- ScriptOtimiza(
          TamanhoArquivDuracoes,
          TamanhoArquivTr,
          ArquivTr,
          ArquivDuracoes,
          IMaxObs
        )$t0otim
        notim <- ScriptOtimiza(
          TamanhoArquivDuracoes,
          TamanhoArquivTr,
          ArquivTr,
          ArquivDuracoes,
          IMaxObs
        )$notim

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

        # Armazenar resultados apenas se a estação tiver uma distribuição ajustada
        tmp[[length(tmp) + 1]] <- TabelaResultados

        # Armazenar dados da estação
        Station[length(tmp)] <- Station[t]
        Latitude[length(tmp)] <- LatSR[t]
        Longitude[length(tmp)] <- LonSR[t]
      }

      cat("\014")
    }

    # Criar a tabela apenas com as estações que têm distribuições ajustadas
    if (length(tmp) > 0) {
      # Converter a lista de resultados ajustados em um data frame
      TableFinal <- do.call(rbind, tmp)

      # Adicionar dados de estação à TableFinal
      TableFinal$Station <- Station[1:length(tmp)]
      TableFinal$Latitude <- Latitude[1:length(tmp)]
      TableFinal$Longitude <- Longitude[1:length(tmp)]

      # Reorganizar as colunas para que as três últimas se tornem as três primeiras
      TableFinal <- TableFinal[, c(10:12, 1:9)]
    } else {
      # Caso contrário, definir TableFinal como NULL
      TableFinal <- NULL
    }

    View(TableFinal)

    TabFinal <- TableFinal %>%
      write.table(
        paste0(Directory, "/", "TabCoefScenario.txt"),
        append = FALSE,
        sep = ";",
        dec = ".",
        row.names = FALSE,
        col.names = TRUE
      )

    mapa <- leaflet(data = TableFinal[, -c(8:12)]) %>%
      addTiles() %>%
      setView(
        lng = -49.265,
        lat = -10.939,
        zoom = 4
      )

    mapa <- mapa %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        label = ~ paste(
          "Station",
          "-",
          Station,
          ":",
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
