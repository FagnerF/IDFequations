#' This script calculates the IDF Equations - CLIMBra
#'
#' @param CLIMBra Input Type GCMs (ssp245 or ssp585)
#' @param Station Input Name station
#' @param LatSR Input Latitude
#' @param LonSR Input Longitude
#' @param Directory Input Location where the IDF equations will be saved#'
#' @param Method Input Disaggregation method
#' @param Isozone Input Isozone
#'
#' @importFrom stats setNames qnorm qlnorm qgamma qexp qweibull pweibull sd nlminb na.omit
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %>%
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
  function(CLIMBra, Station, LatSR, LonSR, Directory, Method, Isozone) {
    tmp <- NULL
    TableFinal <- NULL
    Station <- vector("character", length = length(Station))
    Latitude <- vector("numeric", length = length(Station))
    Longitude <- vector("numeric", length = length(Station))
    Equation <- vector("numeric", length = length(Station))

    ArquivTr <- c(5, 10, 15, 20, 25, 50, 100) # (years)

    Probab <- 1 / ArquivTr

    ArquivDuracoes <-
      c(6, 10, 15, 20, 30, 60, 360, 480, 720, 1440) # (minutes)

    # URL base
    base_url <- "https://download.scidb.cn/download"

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

      # Funcao para converter e formatar datas
      formatar_data <- function(data_input) {
        # Tentativas de conversao usando varios formatos comuns
        formatos <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%d-%m-%Y")

        # Tentativa de conversao da entrada em cada formato
        for (formato in formatos) {
          tryCatch({
            data_convertida <- as.Date(data_input, format = formato)
            if (!is.na(data_convertida)) {
              return(format(data_convertida, "%d/%m/%Y"))
            }
          }, error = function(e) {
            # Ignorar erros e tentar o proximo formato
          })
        }

        # Caso a conversao falhe, retornar uma mensagem de erro
        stop("Formato de data nao reconhecido")
      }

      # Aplicar a função a uma coluna de um dataframe
      padronizar_datas <- function(df) {
        df$Data <- sapply(df$Data, formatar_data)
        return(df)
      }

      # Aplicar a padronizacao ao dataframe
      dados_padronizados <- padronizar_datas(TabPrecDiaria)

      MaxPrec <- dados_padronizados %>%
        dplyr::mutate(Ano = lubridate::year(lubridate::dmy(Data))) %>%
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
        )

        # Calcular valores máximos observados
        IMaxObs <- ScriptPMax24h(
          X,
          ArquivDuracoes,
          ArquivTr,
          Method,
          Isozone
        )

        # Executar a otimização para cada equação e selecionar o melhor resultado
        best_result <- NULL
        best_NS <- -Inf # Inicializar com valor muito baixo para maximização

        for (eq_number in 1:5) {
          result <- ScriptOtimiza(eq_number, ArquivTr, ArquivDuracoes, IMaxObs)

          # Verificar se result não é NULL e se NS é maior que best_NS
          if (!is.null(result) && result$NS > best_NS) {
            best_result <- result
            best_NS <- result$NS
          }
        }

        # Gerar tabela final de resultados
        TabFinal <- ScriptExit(
          ArquivDuracoes,
          ArquivTr,
          best_result,
          IMaxObs,
          dfTAD
        )

        # Armazenar resultados apenas se a estação tiver uma distribuição ajustada
        tmp[[length(tmp) + 1]] <- TabFinal

        # Armazenar dados da estação
        Station[length(tmp)] <- Station[t]
        Latitude[length(tmp)] <- LatSR[t]
        Longitude[length(tmp)] <- LonSR[t]
        Equation[length(tmp)] <- best_result$Equation
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
      TableFinal$Equation <- Equation[1:length(tmp)]

      # Reorganizar as colunas para que as três últimas se tornem as três primeiras
      TableFinal <- TableFinal[, c(12:15, 1:11)]
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

    cat("\014")

    cat("Process finished!", "\n")
  }
