#' This script calculates the IDF Equations for Brazil
#'
#' @param StatesANA Input Data server ANA ("ACRE") or Local ("")
#' @param Directory Input Location where the IDF equations will be saved
#' @param Method Input Disaggregation method
#' @param Isozone Input Isozone
#'
#' @importFrom stats setNames qnorm qlnorm qgamma qexp qweibull pweibull sd nlminb na.omit
#' @import leaflet
#' @import dplyr
#' @importFrom magrittr %>%
#' @import readxl
#' @import fitdistrplus
#' @import e1071
#' @import evd
#' @import hydroGOF
#' @import xml2
#' @import optimx
#' @import lubridate
#' @import utils
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import tibble
#' @importFrom tidyr pivot_longer
#'
#' @export
ScriptIDF <- function(StatesANA, Directory, Method, Isozone) {
  tmp <- NULL

  ArquivTr <- c(5, 10, 15, 20, 25, 50, 100) # (years)

  Probab <- 1 / ArquivTr

  ArquivDuracoes <-
    c(6, 10, 15, 20, 30, 60, 360, 480, 720, 1440) # (minutes)

  stationType <- 2

  if (length(StatesANA) >= 1) {
    CodStat <- ScriptCodStat(StatesANA, stationType)$CodStat

    State <- vector("character", length = nrow(CodStat))
    Name <- vector("character", length = nrow(CodStat))
    Station <- vector("numeric", length = nrow(CodStat))
    Latitude <- vector("numeric", length = nrow(CodStat))
    Longitude <- vector("numeric", length = nrow(CodStat))
    Equation <- vector("numeric", length = nrow(CodStat))

    cat("\n")

    for (i in 1:nrow(CodStat)) {
      cat("== Working on the data ==\n")
      cat("Wait...\n")
      cat("Station", "", CodStat$codstation[i], "\n")

      # Create list to receive results
      serie <- list()
      Data <- list()
      Ano <- list()
      Prec <- list()

      # Wrap the code accessing station data in tryCatch
      tryCatch(
        {
          # Download the historical series from the ANA data-base
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

          # Data-base station
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

          TabPrecDiaria <- tab %>%
            stats::setNames(c("Data", "PrecDiaria")) %>%
            stats::na.omit()

          # Função para converter e formatar datas
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

          # Aplicar a funcao a uma coluna de um dataframe
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

          # Encontrar intervalos consecutivos
          intervalos <- split(MaxPrec$Ano, cumsum(c(1, diff(MaxPrec$Ano) != 1)))

          # Encontrar o intervalo mais recente
          intervalo_mais_recente <- tail(intervalos, 1)[[1]]

          # Filtrar os dados de precipitação com base no intervalo mais recente
          valores_recente <- MaxPrec$PrecDiaria[MaxPrec$Ano %in% intervalo_mais_recente]

          # Criar o dataframe com os anos mais recentes e as precipitações correspondentes
          df_precipitacao_recente <- data.frame(Ano = intervalo_mais_recente, PrecDiaria = valores_recente)

          row_sub <- apply(df_precipitacao_recente, 1, function(row) {
            all(row != 0)
          })

          ArquivPrec <- df_precipitacao_recente[row_sub, ]

          ArquivPrecNumeric <- as.numeric(df_precipitacao_recente[row_sub, 2])

          ArquivPrecCrescente <- sort(c(ArquivPrecNumeric))

          ArquivPrec <- ArquivPrecCrescente

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
            cat(paste0("ATTENTION, no distribution fits the station ", CodStat$codstation[i], "!"))
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
            State[length(tmp)] <- CodStat$State[i]
            Name[length(tmp)] <- CodStat$Name[i]
            Station[length(tmp)] <- CodStat$codstation[i]
            Latitude[length(tmp)] <- CodStat$lat[i]
            Longitude[length(tmp)] <- CodStat$long[i]
            Equation[length(tmp)] <- best_result$Equation

            cat("\014")
          }
        },
        error = function(e) {
          cat("Error occurred while accessing station data:", conditionMessage(e), "\n")
          cat("Moving to the next station...\n")
        }
      )
    }

    # Criar a tabela apenas com as estações que têm distribuições ajustadas
    if (length(tmp) > 0) {
      # Converter a lista de resultados ajustados em um data frame
      TableFinal <- do.call(rbind, tmp)

      # Adicionar dados de estação à TableFinal
      TableFinal$State <- State[1:length(tmp)]
      TableFinal$Name <- Name[1:length(tmp)]
      TableFinal$Station <- Station[1:length(tmp)]
      TableFinal$Latitude <- Latitude[1:length(tmp)]
      TableFinal$Longitude <- Longitude[1:length(tmp)]
      TableFinal$Equation <- Equation[1:length(tmp)]

      # Reorganizar as colunas
      TableFinal <- TableFinal[, c(12:17, 1:11)]
    } else {
      # Caso contrário, definir TableFinal como NULL
      TableFinal <- NULL
    }

    View(TableFinal)

    TabFinal <- TableFinal %>%
      write.table(
        paste0(Directory, "/", "TabCoefANA.txt"),
        append = FALSE,
        sep = ";",
        dec = ".",
        row.names = FALSE,
        col.names = TRUE
      )
  } else {
    # Create list to receive results
    LatLon <- list()

    path <- file.choose()

    NomesAbas <- readxl::excel_sheets(path)

    LatLon <- readxl::read_xlsx(paste0(path), "LatLon") %>%
      stats::na.omit()

    iteration <- (length(NomesAbas) - 1)

    Station <- vector("character", length = iteration)
    Latitude <- vector("numeric", length = iteration)
    Longitude <- vector("numeric", length = iteration)
    Equation <- vector("numeric", length = iteration)

    for (t in 1:iteration) {
      cat("== Working on the data ==", "\n")
      cat("Wait...", "\n")
      cat("Station", "-", NomesAbas[t + 1], "\n")

      # Create list to receive results
      Dados <- list()

      Dados <- readxl::read_xlsx(paste0(path), NomesAbas[t + 1]) %>%
        stats::na.omit()

      Dados <- Dados[Dados[, 2] != "0", ]
      Dados <- Dados[Dados[, 2] != "0.1", ]

      Dadosdf <- data.frame(Dados[, ])

      DadosdfNumeric <- as.numeric(Dadosdf[, 2])

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
        cat(paste0("ATTENTION, no distribution fits the station ", NomesAbas[t + 1], "!"))
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
        Station[length(tmp)] <- NomesAbas[t + 1]
        Latitude[length(tmp)] <- LatLon$Latitude[t]
        Longitude[length(tmp)] <- LatLon$Longitude[t]
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
        paste0(Directory, "/", "TabCoefLocal.txt"),
        append = FALSE,
        sep = ";",
        dec = ".",
        row.names = FALSE,
        col.names = TRUE
      )

    cat("\014")

    cat("Process finished!", "\n")
  }
}
