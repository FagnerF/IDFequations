ScriptCodStat <- function(StatesANA, stationType) {
  # Inicializa a lista CodStat
  CodStat <- list()

  dfTab <- NULL

  cat("Downloading... \n")

  # Configuração de paralelismo
  num_cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores = num_cores)

  # Define a função para recuperar dados das estações
  get_station_data <- function(estadoG, stationType) {
    html_raw <- tryCatch(
      {
        xml2::read_html(
          paste0(
            "http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroInventario?codEstDE=&codEstATE=&tpEst=",
            stationType,
            "&nmEst=&nmRio=&codSubBacia=&codBacia=&nmMunicipio=&nmEstado=",
            estadoG,
            "&sgResp=&sgOper=&telemetrica="
          )
        )
      },
      error = function(e) {
        cat("Error downloading data for", estadoG, ". Moving to the next state.\n")
        return(NULL)
      }
    )
    return(html_raw)
  }

  # Loop para baixar os dados das estações em paralelo
  foreach(u = StatesANA, .combine = "c") %do% {
    estadoG <- gsub(" ", "%20", u)
    html_raw <- get_station_data(estadoG, stationType)

    if (!is.null(html_raw)) {

      # Scrapping info for each station
      estac <- as.data.frame(cbind(
        xml2::xml_text(xml2::xml_contents(xml2::xml_find_all(html_raw, ".//nmestado"))),
        xml2::xml_double(xml2::xml_contents(xml2::xml_find_all(html_raw, ".//codigo"))),
        xml2::xml_text(xml2::xml_contents(xml2::xml_find_all(html_raw, ".//nome"))),
        xml2::xml_double(xml2::xml_contents(xml2::xml_find_all(html_raw, ".//latitude"))),
        xml2::xml_double(xml2::xml_contents(xml2::xml_find_all(html_raw, ".//longitude")))
      ))

      # Save stations by state in list format
      if (nrow(estac) > 0) {
        CodStat <- c(CodStat, list(estac))
        cat("Finished downloading data for", u, "\n")
        list(estac)
      }
    }
  }

  # Bind all states inventory together in single data frame
  if (length(CodStat) > 0) {
    CodStat <- do.call(rbind, CodStat)
    CodStat <- as.data.frame(CodStat) # Convertendo para data frame

    if (is.null(CodStat) || nrow(CodStat) == 0) {
      cat("No station data found.\n")
      return(NULL)
    }

    CodStat <- unique(CodStat) %>%
      data.frame() %>%
      stats::setNames(c("state", "codstation", "name", "lat", "long"))

    cat("\n")

    cat("Please wait, checking the availability of station data...\n")

    # Lista para armazenar os códigos das estações que atendem ao critério
    selected_stations <- list()

    # Lista para armazenar informações de latitude e longitude correspondentes às estações selecionadas
    selected_stations_info <- list()

    # Loop para verificar a disponibilidade dos dados das estações em paralelo
    foreach(i = 1:nrow(CodStat), .combine = "c") %do% {
      station_number <- gsub(" ", "%20", CodStat$codstation[i])

      cat("Station", station_number, "\n")

      html_raw <- tryCatch(
        {
          read_html(
            paste0(
              "http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroSerieHistorica?codEstacao=",
              station_number,
              "&dataInicio=&dataFim=&tipoDados=",
              stationType,
              "&nivelConsistencia="
            )
          )
        },
        error = function(e) {
          cat("Error downloading historical data for station", station_number, ". Moving to the next station.\n")
          return(NULL)
        }
      )

      if (!is.null(html_raw)) {
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
            as_tibble()
        }) %>%
          do.call(what = bind_rows)

        if ("maxima" %in% names(station_df)) {
          station_df <- station_df %>%
            filter(!is.na(maxima))

          anos <- lubridate::year(station_df$datahora) # Garanta que a função year() esteja disponível
          quantidade_anos <- length(unique(anos))

          # Supondo que 'quantidade_anos' é a quantidade de anos de dados que a estação possui
          if (quantidade_anos >= 20) {
            # Adicionar o código da estação à lista de estações selecionadas
            selected_stations <- c(selected_stations, station_number)
            selected_stations_info[[station_number]] <- list(State = CodStat$state[i], Name = CodStat$name[i], lat = CodStat$lat[i], long = CodStat$long[i])
            cat("Station", station_number, "has data for", quantidade_anos, "years. Adding station code to the list.\n")
          } else {
            cat("Data length less than 20 years for station", station_number, ". It has data for", quantidade_anos, "years.\n")
          }
        } else {
          cat("The 'maxima' column does not exist for station", station_number, ".\n")
        }
      }
    }

    # Desregistrar o cluster paralelo
    stopImplicitCluster()

    dfTab <- data.frame(
      codstation = unlist(selected_stations),
      State = sapply(selected_stations_info, function(x) x$State),
      Name = sapply(selected_stations_info, function(x) x$Name),
      lat = sapply(selected_stations_info, function(x) x$lat),
      long = sapply(selected_stations_info, function(x) x$long)
    )

    if (nrow(dfTab) == 0) {
      cat("No station data available.\n")
      return(NULL)
    } else {
      return(list(CodStat=dfTab))
    }
  } else {
    cat("No data was downloaded from any state.\n")
    return(NULL)
  }
}
