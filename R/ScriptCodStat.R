ScriptCodStat <- function(StatesANA, stationType) {
  my_list <- list()
  CodStat <- list()
  estac <- list()
  codstation <- list()
  lat <- list()
  long <- list()
  serief <-  list()

  cat("Downloading... \n")

  for (u in 1:length(StatesANA)) {
    # Adjusting string for API query
    estadoG <- gsub(" ", "%20", StatesANA[u])

    # Raw HTML to retrieve stations inventory
    html_raw <-
      xml2::read_html(
        paste(
          "http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroInventario?codEstDE=&codEstATE=&tpEst=",
          stationType,
          "&nmEst=&nmRio=&codSubBacia=&codBacia=&nmMunicipio=&nmEstado=",
          estadoG,
          "&sgResp=&sgOper=&telemetrica=",
          sep = ""
        )
      )

    # Scrapping info for each station
    estac <- as.data.frame(cbind(
      xml2::xml_double(xml2::xml_contents(
        xml2::xml_find_all(html_raw, ".//codigo")
      )),
      xml2::xml_double(xml2::xml_contents(
        xml2::xml_find_all(html_raw, ".//latitude")
      )),
      xml2::xml_double(xml2::xml_contents(
        xml2::xml_find_all(html_raw, ".//longitude")
      ))
    ))

    # Save stations by state in list format
    serief[[u]] <- estac
    cat(StatesANA[u], " finished. \n")

  }

  # Bind all states inventory together in single data frame
  estac <- do.call(rbind, serief) %>%
    setNames(c("station_code", "lat", "long")) %>%
    data.frame()

  estac <- unique(estac)

  codstation <- matrix(0, nrow(estac), 1)
  lat <- matrix(0, nrow(estac), 1)
  long <- matrix(0, nrow(estac), 1)

  cat("\n")

  cat("Please wait, checking the availability of station data... \n")

  for (i in 1:nrow(estac)) {
    #Download the historical series from the ANA data-base
    station_number <- gsub(" ", "%20", estac$station_code[i])

    #Create list to receive results
    Data <- list()
    Ano <- list()
    Prec <- list()
    serie <- list()

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

    if (ncol(station_df) > 1) {
      #Data-base station?
      station_df <- station_df %>%
        # Convert precipitation columns to numeric and datahora to date format
        dplyr::mutate(
          dplyr::across(dplyr::matches("chuva..$"), as.numeric),
          dplyr::across(dplyr::matches("data"), as.Date)
        ) %>%
        dplyr::rename(data = dplyr::any_of("datahora")) %>%
        dplyr::arrange(dplyr::across('data'))

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

      ArquivPrec <- MaxPrec[row_sub, ]
      ArquivPrecNumeric <- as.numeric(MaxPrec[row_sub, 2])
      ArquivPrecCrescente <- sort(c(ArquivPrecNumeric))
      ArquivPrec <- ArquivPrecCrescente

      if (length(ArquivPrec) >= 30) {
        codstation[i, 1] <- estac$station_code[i]
        lat[i, 1] <- estac$lat[i]
        long[i, 1] <- estac$long[i]

      }

    }

  }

  CodStat <- data.frame(codstation, lat, long)

  CodStat <- CodStat[rowSums(CodStat[]) > 0,]

  my_list <- list(CodStat = CodStat)

  return(my_list)

}
