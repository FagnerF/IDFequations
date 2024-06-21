ScriptTestAd <- function(ArquivPrec) {
  my_list01 <- list()

  dfTAD <- NULL
  modnormal <- NULL
  modln <- NULL
  modgamma <- NULL
  modexpo <- NULL
  modweibull <- NULL
  modgumbel <- NULL
  modGEV <- NULL

  Fr <- seq(1, length(ArquivPrec)) / (length(ArquivPrec) + 1)

  fit_model <- function(distribution, data, ...) {
    tryCatch(
      {
        suppressWarnings(fitdistrplus::fitdist(data, distribution, ...))
      },
      error = function(e) {
        NULL
      }
    )
  }

  fit_gev <- function(data, ...) {
    tryCatch(
      {
        suppressWarnings(evd::fgev(data, ...))
      },
      warning = function(w) {
        NULL
      },
      error = function(e) {
        NULL
      }
    )
  }

  modnormal <- fit_model("norm", ArquivPrec, method = "mle")
  dstnormal <- if (!is.null(modnormal)) {
    tryCatch(
      {
        stats::pnorm(ArquivPrec, modnormal$estimate[1], modnormal$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modln <- fit_model("lnorm", ArquivPrec)
  dstln <- if (!is.null(modln)) {
    tryCatch(
      {
        stats::plnorm(ArquivPrec, mean(log(ArquivPrec)), sd(log(ArquivPrec)))
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgamma <- fit_model("gamma", ArquivPrec)
  dstgamma <- if (!is.null(modgamma)) {
    tryCatch(
      {
        stats::pgamma(ArquivPrec, modgamma$estimate[1], modgamma$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modexpo <- fit_model("exp", ArquivPrec)
  dstexpo <- if (!is.null(modexpo)) {
    tryCatch(
      {
        stats::pexp(ArquivPrec, modexpo$estimate)
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modweibull <- fit_model("weibull", ArquivPrec)
  dstweibull <- if (!is.null(modweibull)) {
    tryCatch(
      {
        stats::pweibull(ArquivPrec, modweibull$estimate[1], modweibull$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgumbel <- fit_gev(ArquivPrec, shape = 0)
  dstgumbel <- if (!is.null(modgumbel)) {
    tryCatch(
      {
        evd::pgev(ArquivPrec, modgumbel$estimate[1], modgumbel$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modGEV <- fit_gev(ArquivPrec)
  dstGEV <- if (!is.null(modGEV)) {
    tryCatch(
      {
        evd::pgev(ArquivPrec, modGEV$estimate[1], modGEV$estimate[2], modGEV$estimate[3])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  # Verificar se todas as variáveis de ajuste são NULL
  if (is.null(modnormal) && is.null(modln) && is.null(modgamma) &&
    is.null(modexpo) && is.null(modweibull) && is.null(modgumbel) && is.null(modGEV)) {
    return(NULL)
  }

  calculate_ks <- function(dst, Fr) {
    if (all(is.na(dst))) {
      return(Inf)
    }
    max(abs(dst - Fr))
  }

  # Calcular estatísticas KS para todas as distribuições ajustadas
  ksfxnormal <- calculate_ks(dstnormal, Fr)
  ksfxln <- calculate_ks(dstln, Fr)
  ksfxgamma <- calculate_ks(dstgamma, Fr)
  ksfxexpo <- calculate_ks(dstexpo, Fr)
  ksweibull <- calculate_ks(dstweibull, Fr)
  ksGumbel <- calculate_ks(dstgumbel, Fr)
  ksGEV <- calculate_ks(dstGEV, Fr)

  # Criar dataframe com resultados KS
  dfTA <- data.frame(
    Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
    KS_Resultados = c(
      round(ksfxnormal, 4), round(ksfxln, 4), round(ksfxgamma, 4),
      round(ksfxexpo, 4), round(ksweibull, 4), round(ksGumbel, 4),
      round(ksGEV, 4)
    )
  )

  Tab_KS <- read.delim("Kolmogorov-Smirnov.txt", sep = ";") %>%
    data.frame() %>%
    setNames(c("n", "20%", "10%", "5%", "2%", "1%"))

  Tab_Qui <- read.delim("Qui-Quadrado.txt", sep = ";") %>%
    data.frame() %>%
    setNames(c("n", "0.005", "0.01", "0.025", "0.05", "0.1", "0.25", "0.5", "0.75", "0.9", "0.95", "0.975", "0.99", "0.995"))

  # Função para obter valor crítico do teste qui-quadrado
  get_critical_value <- function(n, alpha) {
    # Converter df para inteiro
    df <- as.integer(n)

    # Verificar se n está na tabela
    if (df %in% Tab_Qui$n) {
      val_critico <- Tab_Qui[Tab_Qui$n == df, as.character(alpha)]
    } else {
      # Fazer interpolação linear
      # Encontrar os índices inferior e superior
      idx <- which(Tab_Qui$n >= df)[1]
      df_lower <- Tab_Qui$n[idx - 1]
      df_upper <- Tab_Qui$n[idx]

      # Coletar os valores para interpolação
      y_lower <- Tab_Qui[idx - 1, as.character(alpha)]
      y_upper <- Tab_Qui[idx, as.character(alpha)]

      # Interpolação linear entre df_lower e df_upper
      val_critico <- y_lower + (df - df_lower) * (y_upper - y_lower) / (df_upper - df_lower)
    }

    return(val_critico)
  }

  # Valor crítico para o teste KS (alfa = 5%)
  valcrit <- ifelse(length(ArquivPrec) > 40, 1.36 / sqrt(length(ArquivPrec)), Tab_KS$`5%`[length(ArquivPrec)])

  # Identificar as distribuições que empataram em KS_Resultados
  empates_KS <- dfTA[duplicated(dfTA$KS_Resultados) | duplicated(dfTA$KS_Resultados, fromLast = TRUE), ]

  # Verifica se há empate entre as distribuições
  if (nrow(empates_KS) > 1) {
    # Função para calcular a razão da verossimilhança com tratamento de erros
    calculate_likelihood_ratio <- function(fit, distribution, data, null_data) {
      tryCatch(
        {
          log_likelihood_fit <- switch(distribution,
            norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
            lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
            gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
            exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
            weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
          )

          log_likelihood_null <- switch(distribution,
            norm = sum(stats::dnorm(null_data, mean = mean(null_data), sd = sd(null_data), log = TRUE)),
            lnorm = sum(stats::dlnorm(null_data, meanlog = meanlog(null_data), sdlog = sdlog(null_data), log = TRUE)),
            gamma = sum(stats::dgamma(null_data, shape = 1, rate = mean(null_data), log = TRUE)),
            exp = sum(stats::dexp(null_data, rate = mean(null_data), log = TRUE)),
            weibull = sum(stats::dweibull(null_data, shape = 1, scale = mean(null_data), log = TRUE)),
            gumbel = sum(evd::dgumbel(null_data, loc = mean(null_data), scale = sd(null_data), log = TRUE)),
            GEV = sum(evd::dgev(null_data, loc = mean(null_data), scale = sd(null_data), shape = 0, log = TRUE))
          )

          lambda <- -2 * (log_likelihood_fit - log_likelihood_null)
          return(lambda)
        },
        error = function(e) {
          -Inf # Retorna -Inf em caso de erro
        }
      )
    }

    # Calcular a razão da verossimilhança para os modelos empatados
    likelihood_ratios <- c(
      if ("Normal" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modnormal, "norm", ArquivPrec, rnorm(length(ArquivPrec), mean = mean(ArquivPrec), sd = sd(ArquivPrec)))
      } else {
        -Inf
      },
      if ("LogNormal" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modln, "lnorm", ArquivPrec, rlnorm(length(ArquivPrec), meanlog = meanlog(ArquivPrec), sdlog = sdlog(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Gamma" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modgamma, "gamma", ArquivPrec, rgamma(length(ArquivPrec), shape = 1, rate = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Exponencial" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modexpo, "exp", ArquivPrec, rexp(length(ArquivPrec), rate = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Weibull" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modweibull, "weibull", ArquivPrec, rweibull(length(ArquivPrec), shape = 1, scale = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Gumbel" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modgumbel, "gumbel", ArquivPrec, rgumbel(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec)))
      } else {
        -Inf
      },
      if ("GEV" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modGEV, "GEV", ArquivPrec, rgev(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec), shape = 0))
      } else {
        -Inf
      }
    )

    # Criar dataframe com resultados das razões da verossimilhança
    dfTAV <- data.frame(
      Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
      Likelihood_Ratio = likelihood_ratios
    )

    # Selecionar a melhor distribuição pela log-verossimilhança
    dfTAV <- dfTAV %>%
      dplyr::filter(Likelihood_Ratio == max(Likelihood_Ratio, na.rm = TRUE))

    # Calcular o número de graus de liberdade para cada distribuição
    dfTAV <- dfTAV %>%
      dplyr::mutate(
        df = case_when(
          Distribuicoes == "Normal" ~ length(coef(modnormal)),
          Distribuicoes == "LogNormal" ~ length(coef(modln)),
          Distribuicoes == "Gamma" ~ length(coef(modgamma)),
          Distribuicoes == "Exponencial" ~ length(coef(modexpo)),
          Distribuicoes == "Weibull" ~ length(coef(modweibull)),
          Distribuicoes == "Gumbel" ~ length(coef(modgumbel)),
          Distribuicoes == "GEV" ~ length(coef(modGEV)),
          TRUE ~ NA_integer_ # Caso haja alguma distribuição não contemplada
        )
      )

    p_value <- stats::pchisq(dfTAV$Likelihood_Ratio[1], df = dfTAV$df[1], lower.tail = FALSE)

    n <- length(ArquivPrec)
    alpha <- 0.05 # nível de significância

    valcrit_Qui <- get_critical_value(n, alpha)

    # Comparar p-valor
    if (p_value < valcrit_Qui) {
      # Escolher a distribuição com menor KS_Resultados se p-valor for significativo
      best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
    } else {
      # Escolher a distribuição com maior log-verossimilhança se p-valor não for significativo
      best_distribution <- dfTAV[which.max(dfTAV$Likelihood_Ratio), ]
    }

    # Verifica se min(dfTA$KS_Resultados) > valcrit
  } else if (min(dfTA$KS_Resultados) > valcrit) {
    # Função para calcular a razão da verossimilhança com tratamento de erros
    calculate_likelihood_ratio <- function(fit, distribution, data, null_data) {
      tryCatch(
        {
          log_likelihood_fit <- switch(distribution,
            norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
            lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
            gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
            exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
            weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
          )

          log_likelihood_null <- switch(distribution,
            norm = sum(stats::dnorm(null_data, mean = mean(null_data), sd = sd(null_data), log = TRUE)),
            lnorm = sum(stats::dlnorm(null_data, meanlog = meanlog(null_data), sdlog = sdlog(null_data), log = TRUE)),
            gamma = sum(stats::dgamma(null_data, shape = 1, rate = mean(null_data), log = TRUE)),
            exp = sum(stats::dexp(null_data, rate = mean(null_data), log = TRUE)),
            weibull = sum(stats::dweibull(null_data, shape = 1, scale = mean(null_data), log = TRUE)),
            gumbel = sum(evd::dgumbel(null_data, loc = mean(null_data), scale = sd(null_data), log = TRUE)),
            GEV = sum(evd::dgev(null_data, loc = mean(null_data), scale = sd(null_data), shape = 0, log = TRUE))
          )

          lambda <- -2 * (log_likelihood_fit - log_likelihood_null)
          return(lambda)
        },
        error = function(e) {
          -Inf # Retorna -Inf em caso de erro
        }
      )
    }

    # Calcular a razão da verossimilhança para todos os modelos
    likelihood_ratios <- c(
      calculate_likelihood_ratio(modnormal, "norm", ArquivPrec, rnorm(length(ArquivPrec), mean = mean(ArquivPrec), sd = sd(ArquivPrec))),
      calculate_likelihood_ratio(modln, "lnorm", ArquivPrec, rlnorm(length(ArquivPrec), meanlog = meanlog(ArquivPrec), sdlog = sdlog(ArquivPrec))),
      calculate_likelihood_ratio(modgamma, "gamma", ArquivPrec, rgamma(length(ArquivPrec), shape = 1, rate = mean(ArquivPrec))),
      calculate_likelihood_ratio(modexpo, "exp", ArquivPrec, rexp(length(ArquivPrec), rate = mean(ArquivPrec))),
      calculate_likelihood_ratio(modweibull, "weibull", ArquivPrec, rweibull(length(ArquivPrec), shape = 1, scale = mean(ArquivPrec))),
      calculate_likelihood_ratio(modgumbel, "gumbel", ArquivPrec, rgumbel(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec))),
      calculate_likelihood_ratio(modGEV, "GEV", ArquivPrec, rgev(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec), shape = 0))
    )

    # Criar dataframe com resultados das razões da verossimilhança
    dfTAV <- data.frame(
      Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
      Likelihood_Ratio = likelihood_ratios
    )

    # Selecionar a melhor distribuição pela log-verossimilhança
    dfTAV <- dfTAV %>%
      dplyr::filter(Likelihood_Ratio == max(Likelihood_Ratio, na.rm = TRUE))

    # Calcular o número de graus de liberdade para cada distribuição
    dfTAV <- dfTAV %>%
      dplyr::mutate(
        df = case_when(
          Distribuicoes == "Normal" ~ length(coef(modnormal)),
          Distribuicoes == "LogNormal" ~ length(coef(modln)),
          Distribuicoes == "Gamma" ~ length(coef(modgamma)),
          Distribuicoes == "Exponencial" ~ length(coef(modexpo)),
          Distribuicoes == "Weibull" ~ length(coef(modweibull)),
          Distribuicoes == "Gumbel" ~ length(coef(modgumbel)),
          Distribuicoes == "GEV" ~ length(coef(modGEV)),
          TRUE ~ NA_integer_ # Caso haja alguma distribuição não contemplada
        )
      )

    p_value <- stats::pchisq(dfTAV$Likelihood_Ratio[1], df = dfTAV$df[1], lower.tail = FALSE)

    n <- length(ArquivPrec)
    alpha <- 0.05 # nível de significância

    valcrit_Qui <- get_critical_value(n, alpha)

    # Comparar p-valor
    if (p_value < valcrit_Qui) {
      # Escolher a distribuição com menor KS_Resultados se p-valor for significativo
      best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
    } else {
      # Escolher a distribuição com maior log-verossimilhança se p-valor não for significativo
      best_distribution <- dfTAV[which.max(dfTAV$Likelihood_Ratio), ]
    }
  } else {
    # Selecionar a distribuição com menor KS_Resultados
    best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
  }

  my_list01 <- list(
    dfTAD = best_distribution,
    modnormal = modnormal,
    modln = modln,
    modgamma = modgamma,
    modexpo = modexpo,
    modweibull = modweibull,
    modgumbel = modgumbel,
    modGEV = modGEV
  )

  return(my_list01)
}
