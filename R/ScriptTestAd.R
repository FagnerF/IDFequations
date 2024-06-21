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
    tryCatch({
      stats::pnorm(ArquivPrec, modnormal$estimate[1], modnormal$estimate[2])
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modln <- fit_model("lnorm", ArquivPrec)
  dstln <- if (!is.null(modln)) {
    tryCatch({
      stats::plnorm(ArquivPrec, mean(log(ArquivPrec)), sd(log(ArquivPrec)))
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgamma <- fit_model("gamma", ArquivPrec)
  dstgamma <- if (!is.null(modgamma)) {
    tryCatch({
      stats::pgamma(ArquivPrec, modgamma$estimate[1], modgamma$estimate[2])
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modexpo <- fit_model("exp", ArquivPrec)
  dstexpo <- if (!is.null(modexpo)) {
    tryCatch({
      stats::pexp(ArquivPrec, modexpo$estimate)
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modweibull <- fit_model("weibull", ArquivPrec)
  dstweibull <- if (!is.null(modweibull)) {
    tryCatch({
      stats::pweibull(ArquivPrec, modweibull$estimate[1], modweibull$estimate[2])
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgumbel <- fit_gev(ArquivPrec, shape = 0)
  dstgumbel <- if (!is.null(modgumbel)) {
    tryCatch({
      evd::pgev(ArquivPrec, modgumbel$estimate[1], modgumbel$estimate[2])
    }, error = function(e) rep(NA, length(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modGEV <- fit_gev(ArquivPrec)
  dstGEV <- if (!is.null(modGEV)) {
    tryCatch({
      evd::pgev(ArquivPrec, modGEV$estimate[1], modGEV$estimate[2], modGEV$estimate[3])
    }, error = function(e) rep(NA, length(ArquivPrec)))
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

  Tab_KS <- read.delim("Kolmogorov-Smirnov.txt",sep = ";") %>%
    data.frame() %>%
    setNames(c("n", "20%", "10%", "5%", "2%", "1%"))

  Tab_Qui <- read.delim("Qui-Quadrado.txt",sep = ";") %>%
    data.frame() %>%
    setNames(c("n","0.005", "0.01", "0.025", "0.05", "0.1", "0.25", "0.5", "0.75", "0.9", "0.95", "0.975", "0.99", "0.995"))

  # Função para obter valor crítico do teste qui-quadrado
  get_critical_value <- function(n, alpha) {
    # Converter df para inteiro
    df <- as.integer(n)

    # Verificar se n está na tabela
    if (df %in% Tab_Qui$n) {
      val_critico <- Tab_Qui[Tab_Qui$n == df, as.character(alpha)]
    } else {
      #Fazer interpolação linear
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

  # Verificar se há empate entre duas ou mais distribuições
  if (nrow(dfTA) > 1 || min(dfTA$KS_Resultados) > valcrit) {
    # Função para calcular log-verossimilhança
    calculate_log_likelihood <- function(fit, distribution, data) {
      tryCatch({
        switch(distribution,
               norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
               lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
               gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
               exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
               weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
               gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
               GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
        )
      }, error = function(e) {
        -Inf  # Retorna -Inf em caso de erro
      })
    }

    # Calcular as log-verossimilhanças para todos os modelos
    log_likelihoods <- c(
      calculate_log_likelihood(modnormal, "norm", ArquivPrec),
      calculate_log_likelihood(modln, "lnorm", ArquivPrec),
      calculate_log_likelihood(modgamma, "gamma", ArquivPrec),
      calculate_log_likelihood(modexpo, "exp", ArquivPrec),
      calculate_log_likelihood(modweibull, "weibull", ArquivPrec),
      calculate_log_likelihood(modgumbel, "gumbel", ArquivPrec),
      calculate_log_likelihood(modGEV, "GEV", ArquivPrec)
    )

    # Adicionar log-verossimilhanças ao dataframe dfTA
    dfTA$log_likelihoods <- log_likelihoods

    # Selecionar a melhor distribuição pela log-verossimilhança
    dfTA <- dfTA %>%
      dplyr::filter(log_likelihoods == max(log_likelihoods, na.rm = TRUE))

    # Se houver empate, calcular p-valor para desempatar
    if (nrow(dfTA) > 1) {
      log_likelihood_null <- max(log_likelihoods, na.rm = TRUE)  # Log-verossimilhança do modelo nulo
      log_likelihood <- dfTA$log_likelihoods
      lambda <- -2 * (log_likelihood - log_likelihood_null)
      p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)

      n <- length(ArquivPrec)
      alpha <- 0.05 #nível de significância

      valcrit_Qui <- get_critical_value(n, alpha)

      # Comparar p-valor
      if (p_value < valcrit_Qui) {
        # Escolher a distribuição com menor KS_Resultados se p-valor for significativo
        best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
      } else {
        # Escolher a distribuição com maior log-verossimilhança se p-valor não for significativo
        best_distribution <- dfTA[which.max(dfTA$log_likelihoods), ]
      }
    } else {
      # Selecionar diretamente a distribuição única se não houver empate
      best_distribution <- dfTA
    }
  } else {
    # Selecionar a distribuição com menor KS_Resultados se não houver empate
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
