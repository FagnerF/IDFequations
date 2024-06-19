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
    stats::pnorm(ArquivPrec, modnormal$estimate[1], modnormal$estimate[2])
  } else {
    rep(NA, length(ArquivPrec))
  }

  modln <- fit_model("lnorm", ArquivPrec)
  dstln <- if (!is.null(modln)) {
    stats::plnorm(ArquivPrec, mean(log(ArquivPrec)), sd(log(ArquivPrec)))
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgamma <- fit_model("gamma", ArquivPrec)
  dstgamma <- if (!is.null(modgamma)) {
    stats::pgamma(ArquivPrec, modgamma$estimate[1], modgamma$estimate[2])
  } else {
    rep(NA, length(ArquivPrec))
  }

  modexpo <- fit_model("exp", ArquivPrec)
  dstexpo <- if (!is.null(modexpo)) {
    stats::pexp(ArquivPrec, modexpo$estimate)
  } else {
    rep(NA, length(ArquivPrec))
  }

  modweibull <- fit_model("weibull", ArquivPrec)
  dstweibull <- if (!is.null(modweibull)) {
    stats::pweibull(ArquivPrec, modweibull$estimate[1], modweibull$estimate[2])
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgumbel <- fit_gev(ArquivPrec, shape = 0)
  dstgumbel <- if (!is.null(modgumbel)) {
    evd::pgev(ArquivPrec, modgumbel$estimate[1], modgumbel$estimate[2])
  } else {
    rep(NA, length(ArquivPrec))
  }

  modGEV <- fit_gev(ArquivPrec)
  dstGEV <- if (!is.null(modGEV)) {
    evd::pgev(ArquivPrec, modGEV$estimate[1], modGEV$estimate[2], modGEV$estimate[3])
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

  Tab_KS <- data.frame(
    n = 1:40,
    "0.2" = c(
      0.9, 0.684, 0.565, 0.493, 0.447, 0.41, 0.381, 0.358, 0.339, 0.323,
      0.308, 0.296, 0.285, 0.275, 0.266, 0.258, 0.25, 0.244, 0.237, 0.232,
      0.226, 0.221, 0.216, 0.212, 0.208, 0.204, 0.2, 0.197, 0.193, 0.19,
      0.187, 0.184, 0.182, 0.179, 0.177, 0.174, 0.172, 0.17, 0.168, 0.165
    ),
    "0.1" = c(
      0.95, 0.776, 0.636, 0.565, 0.509, 0.468, 0.436, 0.41, 0.387, 0.369,
      0.352, 0.338, 0.325, 0.314, 0.304, 0.295, 0.286, 0.279, 0.271, 0.265,
      0.259, 0.253, 0.247, 0.242, 0.238, 0.233, 0.229, 0.225, 0.221, 0.218,
      0.214, 0.211, 0.208, 0.205, 0.202, 0.199, 0.196, 0.194, 0.191, 0.189
    ),
    "0.05" = c(
      0.975, 0.842, 0.708, 0.624, 0.563, 0.519, 0.483, 0.454, 0.43, 0.409,
      0.391, 0.375, 0.361, 0.349, 0.338, 0.327, 0.318, 0.309, 0.301, 0.294,
      0.287, 0.281, 0.275, 0.269, 0.264, 0.259, 0.254, 0.25, 0.246, 0.242,
      0.238, 0.234, 0.231, 0.227, 0.224, 0.221, 0.218, 0.215, 0.213, 0.21
    ),
    "0.02" = c(
      0.99, 0.9, 0.785, 0.689, 0.627, 0.577, 0.538, 0.507, 0.48, 0.457,
      0.437, 0.419, 0.404, 0.39, 0.377, 0.366, 0.355, 0.346, 0.337, 0.329,
      0.321, 0.314, 0.307, 0.301, 0.295, 0.29, 0.284, 0.279, 0.275, 0.27,
      0.266, 0.262, 0.258, 0.254, 0.251, 0.247, 0.244, 0.241, 0.238, 0.235
    ),
    "0.01" = c(
      0.995, 0.929, 0.829, 0.734, 0.669, 0.617, 0.576, 0.542, 0.513, 0.489,
      0.468, 0.449, 0.432, 0.418, 0.404, 0.392, 0.381, 0.371, 0.361, 0.352,
      0.344, 0.337, 0.33, 0.323, 0.317, 0.311, 0.305, 0.3, 0.295, 0.29,
      0.285, 0.181, 0.277, 0.273, 0.269, 0.265, 0.262, 0.258, 0.255, 0.252
    )
  ) %>%
    setNames(c("N. amostras", "20%", "10%", "5%", "2%", "1%"))

  # Valor crítico para o teste KS (alfa = 5%)
  valcrit <- ifelse(length(ArquivPrec) > 40, 1.36 / sqrt(length(ArquivPrec)), Tab_KS$`5%`[length(ArquivPrec)])

  # Verificar se há empate entre duas ou mais distribuições
  if (nrow(dfTA) > 1 || min(dfTA$KS_Resultados) > valcrit) {
    # Função para calcular log-verossimilhança
    calculate_log_likelihood <- function(fit, distribution, data) {
      switch(distribution,
        norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
        lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
        gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
        exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
        weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
        gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
        GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
      )
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
      log_likelihood <- dfTA$log_likelihoods
      lambda <- -2 * (log_likelihood - log_likelihood)
      p_value <- pchisq(lambda, df = 0, lower.tail = FALSE)

      # Comparar p-valor com nível de significância (ex.: 0.05)
      if (p_value < 0.05) {
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
