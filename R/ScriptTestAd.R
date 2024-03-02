#Identifies the probability distribution that best fits the historical series
ScriptTestAd <- function(ArquivPrec) {
  my_list01 <- list()
  Resultados <- list()

  Fr <- matrix(0, length(ArquivPrec), 1)

  for (j in 1:length(ArquivPrec)) {
    Fr[j] <- j / (length(ArquivPrec) + 1)

  }

  modnormal <- fitdist(ArquivPrec, "norm", method = "mle")
  dstnormal <-
    stats::pnorm(ArquivPrec,
                 modnormal$estimate[1],
                 modnormal$estimate[2])

  modln <- fitdist(ArquivPrec, "lnorm")
  dstln <-
    stats::plnorm(ArquivPrec,
                  mean(log(ArquivPrec)),
                  sd(log(ArquivPrec)))

  modgamma <- fitdist(ArquivPrec, "gamma")
  dstgamma <-
    stats::pgamma(ArquivPrec,
                  modgamma$estimate[1],
                  modgamma$estimate[2])

  modexpo <- fitdist(ArquivPrec, "exp")
  dstexpo <-
    stats::pexp(ArquivPrec,
                modexpo$estimate)

  modweibull <- fitdist(ArquivPrec, "weibull")
  dstweibull <-
    stats::pweibull(ArquivPrec,
                    modweibull$estimate[1],
                    modweibull$estimate[2])

  modgumbel <- evd::fgev(ArquivPrec, shape = 0)
  dstgumbel <-
    evd::pgev(ArquivPrec,
              modgumbel$estimate[1],
              modgumbel$estimate[2])

  modGEV <- evd::fgev(ArquivPrec)
  dstGEV <-
    evd::pgev(ArquivPrec,
              modGEV$estimate[1],
              modGEV$estimate[2],
              modGEV$estimate[3])

  dstpearson <-
    PearsonDS::ppearson0(ArquivPrec,
                         mean(ArquivPrec),
                         stats::sd(ArquivPrec))
  dstpearsonIII <-
    smwrBase::ppearsonIII(
      ArquivPrec,
      mean(ArquivPrec),
      stats::sd(ArquivPrec),
      e1071::skewness(ArquivPrec, type =
                        1)
    )
  dstlpearsonIII <-
    smwrBase::plpearsonIII(ArquivPrec,
                           mean(log(ArquivPrec)),
                           stats::sd(log(ArquivPrec)),
                           e1071::skewness(log(ArquivPrec), type =
                                             1))

  fxnormal <- matrix(0, length(ArquivPrec), 1)
  fxln <- matrix(0, length(ArquivPrec), 1)
  fxgamma <- matrix(0, length(ArquivPrec), 1)
  fxexpo <- matrix(0, length(ArquivPrec), 1)
  fxweibull <- matrix(0, length(ArquivPrec), 1)
  fxgumbel <- matrix(0, length(ArquivPrec), 1)
  fxGEV <- matrix(0, length(ArquivPrec), 1)
  fxpearson <- matrix(0, length(ArquivPrec), 1)
  fxpearsonIII <- matrix(0, length(ArquivPrec), 1)
  fxlpearsonIII <- matrix(0, length(ArquivPrec), 1)

  for (l in 1:length(ArquivPrec)) {
    fxnormal[l] <- abs(dstnormal[l] - Fr[l])
    fxln[l] <- abs(dstln[l] - Fr[l])
    fxgamma[l] <- abs(dstgamma[l] - Fr[l])
    fxexpo[l] <- abs(dstexpo[l] - Fr[l])
    fxweibull[l] <- abs(dstweibull[l] - Fr[l])
    fxgumbel[l] <- abs(dstgumbel[l] - Fr[l])
    fxGEV[l] <- abs(dstGEV[l] - Fr[l])
    fxpearson[l] <- abs(dstpearson[l] - Fr[l])
    fxpearsonIII[l] <- abs(dstpearsonIII[l] - Fr[l])
    fxlpearsonIII[l] <- abs(dstlpearsonIII[l] - Fr[l])

  }

  ksfxnormal <- max(fxnormal)
  ksfxln <- max(fxln)
  ksfxgamma <- max(fxgamma)
  ksfxexpo <- max(fxexpo)
  ksweibull <- max(fxweibull)
  ksGumbel <- max(fxgumbel)
  ksGEV <- max(fxGEV)
  ksPearson <- max(fxpearson)
  ksPearsonIII <- max(fxpearsonIII)
  ksLogPearsonIII <- max(fxlpearsonIII)

  dfTestAdr <- data.frame(
    Distribuicoes = c(
      "Normal",
      "LogNormal",
      "Gamma",
      "Exponencial",
      "Weibull",
      "Gumbel",
      "GEV",
      "Pearson",
      "PearsonIII",
      "LogPearsonIII"
    ),
    Resultados = c(
      round(ksfxnormal, 4),
      round(ksfxln, 4),
      round(ksfxgamma, 4),
      round(ksfxexpo, 4),
      round(ksweibull, 4),
      round(ksGumbel, 4),
      round(ksGEV, 4),
      round(ksPearson, 4),
      round(ksPearsonIII, 4),
      round(ksLogPearsonIII, 4)
    )
  )

  dfTA <- dfTestAdr %>%
    dplyr::filter(Resultados == min(Resultados))

  my_list01 <-
    list(
      dfTA = dfTA,
      modnormal = modnormal,
      modln = modln,
      modgamma = modgamma,
      modexpo = modexpo ,
      modweibull = modweibull,
      modgumbel = modgumbel,
      modGEV = modGEV,
      modPerason = dstpearson,
      modPersoanIII = dstpearsonIII,
      modLogPearsonIII = dstlpearsonIII
    )

  return(my_list01)

}
