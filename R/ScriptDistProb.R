#Determine the maximum observed intensities from the probabilistic distribution that best fit the historical series
ScriptDistProb <- function(dfTAD,
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
                           ArquivDuracoes) {
  my_list02 <- list()

  result <- switch(
    dfTAD,
    "Normal" = cat(X <-
                     c(
                       stats::qnorm(
                         p = 1 - Probab,
                         mean = modnormal$estimate[1],
                         sd = modnormal$estimate[2],
                         lower.tail = TRUE,
                         log.p = FALSE
                       )
                     )),
    "LogNormal" = cat(X <-
                        c(
                          stats::qlnorm(
                            p = 1 - Probab,
                            meanlog = modln$estimate[1],
                            sdlog = modln$estimate[2],
                            lower.tail = TRUE,
                            log.p = FALSE
                          )
                        )),
    "Gamma" = cat(X <-
                    c(
                      stats::qgamma(
                        p = 1 - Probab,
                        shape = modgamma$estimate[1],
                        rate = modgamma$estimate[2],
                        lower.tail = TRUE,
                        log.p = FALSE
                      )
                    )),
    "Exponencial" = cat(X <-
                          c(
                            stats::qexp(p = 1 - Probab, rate = modexpo$estimate)
                          )),
    "Weibull" = cat(X <-
                      c(
                        stats::qweibull(1 - Probab,
                                        modweibull$estimate[1],
                                        modweibull$estimate[2])
                      )),
    "Gumbel" = cat(X <-
                     c(
                       evd::qgev(1 - Probab, modgumbel$estimate[1], modgumbel$estimate[2])
                     )),
    "GEV" = cat(X <-
                  c(
                    evd::qgev(
                      1 - Probab,
                      modGEV$estimate[1],
                      modGEV$estimate[2],
                      modGEV$estimate[3]
                    )
                  )),
    "Pearson" = cat(X <-
                      c(
                        PearsonDS::qpearson0(1 - Probab,
                                             mean(ArquivPrec),
                                             stats::sd(ArquivPrec))
                      )),
    "PearsonIII" = cat(X <-
                         c(
                           smwrBase::qpearsonIII(
                             1 - Probab,
                             mean(ArquivPrec),
                             stats::sd(ArquivPrec),
                             e1071::skewness(ArquivPrec, type =
                                               1)
                           )
                         )),
    "LogPearsonIII" = cat(X <-
                            c(
                              smwrBase::qlpearsonIII(
                                1 - Probab,
                                mean(log(ArquivPrec)),
                                stats::sd(log(ArquivPrec)),
                                e1071::skewness(log(ArquivPrec), type =
                                                  1)
                              )
                            ))
  )

  my_list02 <- list(X = X)

  return(my_list02)

}
