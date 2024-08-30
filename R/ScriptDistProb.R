# Determine the maximum observed intensities from the probabilistic distribution that best fit the historical series
ScriptDistProb <- function(dfTAD,
                           Probab,
                           modnormal,
                           modln,
                           modgamma,
                           modexpo,
                           modweibull,
                           modgumbel,
                           modGEV) {
  
  X <- NULL

  result <- switch(dfTAD,
    "Normal" = cat(X <-
      c(
        stats::qnorm(
          p = 1 - Probab,
          mean = modnormal$estimate[1],
          sd = modnormal$estimate[2],
          lower.tail = TRUE
        )
      )),
    "LogNormal" = cat(X <-
      c(
        stats::qlnorm(
          p = 1 - Probab,
          meanlog = modln$estimate[1],
          sdlog = modln$estimate[2],
          lower.tail = TRUE
        )
      )),
    "Gamma" = cat(X <-
      c(
        stats::qgamma(
          p = 1 - Probab,
          shape = modgamma$estimate[1],
          rate = modgamma$estimate[2],
          lower.tail = TRUE
        )
      )),
    "Exponencial" = cat(X <-
      c(
        stats::qexp(
          p = 1 - Probab,
          rate = modexpo$estimate,
          lower.tail = TRUE
        )
      )),
    "Weibull" = cat(X <-
      c(
        stats::qweibull(
          p = 1 - Probab,
          shape = modweibull$estimate[1],
          scale = modweibull$estimate[2],
          lower.tail = TRUE
        )
      )),
    "Gumbel" = cat(X <-
      c(
        evd::qgumbel(
          p = 1 - Probab,
          loc = modgumbel$estimate[1],
          scale = modgumbel$estimate[2],
          lower.tail = TRUE
        )
      )),
    "GEV" = cat(
      X <-
        c(evd::qgev(
          p = 1 - Probab,
          loc = modGEV$estimate[1],
          scale = modGEV$estimate[2],
          shape = modGEV$estimate[3],
          lower.tail = TRUE
        ))
    )
  )  

  return(X)
}
