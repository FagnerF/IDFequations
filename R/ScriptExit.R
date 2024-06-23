# Generates a table containing the main results
ScriptExit <- function(ArquivDuracoes,
                       ArquivTr,
                       best_result,
                       IMaxObs,
                       dfTAD) {
  # Cria o dataframe TabFinal
  TabFinal <- matrix(NA, nrow = 1, ncol = 13)

  IMaxSim <- matrix(0, nrow = length(ArquivDuracoes), ncol = length(ArquivTr))

  eq <- best_result$Equation

  if (eq == 1) {
    for (j in 1:length(ArquivTr)) {
      for (i in 1:length(ArquivDuracoes)) {
        IMaxSim[i, j] <- (best_result$aotim * (ArquivTr[j])^best_result$botim) / ((ArquivDuracoes[i] + best_result$cotim)^best_result$dotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
  } else if (eq == 2) {
    for (j in 1:length(ArquivTr)) {
      for (i in 1:length(ArquivDuracoes)) {
        IMaxSim[i, j] <- (best_result$aotim * (ArquivTr[j])^best_result$botim) / ((ArquivDuracoes[i] + best_result$cotim)^(best_result$dotim * (ArquivTr[j])^best_result$eotim))
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))
  } else if (eq == 3) {
    for (j in 1:length(ArquivTr)) {
      for (i in 1:length(ArquivDuracoes)) {
        IMaxSim[i, j] <- ((best_result$aotim * (ArquivDuracoes[j] + best_result$botim)^best_result$cotim) + (best_result$dotim * (ArquivDuracoes[j] + best_result$eotim)^best_result$fotim) * (best_result$gotim + best_result$hotim * log(log(ArquivTr[i] / (ArquivTr[i] - 1))))) * 60
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 2))
    f <- c(round(best_result$fotim, 3))
    g <- c(round(best_result$gotim, 2))
    h <- c(round(best_result$hotim, 3))
  } else if (eq == 4) {
    for (j in 1:length(ArquivTr)) {
      for (i in 1:length(ArquivDuracoes)) {
        IMaxSim[i, j] <- (best_result$aotim * (ArquivTr[j] - best_result$botim)^best_result$cotim) / (((best_result$dotim * ArquivDuracoes[i]) + best_result$eotim)^best_result$fotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))
    f <- c(round(best_result$fotim, 3))
  } else if (eq == 5) {
    for (j in 1:length(ArquivTr)) {
      for (i in 1:length(ArquivDuracoes)) {
        IMaxSim[i, j] <- (best_result$aotim * (ArquivTr[j] + best_result$botim)^best_result$cotim) / ((ArquivDuracoes[i] + best_result$dotim)^best_result$eotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))
  }

  NS <- c(abs(round(min(
    hydroGOF::NSE(IMaxSim, IMaxObs)
  ), 3)))
  R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
  RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
  MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))

  # Cria o dataframe TabFinal, definindo NA explicitamente
  TabFinal <- data.frame(
    a = a,
    b = b,
    c = c,
    d = d,
    e = e,
    f = f,
    g = g,
    h = h,
    Distribuicao = dfTAD,
    NS = NS,
    R2 = R2,
    RMSE = RMSE,
    MAE = MAE,
    stringsAsFactors = FALSE
  )

  return(TabFinal)
}
