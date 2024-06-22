# Generates a table containing the main results
ScriptExit <- function(TamanhoArquivDuracoes,
                       TamanhoArquivTr,
                       best_result,
                       ArquivPrec,
                       IMaxObs,
                       dfTAD) {
  my_list05 <- list()

  TabFinal <- NULL
  a <- NA_real_
  b <- NA_real_
  c <- NA_real_
  d <- NA_real_
  e <- NA_real_
  f <- NA_real_
  g <- NA_real_
  h <- NA_real_

  IMaxSim <- matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)

  eq_number <- best_result$Equation

  if (eq_number == 1) {
    for (j in 1:TamanhoArquivTr) {
      for (i in 1:TamanhoArquivDuracoes) {
        IMaxSim[i, j] <- (best_result$aotim * (ArquivTr[j])^best_result$botim) / ((ArquivDuracoes[i] + best_result$cotim)^best_result$dotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))

    NS <- c(abs(round(min(
      hydroGOF::NSE(IMaxSim, IMaxObs)
    ), 3)))
    R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
    RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
    MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))
  } else if (eq_number == 2) {
    for (m in 1:TamanhoArquivTr) {
      for (n in 1:TamanhoArquivDuracoes) {
        IMaxSim[n, m] <- (best_result$aotim * (ArquivTr[m])^best_result$botim) / ((ArquivDuracoes[n] + best_result$cotim)^(best_result$dotim * (ArquivTr[m])^best_result$eotim))
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))

    NS <- c(abs(round(min(
      hydroGOF::NSE(IMaxSim, IMaxObs)
    ), 3)))
    R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
    RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
    MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))
  } else if (eq_number == 3) {
    for (s in 1:TamanhoArquivTr) {
      for (v in 1:TamanhoArquivDuracoes) {
        IMaxSim[v, s] <- ((best_result$aotim * (ArquivDuracoes[v] + best_result$botim)^best_result$cotim) + (best_result$dotim * (ArquivDuracoes[v] + best_result$eotim)^best_result$fotim) * (best_result$gotim + best_result$hotim * log(log(ArquivTr[s] / (ArquivTr[s] - 1))))) * 60
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

    NS <- c(abs(round(min(
      hydroGOF::NSE(IMaxSim, IMaxObs)
    ), 3)))
    R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
    RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
    MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))
  } else if (eq_number == 4) {
    for (b in 1:TamanhoArquivTr) {
      for (a in 1:TamanhoArquivDuracoes) {
        IMaxSim[a, b] <- (best_result$aotim * (ArquivTr[b] - best_result$botim)^best_result$cotim) / (((best_result$dotim * ArquivDuracoes[a]) + best_result$eotim)^best_result$fotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))
    f <- c(round(best_result$fotim, 3))

    NS <- c(abs(round(min(
      hydroGOF::NSE(IMaxSim, IMaxObs)
    ), 3)))
    R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
    RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
    MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))
  } else if (eq_number == 5) {
    for (g in 1:TamanhoArquivTr) {
      for (r in 1:TamanhoArquivDuracoes) {
        IMaxSim[r, g] <- (best_result$aotim * (ArquivTr[g]+best_result$botim)^best_result$cotim) / ((ArquivDuracoes[r] + best_result$dotim)^best_result$eotim)
      }
    }

    a <- c(round(best_result$aotim, 2))
    b <- c(round(best_result$botim, 3))
    c <- c(round(best_result$cotim, 2))
    d <- c(round(best_result$dotim, 3))
    e <- c(round(best_result$eotim, 3))

    NS <- c(abs(round(min(
      hydroGOF::NSE(IMaxSim, IMaxObs)
    ), 3)))
    R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
    RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
    MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))
  }

  # Criar TabFinal, definindo NA explicitamente onde necessário para a, b, c, d, e, f, g, h
  TabFinal <- data.frame(
    a = ifelse(all(is.na(a)), NA, a),
    b = ifelse(all(is.na(b)), NA, b),
    c = ifelse(all(is.na(c)), NA, c),
    d = ifelse(all(is.na(d)), NA, d),
    e = ifelse(all(is.na(e)), NA, e),
    f = ifelse(all(is.na(f)), NA, f),
    g = ifelse(all(is.na(g)), NA, g),
    h = ifelse(all(is.na(h)), NA, h),
    Distribuicao = c(dfTAD),
    NS = ifelse(is.na(NS), NA, NS),
    R2 = ifelse(is.na(R2), NA, R2),
    RMSE = ifelse(is.na(RMSE), NA, RMSE),
    MAE = ifelse(is.na(MAE), NA, MAE)
  )

  my_list05 <- list(TabFinal = TabFinal)

  return(my_list05)
}
