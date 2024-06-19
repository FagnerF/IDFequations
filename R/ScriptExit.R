# Generates a table containing the main results
ScriptExit <- function(TamanhoArquivDuracoes,
                       TamanhoArquivTr,
                       kotim,
                       ArquivTr,
                       motim,
                       ArquivDuracoes,
                       t0otim,
                       notim,
                       ArquivPrec,
                       IMaxObs,
                       dfTAD) {
  my_list05 <- list()

  IMaxSim <- matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)

  for (l in 1:TamanhoArquivTr) {
    for (w in 1:TamanhoArquivDuracoes) {
      IMaxSim[w, l] <-
        ((kotim) * (ArquivTr[l])^(motim)) / ((ArquivDuracoes[w] + (t0otim))^
          (notim))
    }
  }

  K <- c(round(kotim, 2))
  m <- c(round(motim, 3))
  t0 <- c(round(t0otim, 2))
  n <- c(round(notim, 3))

  Distribuicao <- c(dfTAD)

  NS <- c(abs(round(min(
    hydroGOF::NSE(IMaxSim, IMaxObs)
  ), 3)))
  R2 <- c(round(min(hydroGOF::br2(IMaxSim, IMaxObs)), 3))
  RMSE <- c(round(max(hydroGOF::rmse(IMaxSim, IMaxObs)), 3))
  MAE <- c(round(max(hydroGOF::mae(IMaxSim, IMaxObs)), 3))

  TabFinal <-
    data.frame(K, m, t0, n, Distribuicao, NS, R2, RMSE, MAE)

  my_list05 <- list(TabFinal = TabFinal)

  return(my_list05)
}
