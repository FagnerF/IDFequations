ScriptPMax24h <- function(X,
                          TamanhoArquivDuracoes,
                          TamanhoArquivTr,
                          ArquivDuracoes) {

  my_list03 <- list()

  PMax <- matrix(0, TamanhoArquivDuracoes - 1, TamanhoArquivTr)
  PMax24h <- matrix(0, 1, TamanhoArquivTr)

  for (b in 1:TamanhoArquivTr) {
    PMax24h[b] <- X[b] * 1.14

    for (f in 1:(TamanhoArquivDuracoes - 1)) {
      PMax[f, b] <-
        PMax24h[b] * (exp(1.5 * log((
          log(ArquivDuracoes[f]) / 7.3
        ))))

    }

  }

  PMaxCorrig <-
    matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)

  for (l in 1:TamanhoArquivTr) {
    PMaxCorrig[1, l] <- ifelse(PMax[1, l] < 8, 8, PMax[1, l])
    PMaxCorrig[2, l] <- ifelse(PMax[2, l] < 10, 10, PMax[2, l])
    PMaxCorrig[3, l] <- ifelse(PMax[3, l] < 15, 15, PMax[3, l])
    PMaxCorrig[4, l] <- ifelse(PMax[4, l] < 15, 15, PMax[4, l])
    PMaxCorrig[5, l] <- ifelse(PMax[5, l] < 20, 20, PMax[5, l])
    PMaxCorrig[6, l] <- ifelse(PMax[6, l] < 25, 25, PMax[6, l])
    PMaxCorrig[7, l] <- ifelse(PMax[7, l] < 40, 40, PMax[7, l])
    PMaxCorrig[8, l] <- ifelse(PMax[8, l] < 40, 40, PMax[8, l])
    PMaxCorrig[9, l] <- ifelse(PMax[9, l] < 47, 47, PMax[9, l])
    PMaxCorrig[10, l] <- ifelse(PMax24h[l] < 55, 55, PMax24h[l])

  }

  IMaxObs <- matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)

  for (f in 1:TamanhoArquivTr) {
    for (g in 1:TamanhoArquivDuracoes) {
      IMaxObs[g, f] <- (PMaxCorrig[g, f] / ArquivDuracoes[g]) * 60

    }

  }

  my_list03 <- list(IMaxObs = IMaxObs)

  return(my_list03)

}
