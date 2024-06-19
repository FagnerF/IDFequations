# Optimizes the coefficients of the IDF equations
ScriptOtimiza <- function(TamanhoArquivDuracoes,
                          TamanhoArquivTr,
                          ArquivTr,
                          ArquivDuracoes,
                          IMaxObs) {
  my_list04 <- list()

  kotim <- NULL
  motim <- NULL
  t0otim <- NULL
  notim <- NULL

  IMaxSim <- matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)
  erro <- matrix(0, TamanhoArquivDuracoes, TamanhoArquivTr)
  Sum.erroInicial <- matrix(0, 1, TamanhoArquivTr)

  function.min <- function(par) {
    for (i in 1:TamanhoArquivTr) {
      for (j in 1:TamanhoArquivDuracoes) {
        IMaxSim[j, i] <-
          (par[1] * (ArquivTr[i])^par[2]) / ((ArquivDuracoes[j] + par[3])^par[4])
        erro[j, i] <-
          abs(IMaxSim[j, i] - IMaxObs[j, i]) / IMaxObs[j, i]
      }
      Sum.erroInicial[1, i] <- sum(erro[, i])
    }
    sum(Sum.erroInicial[1, ])
  }

  optmin <-
    stats::nlminb(
      c(0, 0, 0, 0),
      function.min,
      control = list(
        trace = FALSE,
        iter.max =
          100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf)
    )

  kinicial <- optmin$par[1]
  minicial <- optmin$par[2]
  t0inicial <- optmin$par[3]
  ninicial <- optmin$par[4]

  erro1 <-
    matrix(0, TamanhoArquivDuracoes * TamanhoArquivTr, 1)
  erro2 <-
    matrix(0, TamanhoArquivDuracoes * TamanhoArquivTr, 1)

  iter <- 0

  function.NS <- function(par) {
    for (k in 1:TamanhoArquivTr) {
      for (l in 1:TamanhoArquivDuracoes) {
        IMaxSim[l, k] <-
          (par[1] * (ArquivTr[k])^par[2]) / ((ArquivDuracoes[l] + par[3])^par[4])
        erro1[l + iter, 1] <-
          (IMaxSim[l, k] - IMaxObs[l, k])^
          2
        erro2[l + iter, 1] <-
          (IMaxObs[l, k] - mean(IMaxObs[, k]))^2
      }
      iter <- iter + 1
    }
    (1 - (sum(erro1) / sum(erro2)))
  }

  function.max <- function(par) {
    -function.NS(par)
  }

  optmax <-
    stats::nlminb(
      c(kinicial, minicial, t0inicial, ninicial),
      function.max,
      control = list(
        trace = FALSE,
        iter.max =
          100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf)
    )

  kotim <- optmax$par[1]
  motim <- optmax$par[2]
  t0otim <- optmax$par[3]
  notim <- optmax$par[4]

  my_list04 <-
    list(
      kotim = kotim,
      motim = motim,
      t0otim = t0otim,
      notim = notim
    )

  return(my_list04)
}
