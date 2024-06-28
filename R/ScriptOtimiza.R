# Define function to run optimization for each equation and return the max NS value and coefficients
ScriptOtimiza <- function(eq_number, ArquivTr, ArquivDuracoes, IMaxObs) {

  set.seed(123)

  # Initialize variables to store optimal parameters
  aotim <- NA_real_
  botim <- NA_real_
  cotim <- NA_real_
  dotim <- NA_real_
  eotim <- NA_real_
  fotim <- NA_real_

  IMaxSim <- matrix(0, length(ArquivDuracoes), length(ArquivTr))
  erro <- matrix(0, length(ArquivDuracoes), length(ArquivTr))
  Sum.erroInicial <- matrix(0, 1, length(ArquivTr))

  # Define function.min and function.NS based on the equation number
  if (eq_number == 1) {
    function.min <- function(par) {
      for (i in 1:length(ArquivTr)) {
        for (j in 1:length(ArquivDuracoes)) {
          IMaxSim[j, i] <- (par[1] * (ArquivTr[i])^par[2]) / ((ArquivDuracoes[j] + par[3])^par[4])
          erro[j, i] <- abs(IMaxSim[j, i] - IMaxObs[j, i]) / IMaxObs[j, i]
        }
        Sum.erroInicial[1, i] <- sum(erro[, i])
      }
      sum(Sum.erroInicial[1, ])
    }

    optmin <- stats::nlminb(
      c(0, 0, 0, 0),
      function.min,
      control = list(
        trace = FALSE,
        iter.max = 100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf)
    )

    ainicial <- optmin$par[1]
    binicial <- optmin$par[2]
    cinicial <- optmin$par[3]
    dinicial <- optmin$par[4]

    erro1 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)
    erro2 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)

    iter <- 0

    function.NS <- function(par) {
      for (k in 1:length(ArquivTr)) {
        for (l in 1:length(ArquivDuracoes)) {
          IMaxSim[l, k] <- (par[1] * (ArquivTr[k])^par[2]) / ((ArquivDuracoes[l] + par[3])^par[4])
          erro1[l + iter, 1] <- (IMaxSim[l, k] - IMaxObs[l, k])^2
          erro2[l + iter, 1] <- (IMaxObs[l, k] - mean(IMaxObs[, k]))^2
        }
        iter <- iter + 1
      }
      (1 - (sum(erro1) / sum(erro2)))
    }

    function.max <- function(par) {
      -function.NS(par)
    }

    optmax <- stats::nlminb(
      c(ainicial, binicial, cinicial, dinicial),
      function.max,
      control = list(
        trace = FALSE,
        iter.max = 100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf)
    )

    aotim <- optmax$par[1]
    botim <- optmax$par[2]
    cotim <- optmax$par[3]
    dotim <- optmax$par[4]
  } else if (eq_number == 2) {
    function.min <- function(par) {
      for (m in 1:length(ArquivTr)) {
        for (n in 1:length(ArquivDuracoes)) {
          IMaxSim[n, m] <- (par[1] * (ArquivTr[m])^par[2]) / ((ArquivDuracoes[n] + par[3])^(par[4] * (ArquivTr[m])^par[5]))
          erro[n, m] <- abs(IMaxSim[n, m] - IMaxObs[n, m]) / IMaxObs[n, m]
        }
        Sum.erroInicial[1, m] <- sum(erro[, m])
      }
      sum(Sum.erroInicial[1, ])
    }

    optmin <- stats::nlminb(
      c(0, 0, 0, 0, 0),
      function.min,
      control = list(
        trace = FALSE,
        iter.max = 100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf, Inf)
    )

    ainicial <- optmin$par[1]
    binicial <- optmin$par[2]
    cinicial <- optmin$par[3]
    dinicial <- optmin$par[4]
    einicial <- optmin$par[5]

    erro1 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)
    erro2 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)

    iter <- 0

    function.NS <- function(par) {
      for (o in 1:length(ArquivTr)) {
        for (p in 1:length(ArquivDuracoes)) {
          IMaxSim[p, o] <- (par[1] * (ArquivTr[o])^par[2]) / ((ArquivDuracoes[p] + par[3])^(par[4] * (ArquivTr[o])^par[5]))
          erro1[p + iter, 1] <- (IMaxSim[p, o] - IMaxObs[p, o])^2
          erro2[p + iter, 1] <- (IMaxObs[p, o] - mean(IMaxObs[, o]))^2
        }
        iter <- iter + 1
      }
      (1 - (sum(erro1) / sum(erro2)))
    }

    function.max <- function(par) {
      -function.NS(par)
    }

    optmax <- stats::nlminb(
      c(ainicial, binicial, cinicial, dinicial, einicial),
      function.max,
      control = list(
        trace = FALSE,
        iter.max = 100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf, Inf)
    )

    aotim <- optmax$par[1]
    botim <- optmax$par[2]
    cotim <- optmax$par[3]
    dotim <- optmax$par[4]
    eotim <- optmax$par[5]
  } else if (eq_number == 3) {
    function.min <- function(par) {
      for (b in 1:length(ArquivTr)) {
        for (a in 1:length(ArquivDuracoes)) {
          IMaxSim[a, b] <-
            (par[1] * (ArquivTr[b] - par[2])^par[3]) / (((par[4] * ArquivDuracoes[a]) + par[5])^par[6])
          erro[a, b] <-
            abs(IMaxSim[a, b] - IMaxObs[a, b]) / IMaxObs[a, b]
        }
        Sum.erroInicial[1, b] <- sum(erro[, b])
      }
      sum(Sum.erroInicial[1, ])
    }

    optmin <-
      stats::nlminb(
        c(0, 0, 0, 0, 0, 0),
        function.min,
        control = list(
          trace = FALSE,
          iter.max =
            100000,
          eval.max = 20000
        ),
        lower = c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
        upper = c(Inf, Inf, Inf, Inf, Inf, Inf)
      )

    ainicial <- optmin$par[1]
    binicial <- optmin$par[2]
    cinicial <- optmin$par[3]
    dinicial <- optmin$par[4]
    einicial <- optmin$par[5]
    finicial <- optmin$par[6]

    erro1 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)
    erro2 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)

    iter <- 0

    function.NS <- function(par) {
      for (x in 1:length(ArquivTr)) {
        for (z in 1:length(ArquivDuracoes)) {
          IMaxSim[z, x] <-
            (par[1] * (ArquivTr[x] - par[2])^par[3]) / (((par[4] * ArquivDuracoes[z]) + par[5])^par[6])
          erro1[z + iter, 1] <-
            (IMaxSim[z, x] - IMaxObs[z, x])^
            2
          erro2[z + iter, 1] <-
            (IMaxObs[z, x] - mean(IMaxObs[, x]))^2
        }
        iter <- iter + 1
      }
      (1 - (sum(erro1) / sum(erro2)))
    }

    function.max <- function(par) {
      -function.NS(par)
    }

    optmax <- stats::nlminb(
      c(ainicial, binicial, cinicial, dinicial, einicial, finicial),
      function.max,
      control = list(
        trace = FALSE,
        iter.max = 100000,
        eval.max = 20000
      ),
      lower = c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
      upper = c(Inf, Inf, Inf, Inf, Inf, Inf)
    )

    aotim <- optmax$par[1]
    botim <- optmax$par[2]
    cotim <- optmax$par[3]
    dotim <- optmax$par[4]
    eotim <- optmax$par[5]
    fotim <- optmax$par[6]
  } else if (eq_number == 4) {
    function.min <- function(par) {
      for (g in 1:length(ArquivTr)) {
        for (r in 1:length(ArquivDuracoes)) {
          IMaxSim[r, g] <-
            (par[1] * (ArquivTr[g] + par[2])^par[3]) / ((ArquivDuracoes[r] + par[4])^par[5])
          erro[r, g] <-
            abs(IMaxSim[r, g] - IMaxObs[r, g]) / IMaxObs[r, g]
        }
        Sum.erroInicial[1, g] <- sum(erro[, g])
      }
      sum(Sum.erroInicial[1, ])
    }

    optmin <-
      stats::nlminb(
        c(0, 0, 0, 0, 0),
        function.min,
        control = list(
          trace = FALSE,
          iter.max =
            100000,
          eval.max = 20000
        ),
        lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
        upper = c(Inf, Inf, Inf, Inf, Inf)
      )

    ainicial <- optmin$par[1]
    binicial <- optmin$par[2]
    cinicial <- optmin$par[3]
    dinicial <- optmin$par[4]
    einicial <- optmin$par[5]

    erro1 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)
    erro2 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)

    iter <- 0

    function.NS <- function(par) {
      for (u in 1:length(ArquivTr)) {
        for (z in 1:length(ArquivDuracoes)) {
          IMaxSim[z, u] <-
            (par[1] * (ArquivTr[u] + par[2])^par[3]) / ((ArquivDuracoes[z] + par[4])^par[5])
          erro1[z + iter, 1] <-
            (IMaxSim[z, u] - IMaxObs[z, u])^
            2
          erro2[z + iter, 1] <-
            (IMaxObs[z, u] - mean(IMaxObs[, u]))^2
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
        c(ainicial, binicial, cinicial, dinicial, einicial),
        function.max,
        control = list(
          trace = FALSE,
          iter.max =
            100000,
          eval.max = 20000
        ),
        lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
        upper = c(Inf, Inf, Inf, Inf, Inf)
      )

    aotim <- optmax$par[1]
    botim <- optmax$par[2]
    cotim <- optmax$par[3]
    dotim <- optmax$par[4]
    eotim <- optmax$par[5]
  } else if (eq_number == 5) {
    function.min <- function(par) {
      for (c in 1:length(ArquivTr)) {
        for (d in 1:length(ArquivDuracoes)) {
          IMaxSim[d, c] <-
            (par[1] * (ArquivTr[c])^par[2]) / ((ArquivDuracoes[d])^par[3])
          erro[d, c] <-
            abs(IMaxSim[d, c] - IMaxObs[d, c]) / IMaxObs[d, c]
        }
        Sum.erroInicial[1, c] <- sum(erro[, c])
      }
      sum(Sum.erroInicial[1, ])
    }

    optmin <-
      stats::nlminb(
        c(0, 0, 0),
        function.min,
        control = list(
          trace = FALSE,
          iter.max =
            100000,
          eval.max = 20000
        ),
        lower = c(0.001, 0.001, 0.001),
        upper = c(Inf, Inf, Inf)
      )

    ainicial <- optmin$par[1]
    binicial <- optmin$par[2]
    cinicial <- optmin$par[3]

    erro1 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)
    erro2 <- matrix(0, length(ArquivDuracoes) * length(ArquivTr), 1)

    iter <- 0

    function.NS <- function(par) {
      for (c in 1:length(ArquivTr)) {
        for (d in 1:length(ArquivDuracoes)) {
          IMaxSim[d, c] <-
            (par[1] * (ArquivTr[c])^par[2]) / ((ArquivDuracoes[d])^par[3])
          erro1[d + iter, 1] <-
            (IMaxSim[d, c] - IMaxObs[d, c])^
            2
          erro2[d + iter, 1] <-
            (IMaxObs[d, c] - mean(IMaxObs[, c]))^2
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
        c(ainicial, binicial, cinicial),
        function.max,
        control = list(
          trace = FALSE,
          iter.max =
            100000,
          eval.max = 20000
        ),
        lower = c(0.001, 0.001, 0.001),
        upper = c(Inf, Inf, Inf)
      )

    aotim <- optmax$par[1]
    botim <- optmax$par[2]
    cotim <- optmax$par[3]
  } else {
    stop("Invalid equation number.")
  }

  # Calcular o valor de NS após otimização
  NS_value <- function.NS(c(aotim, botim, cotim, dotim, eotim, fotim))

  # Criar um row para o dataframe de resultados
  result <- data.frame(
    Equation = eq_number,
    aotim = aotim,
    botim = botim,
    cotim = cotim,
    dotim = dotim,
    eotim = eotim,
    fotim = fotim,
    NS = NS_value,
    stringsAsFactors = FALSE
  )

  return(result)
}
