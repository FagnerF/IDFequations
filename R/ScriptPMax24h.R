ScriptPMax24h <- function(X,
                          ArquivDuracoes,
                          ArquivTr,
                          Method,
                          Isozone) {
  IMaxObs <- matrix(0, length(ArquivDuracoes), length(ArquivTr))

  if (Method == "Isozone") {
    # Table KS
    dados <- matrix(c(
      0.362, 0.358, 0.356, 0.355, 0.354, 0.353, 0.350, 0.347, 0.336, 0.336, 0.07, 0.063,
      0.381, 0.378, 0.375, 0.374, 0.373, 0.372, 0.369, 0.366, 0.354, 0.343, 0.084, 0.075,
      0.401, 0.397, 0.395, 0.393, 0.392, 0.391, 0.388, 0.384, 0.372, 0.360, 0.098, 0.088,
      0.420, 0.416, 0.414, 0.412, 0.411, 0.410, 0.407, 0.403, 0.390, 0.378, 0.112, 0.10,
      0.440, 0.436, 0.433, 0.432, 0.430, 0.429, 0.426, 0.422, 0.409, 0.396, 0.126, 0.112,
      0.460, 0.455, 0.453, 0.451, 0.449, 0.448, 0.445, 0.441, 0.427, 0.413, 0.139, 0.124,
      0.479, 0.474, 0.472, 0.470, 0.468, 0.467, 0.464, 0.459, 0.445, 0.431, 0.154, 0.137,
      0.499, 0.494, 0.491, 0.489, 0.486, 0.486, 0.483, 0.478, 0.463, 0.448, 0.167, 0.149
    ), nrow = 8, byrow = TRUE)

    # Criando o dataframe
    df <- as.data.frame(dados) %>%
      setNames(c("5", "10", "15", "20", "25", "30", "50", "100", "1000", "10000", "Tr_5-50", "Tr_100"))

    # Inserindo a coluna das isozonas
    df$Isozone <- c("A", "B", "C", "D", "E", "F", "G", "H")

    # Reorganizar as colunas
    df <- df[, c(13, 1:12)]

    PMax24h <- matrix(0, 1, length(ArquivTr))
    PMax1h <- matrix(0, 1, length(ArquivTr))
    PMax6min <- matrix(0, 1, length(ArquivTr))

    for (w in 1:length(ArquivTr)) {
      PMax24h[w] <- X[w] * 1.095
      if (w < length(ArquivTr)) {
        PMax6min[w] <- PMax24h[w] * df[df$Isozone == Isozone, as.character("Tr_5-50")]
      } else {
        PMax6min[w] <- PMax24h[w] * df[df$Isozone == Isozone, as.character("Tr_100")]
      }

      PMax1h[w] <- PMax24h[w] * df[df$Isozone == Isozone, as.character(ArquivTr[w])]
    }

    # Dados fornecidos na forma de matriz
    matriz <- matrix(c(
      PMax6min,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      PMax1h,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA, NA, NA,
      PMax24h
    ), nrow = 10, byrow = TRUE)

    # Nomes das colunas e linhas
    colnames(matriz) <- ArquivTr
    rownames(matriz) <- ArquivDuracoes

    # Função para interpolar valores NA usando a interpolação linear
    interpolar_na <- function(matriz) {
      for (i in 1:ncol(matriz)) {
        # Encontrar índices dos valores não NA
        indices <- which(!is.na(matriz[, i]))
        if (length(indices) > 1) {
          # Interpolação linear usando a função approx
          matriz[, i] <- approx(rownames(matriz)[indices], matriz[indices, i], xout = rownames(matriz), method = "linear")$y
        }
      }
      return(matriz)
    }

    # Chamar a função para interpolar os valores NA na matriz
    matriz_interp <- interpolar_na(matriz)

    PMax <- round(matriz_interp, 1)
  } else if (Method == "CETESB") {
    PMax <- matrix(0, length(ArquivDuracoes) - 1, length(ArquivTr))
    PMax24h <- matrix(0, 1, length(ArquivTr))

    for (b in 1:length(ArquivTr)) {
      PMax24h[b] <- X[b] * 1.14

      for (f in 1:(length(ArquivDuracoes) - 1)) {
        PMax[f, b] <-
          PMax24h[b] * (exp(1.5 * log((
            log(ArquivDuracoes[f]) / 7.3
          ))))
      }
    }
  } else if (Method == "Back") {
    PMax <- matrix(0, length(ArquivDuracoes) - 1, length(ArquivTr))
    PMax24h <- matrix(0, 1, length(ArquivTr))

    for (b in 1:length(ArquivTr)) {
      PMax24h[b] <- X[b] * 1.14

      for (f in 1:(length(ArquivDuracoes) - 1)) {
        PMax[f, b] <-
          (ArquivDuracoes[f] / (27.9327 + (3.8346 * (ArquivDuracoes[f])^0.7924))) * PMax24h[b]
      }
    }
  }

  PMaxCorrig <-
    matrix(0, length(ArquivDuracoes), length(ArquivTr))

  for (l in 1:length(ArquivTr)) {
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

  for (f in 1:length(ArquivTr)) {
    for (g in 1:length(ArquivDuracoes)) {
      IMaxObs[g, f] <- (PMaxCorrig[g, f] / ArquivDuracoes[g]) * 60
    }
  }

  return(IMaxObs)
}
