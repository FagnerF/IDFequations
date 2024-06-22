ScriptTestAd <- function(ArquivPrec) {
  my_list01 <- list()

  dfTAD <- NULL
  modnormal <- NULL
  modln <- NULL
  modgamma <- NULL
  modexpo <- NULL
  modweibull <- NULL
  modgumbel <- NULL
  modGEV <- NULL

  Fr <- seq(1, length(ArquivPrec)) / (length(ArquivPrec) + 1)

  fit_model <- function(distribution, data, ...) {
    tryCatch(
      {
        suppressWarnings(fitdistrplus::fitdist(data, distribution, ...))
      },
      error = function(e) {
        NULL
      }
    )
  }

  fit_gev <- function(data, ...) {
    tryCatch(
      {
        suppressWarnings(evd::fgev(data, ...))
      },
      warning = function(w) {
        NULL
      },
      error = function(e) {
        NULL
      }
    )
  }

  modnormal <- fit_model("norm", ArquivPrec, method = "mle")
  dstnormal <- if (!is.null(modnormal)) {
    tryCatch(
      {
        stats::pnorm(ArquivPrec, modnormal$estimate[1], modnormal$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modln <- fit_model("lnorm", ArquivPrec)
  dstln <- if (!is.null(modln)) {
    tryCatch(
      {
        stats::plnorm(ArquivPrec, mean(log(ArquivPrec)), sd(log(ArquivPrec)))
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgamma <- fit_model("gamma", ArquivPrec)
  dstgamma <- if (!is.null(modgamma)) {
    tryCatch(
      {
        stats::pgamma(ArquivPrec, modgamma$estimate[1], modgamma$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modexpo <- fit_model("exp", ArquivPrec)
  dstexpo <- if (!is.null(modexpo)) {
    tryCatch(
      {
        stats::pexp(ArquivPrec, modexpo$estimate)
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modweibull <- fit_model("weibull", ArquivPrec)
  dstweibull <- if (!is.null(modweibull)) {
    tryCatch(
      {
        stats::pweibull(ArquivPrec, modweibull$estimate[1], modweibull$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modgumbel <- fit_gev(ArquivPrec, shape = 0)
  dstgumbel <- if (!is.null(modgumbel)) {
    tryCatch(
      {
        evd::pgev(ArquivPrec, modgumbel$estimate[1], modgumbel$estimate[2])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  modGEV <- fit_gev(ArquivPrec)
  dstGEV <- if (!is.null(modGEV)) {
    tryCatch(
      {
        evd::pgev(ArquivPrec, modGEV$estimate[1], modGEV$estimate[2], modGEV$estimate[3])
      },
      error = function(e) rep(NA, length(ArquivPrec))
    )
  } else {
    rep(NA, length(ArquivPrec))
  }

  # Verificar se todas as variáveis de ajuste são NULL
  if (is.null(modnormal) && is.null(modln) && is.null(modgamma) &&
    is.null(modexpo) && is.null(modweibull) && is.null(modgumbel) && is.null(modGEV)) {
    return(NULL)
  }

  calculate_ks <- function(dst, Fr) {
    if (all(is.na(dst))) {
      return(Inf)
    }
    max(abs(dst - Fr))
  }

  # Calcular estatísticas KS para todas as distribuições ajustadas
  ksfxnormal <- calculate_ks(dstnormal, Fr)
  ksfxln <- calculate_ks(dstln, Fr)
  ksfxgamma <- calculate_ks(dstgamma, Fr)
  ksfxexpo <- calculate_ks(dstexpo, Fr)
  ksweibull <- calculate_ks(dstweibull, Fr)
  ksGumbel <- calculate_ks(dstgumbel, Fr)
  ksGEV <- calculate_ks(dstGEV, Fr)

  # Criar dataframe com resultados KS
  dfTA <- data.frame(
    Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
    KS_Resultados = c(
      round(ksfxnormal, 4), round(ksfxln, 4), round(ksfxgamma, 4),
      round(ksfxexpo, 4), round(ksweibull, 4), round(ksGumbel, 4),
      round(ksGEV, 4)
    )
  )

  #Valores tabelado para o teste KS
  dados_KS <- "
n;0.2;0.1;0.05;0.02;0.01
1;0.9;0.95;0.975;0.99;0.995
2;0.684;0.776;0.842;0.9;0.929
3;0.565;0.636;0.708;0.785;0.829
4;0.493;0.565;0.624;0.689;0.734
5;0.447;0.509;0.563;0.627;0.669
6;0.41;0.468;0.519;0.577;0.617
7;0.381;0.436;0.483;0.538;0.576
8;0.358;0.41;0.454;0.407;0.542
9;0.339;0.387;0.43;0.48;0.513
10;0.323;0.369;0.409;0.457;0.489
11;0.308;0.352;0.391;0.437;0.468
12;0.296;0.338;0.375;0.419;0.449
13;0.285;0.325;0.361;0.404;0.432
14;0.275;0.314;0.349;0.39;0.418
15;0.266;0.304;0.338;0.377;0.404
16;0.258;0.295;0.327;0.366;0.392
17;0.25;0.286;0.318;0.355;0.381
18;0.244;0.279;0.309;0.346;0.371
19;0.237;0.271;0.301;0.337;0.361
20;0.232;0.265;0.294;0.329;0.352
21;0.226;0.259;0.287;0.321;0.344
22;0.221;0.253;0.281;0.314;0.337
23;0.216;0.247;0.275;0.307;0.33
24;0.212;0.242;0.269;0.301;0.323
25;0.208;0.238;0.264;0.295;0.317
26;0.204;0.233;0.259;0.29;0.311
27;0.2;0.229;0.254;0.284;0.305
28;0.197;0.225;0.25;0.279;0.3
29;0.193;0.221;0.246;0.275;0.295
30;0.19;0.218;0.242;0.27;0.29
31;0.187;0.214;0.238;0.266;0.285
32;0.184;0.211;0.234;0.262;0.181
33;0.182;0.208;0.231;0.258;0.277
34;0.179;0.205;0.227;0.254;0.273
35;0.177;0.202;0.224;0.251;0.269
36;0.174;0.199;0.221;0.247;0.265
37;0.172;0.196;0.218;0.244;0.262
38;0.17;0.194;0.215;0.241;0.258
39;0.168;0.191;0.213;0.238;0.255
40;0.165;0.189;0.21;0.235;0.252
"

  # Convertendo para dataframe
  Tab_KS <- read.table(text = dados_KS, sep = ";", header = TRUE) %>%
    data.frame() %>%
    setNames(c("n", "20%", "10%", "5%", "2%", "1%"))

  #Valores tabelado para o teste Qui-quadrado
  dados_Qui <- "
n;0.005;0.01;0.025;0.05;0.1;0.25;0.5;0.75;0.9;0.95;0.975;0.99;0.995
1;3.93E-05;0.000157;0.000982;0.003932;0.016;0.102;0.455;1.323;2.706;3.841;5.024;6.635;7.879
2;0.010;0.020;0.051;0.103;0.211;0.575;1.386;2.773;4.605;5.991;7.378;9.210;10.597
3;0.072;0.115;0.216;0.352;0.584;1.213;2.366;4.108;6.251;7.815;9.348;11.345;12.838
4;0.207;0.297;0.484;0.711;1.064;1.923;3.357;5.385;7.779;9.488;11.143;13.277;14.860
5;0.412;0.554;0.831;1.145;1.610;2.675;4.351;6.626;9.236;11.070;12.832;15.086;16.750
6;0.676;0.872;1.237;1.635;2.204;3.455;5.348;7.841;10.645;12.592;14.449;16.812;18.548
7;0.989;1.239;1.690;2.167;2.833;4.255;6.346;9.037;12.017;14.067;16.013;18.475;20.278
8;1.344;1.647;2.180;2.733;3.490;5.071;7.344;10.219;13.362;15.507;17.535;20.090;21.955
9;1.735;2.088;2.700;3.325;4.168;5.899;8.343;11.389;14.684;16.919;19.023;21.666;23.589
10;2.156;2.558;3.247;3.940;4.865;6.737;9.342;12.549;15.987;18.307;20.483;23.209;25.188
11;2.603;3.053;3.816;4.575;5.578;7.584;10.341;13.701;17.275;19.675;21.920;24.725;26.757
12;3.074;3.571;4.404;5.226;6.304;8.438;11.340;14.845;18.549;21.026;23.337;26.217;28.300
13;3.565;4.107;5.009;5.892;7.041;9.299;12.340;15.984;19.812;22.362;24.736;27.688;29.819
14;4.075;4.660;5.629;6.571;7.790;10.165;13.339;17.117;21.064;23.685;26.119;29.141;31.319
15;4.601;5.229;6.262;7.261;8.547;11.037;14.339;18.245;22.307;24.996;27.488;30.578;32.801
16;5.142;5.812;6.908;7.962;9.312;11.912;15.338;19.369;23.542;26.296;28.845;32.000;34.267
17;5.697;6.408;7.564;8.672;10.085;12.792;16.338;20.489;24.769;27.587;30.191;33.409;35.718
18;6.265;7.015;8.231;9.390;10.865;13.675;17.338;21.605;25.989;28.869;31.526;34.805;37.156
19;6.844;7.633;8.907;10.117;11.651;14.562;18.338;22.718;27.204;30.144;32.852;36.191;38.582
20;7.434;8.260;9.591;10.851;12.443;15.452;19.337;23.828;28.412;31.410;34.170;37.566;39.997
21;8.034;8.897;10.283;11.591;13.240;16.344;20.337;24.935;29.615;32.671;35.479;38.932;41.401
22;8.643;9.542;10.982;12.338;14.041;17.240;21.337;26.039;30.813;33.924;36.781;40.289;42.796
23;9.260;10.196;11.689;13.091;14.848;18.137;22.337;27.141;32.007;35.172;38.076;41.638;44.181
24;9.886;10.856;12.401;13.848;15.659;19.037;23.337;28.241;33.196;36.415;39.364;42.980;45.558
25;10.520;11.524;13.120;14.611;16.473;19.939;24.337;29.339;34.382;37.652;40.646;44.314;46.928
26;11.160;12.198;13.844;15.379;17.292;20.843;25.336;30.435;35.563;38.885;41.923;45.642;48.290
27;11.808;12.878;14.573;16.151;18.114;21.749;26.336;31.528;36.741;40.113;43.195;46.963;49.645
28;12.461;13.565;15.308;16.928;18.939;22.657;27.336;32.620;37.916;41.337;44.461;48.278;50.994
29;13.121;14.256;16.047;17.708;19.768;23.567;28.336;33.711;39.087;42.557;45.722;49.588;52.335
30;13.787;14.953;16.791;18.493;20.599;24.478;29.336;34.800;40.256;43.773;46.979;50.892;53.672
40;20.707;22.164;24.433;26.509;29.051;33.660;39.335;45.616;51.805;55.758;59.342;63.691;66.766
50;27.991;29.707;32.357;34.764;37.689;42.942;49.335;56.334;63.167;67.505;71.420;76.154;79.490
60;35.534;37.485;40.482;43.188;46.459;52.294;59.335;66.981;74.397;79.082;83.298;88.379;91.952
70;43.275;45.442;48.758;51.739;55.329;61.698;69.334;77.577;85.527;90.531;95.023;100.425;104.215
80;51.172;53.540;57.153;60.391;64.278;71.145;79.334;88.130;96.578;101.879;106.629;112.329;116.321
90;59.196;61.754;65.647;69.126;73.291;80.625;89.334;98.650;107.565;113.145;118.136;124.116;128.299
100;67.328;70.065;74.222;77.929;82.358;90.133;99.334;109.141;118.498;124.342;129.561;135.807;140.170
"

  # Convertendo para dataframe
  Tab_Qui <- read.table(text = dados_Qui, sep = ";", header = TRUE) %>%
    data.frame() %>%
    setNames(c("n", "0.005", "0.01", "0.025", "0.05", "0.1", "0.25", "0.5", "0.75", "0.9", "0.95", "0.975", "0.99", "0.995"))

  # Função para obter valor crítico do teste qui-quadrado
  get_critical_value <- function(n, alpha) {
    # Converter df para inteiro
    df <- as.integer(n)

    # Verificar se n está na tabela
    if (df %in% Tab_Qui$n) {
      val_critico <- Tab_Qui[Tab_Qui$n == df, as.character(alpha)]
    } else {
      # Fazer interpolação linear
      # Encontrar os índices inferior e superior
      idx <- which(Tab_Qui$n >= df)[1]
      df_lower <- Tab_Qui$n[idx - 1]
      df_upper <- Tab_Qui$n[idx]

      # Coletar os valores para interpolação
      y_lower <- Tab_Qui[idx - 1, as.character(alpha)]
      y_upper <- Tab_Qui[idx, as.character(alpha)]

      # Interpolação linear entre df_lower e df_upper
      val_critico <- y_lower + (df - df_lower) * (y_upper - y_lower) / (df_upper - df_lower)
    }

    return(val_critico)
  }

  # Valor crítico para o teste KS (alfa = 5%)
  valcrit <- ifelse(length(ArquivPrec) > 40, 1.36 / sqrt(length(ArquivPrec)), Tab_KS$`5%`[length(ArquivPrec)])

  # Identificar as distribuições que empataram em KS_Resultados
  empates_KS <- dfTA[duplicated(dfTA$KS_Resultados) | duplicated(dfTA$KS_Resultados, fromLast = TRUE), ]

  # Verifica se há empate entre as distribuições
  if (nrow(empates_KS) > 1) {
    # Função para calcular a razão da verossimilhança com tratamento de erros
    calculate_likelihood_ratio <- function(fit, distribution, data, null_data) {
      tryCatch(
        {
          log_likelihood_fit <- switch(distribution,
            norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
            lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
            gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
            exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
            weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
          )

          log_likelihood_null <- switch(distribution,
            norm = sum(stats::dnorm(null_data, mean = mean(null_data), sd = sd(null_data), log = TRUE)),
            lnorm = sum(stats::dlnorm(null_data, meanlog = meanlog(null_data), sdlog = sdlog(null_data), log = TRUE)),
            gamma = sum(stats::dgamma(null_data, shape = 1, rate = mean(null_data), log = TRUE)),
            exp = sum(stats::dexp(null_data, rate = mean(null_data), log = TRUE)),
            weibull = sum(stats::dweibull(null_data, shape = 1, scale = mean(null_data), log = TRUE)),
            gumbel = sum(evd::dgumbel(null_data, loc = mean(null_data), scale = sd(null_data), log = TRUE)),
            GEV = sum(evd::dgev(null_data, loc = mean(null_data), scale = sd(null_data), shape = 0, log = TRUE))
          )

          lambda <- -2 * (log_likelihood_fit - log_likelihood_null)
          return(lambda)
        },
        error = function(e) {
          -Inf # Retorna -Inf em caso de erro
        }
      )
    }

    # Calcular a razão da verossimilhança para os modelos empatados
    likelihood_ratios <- c(
      if ("Normal" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modnormal, "norm", ArquivPrec, rnorm(length(ArquivPrec), mean = mean(ArquivPrec), sd = sd(ArquivPrec)))
      } else {
        -Inf
      },
      if ("LogNormal" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modln, "lnorm", ArquivPrec, rlnorm(length(ArquivPrec), meanlog = meanlog(ArquivPrec), sdlog = sdlog(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Gamma" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modgamma, "gamma", ArquivPrec, rgamma(length(ArquivPrec), shape = 1, rate = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Exponencial" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modexpo, "exp", ArquivPrec, rexp(length(ArquivPrec), rate = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Weibull" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modweibull, "weibull", ArquivPrec, rweibull(length(ArquivPrec), shape = 1, scale = mean(ArquivPrec)))
      } else {
        -Inf
      },
      if ("Gumbel" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modgumbel, "gumbel", ArquivPrec, rgumbel(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec)))
      } else {
        -Inf
      },
      if ("GEV" %in% empates_KS$Distribuicoes) {
        calculate_likelihood_ratio(modGEV, "GEV", ArquivPrec, rgev(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec), shape = 0))
      } else {
        -Inf
      }
    )

    # Criar dataframe com resultados das razões da verossimilhança
    dfTAV <- data.frame(
      Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
      Likelihood_Ratio = likelihood_ratios
    )

    # Selecionar a melhor distribuição pela log-verossimilhança
    dfTAV <- dfTAV %>%
      dplyr::filter(Likelihood_Ratio == max(Likelihood_Ratio, na.rm = TRUE))

    # Calcular o número de graus de liberdade para cada distribuição
    dfTAV <- dfTAV %>%
      dplyr::mutate(
        df = case_when(
          Distribuicoes == "Normal" ~ length(coef(modnormal)),
          Distribuicoes == "LogNormal" ~ length(coef(modln)),
          Distribuicoes == "Gamma" ~ length(coef(modgamma)),
          Distribuicoes == "Exponencial" ~ length(coef(modexpo)),
          Distribuicoes == "Weibull" ~ length(coef(modweibull)),
          Distribuicoes == "Gumbel" ~ length(coef(modgumbel)),
          Distribuicoes == "GEV" ~ length(coef(modGEV)),
          TRUE ~ NA_integer_ # Caso haja alguma distribuição não contemplada
        )
      )

    p_value <- stats::pchisq(dfTAV$Likelihood_Ratio[1], df = dfTAV$df[1], lower.tail = FALSE)

    n <- length(ArquivPrec)
    alpha <- 0.05 # nível de significância

    valcrit_Qui <- get_critical_value(n, alpha)

    # Comparar p-valor
    if (p_value < valcrit_Qui) {
      # Escolher a distribuição com menor KS_Resultados se p-valor for significativo
      best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
    } else {
      # Escolher a distribuição com maior log-verossimilhança se p-valor não for significativo
      best_distribution <- dfTAV[which.max(dfTAV$Likelihood_Ratio), ]
    }

    # Verifica se min(dfTA$KS_Resultados) > valcrit
  } else if (min(dfTA$KS_Resultados) > valcrit) {
    # Função para calcular a razão da verossimilhança com tratamento de erros
    calculate_likelihood_ratio <- function(fit, distribution, data, null_data) {
      tryCatch(
        {
          log_likelihood_fit <- switch(distribution,
            norm = sum(stats::dnorm(data, mean = fit$estimate[1], sd = fit$estimate[2], log = TRUE)),
            lnorm = sum(stats::dlnorm(data, meanlog = fit$estimate[1], sdlog = fit$estimate[2], log = TRUE)),
            gamma = sum(stats::dgamma(data, shape = fit$estimate[1], rate = fit$estimate[2], log = TRUE)),
            exp = sum(stats::dexp(data, rate = fit$estimate, log = TRUE)),
            weibull = sum(stats::dweibull(data, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            gumbel = sum(evd::dgumbel(data, loc = fit$estimate[1], scale = fit$estimate[2], log = TRUE)),
            GEV = sum(evd::dgev(data, loc = fit$estimate[1], scale = fit$estimate[2], shape = fit$estimate[3], log = TRUE))
          )

          log_likelihood_null <- switch(distribution,
            norm = sum(stats::dnorm(null_data, mean = mean(null_data), sd = sd(null_data), log = TRUE)),
            lnorm = sum(stats::dlnorm(null_data, meanlog = meanlog(null_data), sdlog = sdlog(null_data), log = TRUE)),
            gamma = sum(stats::dgamma(null_data, shape = 1, rate = mean(null_data), log = TRUE)),
            exp = sum(stats::dexp(null_data, rate = mean(null_data), log = TRUE)),
            weibull = sum(stats::dweibull(null_data, shape = 1, scale = mean(null_data), log = TRUE)),
            gumbel = sum(evd::dgumbel(null_data, loc = mean(null_data), scale = sd(null_data), log = TRUE)),
            GEV = sum(evd::dgev(null_data, loc = mean(null_data), scale = sd(null_data), shape = 0, log = TRUE))
          )

          lambda <- -2 * (log_likelihood_fit - log_likelihood_null)
          return(lambda)
        },
        error = function(e) {
          -Inf # Retorna -Inf em caso de erro
        }
      )
    }

    # Calcular a razão da verossimilhança para todos os modelos
    likelihood_ratios <- c(
      calculate_likelihood_ratio(modnormal, "norm", ArquivPrec, rnorm(length(ArquivPrec), mean = mean(ArquivPrec), sd = sd(ArquivPrec))),
      calculate_likelihood_ratio(modln, "lnorm", ArquivPrec, rlnorm(length(ArquivPrec), meanlog = meanlog(ArquivPrec), sdlog = sdlog(ArquivPrec))),
      calculate_likelihood_ratio(modgamma, "gamma", ArquivPrec, rgamma(length(ArquivPrec), shape = 1, rate = mean(ArquivPrec))),
      calculate_likelihood_ratio(modexpo, "exp", ArquivPrec, rexp(length(ArquivPrec), rate = mean(ArquivPrec))),
      calculate_likelihood_ratio(modweibull, "weibull", ArquivPrec, rweibull(length(ArquivPrec), shape = 1, scale = mean(ArquivPrec))),
      calculate_likelihood_ratio(modgumbel, "gumbel", ArquivPrec, rgumbel(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec))),
      calculate_likelihood_ratio(modGEV, "GEV", ArquivPrec, rgev(length(ArquivPrec), loc = mean(ArquivPrec), scale = sd(ArquivPrec), shape = 0))
    )

    # Criar dataframe com resultados das razões da verossimilhança
    dfTAV <- data.frame(
      Distribuicoes = c("Normal", "LogNormal", "Gamma", "Exponencial", "Weibull", "Gumbel", "GEV"),
      Likelihood_Ratio = likelihood_ratios
    )

    # Selecionar a melhor distribuição pela log-verossimilhança
    dfTAV <- dfTAV %>%
      dplyr::filter(Likelihood_Ratio == max(Likelihood_Ratio, na.rm = TRUE))

    # Calcular o número de graus de liberdade para cada distribuição
    dfTAV <- dfTAV %>%
      dplyr::mutate(
        df = case_when(
          Distribuicoes == "Normal" ~ length(coef(modnormal)),
          Distribuicoes == "LogNormal" ~ length(coef(modln)),
          Distribuicoes == "Gamma" ~ length(coef(modgamma)),
          Distribuicoes == "Exponencial" ~ length(coef(modexpo)),
          Distribuicoes == "Weibull" ~ length(coef(modweibull)),
          Distribuicoes == "Gumbel" ~ length(coef(modgumbel)),
          Distribuicoes == "GEV" ~ length(coef(modGEV)),
          TRUE ~ NA_integer_ # Caso haja alguma distribuição não contemplada
        )
      )

    p_value <- stats::pchisq(dfTAV$Likelihood_Ratio[1], df = dfTAV$df[1], lower.tail = FALSE)

    n <- length(ArquivPrec)
    alpha <- 0.05 # nível de significância

    valcrit_Qui <- get_critical_value(n, alpha)

    # Comparar p-valor
    if (p_value < valcrit_Qui) {
      # Escolher a distribuição com menor KS_Resultados se p-valor for significativo
      best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
    } else {
      # Escolher a distribuição com maior log-verossimilhança se p-valor não for significativo
      best_distribution <- dfTAV[which.max(dfTAV$Likelihood_Ratio), ]
    }
  } else {
    # Selecionar a distribuição com menor KS_Resultados
    best_distribution <- dfTA[which.min(dfTA$KS_Resultados), ]
  }

  my_list01 <- list(
    dfTAD = best_distribution,
    modnormal = modnormal,
    modln = modln,
    modgamma = modgamma,
    modexpo = modexpo,
    modweibull = modweibull,
    modgumbel = modgumbel,
    modGEV = modGEV
  )

  return(my_list01)
}
