#' This script generates the graph of the intensities
#'
#' @param Directory Parameters IDF Equation
#'
#' @import utils
#' @import graphics
#' @import grDevices
#' @import tiff
#'
#' @export
ScriptGraphic <- function(Directory) {
  path <- file.choose()

  TabCoef <- read.delim(path, header = TRUE, sep = ";")

  #Estimating intensities for different durations and periods of returns
  D <- rep(1:1440, 1)
  Tr <- c("2", "5", "10", "25", "50", "100")

  Int <- matrix(0, length(D), length(Tr))
  for (y in 1:length(TabCoef[, 1])) {
    cat("== Working on the data ==", "\n")
    cat("Wait...", "\n")
    cat("Station", "", TabCoef[y, 1], "\n")

    for (i in 1:length(Tr)) {
      for (j in 1:length(D)) {
        Int[j, i] <-
          (TabCoef[y, 4] * as.numeric(Tr[i]) ^ TabCoef[y, 5]) / ((as.numeric(D[j]) +
                                                                    TabCoef[y, 6]) ^ TabCoef[y, 7])
      }
    }

    tiff(
      paste0(Directory, "/", TabCoef[y, 1], ".tiff"),
      width = 8,
      height = 8,
      units = "in",
      res = 1200
    )

    windowsFonts(A = windowsFont("Rockwell"))

    #Plotting the graphs
    plot(
      Int[, 1],
      family = "A",
      cex.axis = 1.5,
      cex.lab = 1.5,
      cex.main = 1.5,
      col = "black",
      xlim = c(0, 420),
      main = TabCoef$Station[y],
      xlab = "Time (min)",
      ylab = "Intensity (mm/h)",
      type = "l",
      lwd = 1
    )
    lines(Int[, 2], col = "yellow")
    lines(Int[, 3], col = "red")
    lines(Int[, 4], col = "blue")
    lines(Int[, 5], col = "orange")
    lines(Int[, 6], col = "gray")
    legend(
      x = "topright",
      legend = c(
        "Tr= 2 years",
        "Tr= 5 years",
        "Tr= 10 years",
        "Tr= 25 years",
        "Tr= 50 years",
        "Tr= 100 years"
      ),
      col = c("black", "yellow", "red", "blue", "orange", "gray"),
      lty = 1,
      lwd = 1,
      xpd = TRUE,
      bty = "o",
      cex = 1.5,
      horiz = FALSE
    )
    grid(
      nx = NULL,
      ny = NULL,
      lty = 4,
      lwd = 0.5,
      col = "gray"
    )

    dev.off()

    cat("\014")

  }

  cat("\014")

  cat("Process finished!", "\n")

}
