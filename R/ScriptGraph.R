#' This script calculates the IDF-Graph
#'
#' @param Directory Input Location where the IDF equations
#'
#' @import grDevices
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
#' @import utils
#'
#' @export
ScriptGraph <- function(Directory) {
  # Choose the file
  path <- file.choose()

  # Read data
  TabCoef <- utils::read.delim(path, header = TRUE, sep = ";") %>%
    stats::na.omit()

  # Time and duration
  Time <- c(2, 50, 100) # Time (years)
  Duration <- seq(6, 1440, by = 1) # Duration (minutes)

  # Total number of figures
  num_figures <- nrow(TabCoef)

  # Set up parameters for arranging figures
  figures_per_page <- 4
  num_pages <- ceiling(num_figures / figures_per_page)

  # Loop over pages
  for (page in 1:num_pages) {
    # Set up plot parameters
    tiff(
      filename = paste0(Directory,"IDF-Graph_Page_", page, ".tiff"),
      width = 8,
      height = 8,
      units = "in",
      res = 1200
    )
    par(mfrow = c(2, 2), pty = "s")

    # Determine the figures to plot on this page
    start_figure <- (page - 1) * figures_per_page + 1
    end_figure <- min(start_figure + figures_per_page - 1, num_figures)

    # Plotting each graph for this page
    for (k in start_figure:end_figure) {
      Intensity <- matrix(0, nrow = length(Duration), ncol = length(Time))
      for (i in 1:length(Time)) {
        for (j in 1:length(Duration)) {
          Intensity[j, i] <- (TabCoef$K[k] * Time[i]^TabCoef$m[k]) / ((Duration[j] + TabCoef$t0[k])^TabCoef$n[k])
        }
      }
      # Plotting each graph
      plot(
        Intensity[, 1],
        col = "black",
        xlim = c(0, 400),
        family = "Arial",
        cex.axis = 1,
        cex.lab = 1,
        cex.main = 1,
        main = paste(TabCoef$Station[k]),
        xlab = "Time (min)",
        ylab = expression("Intensity (mm." ~ h^-1 * ")"),
        type = "l",
        lwd = 2
      )
      lines(
        Intensity[, 2],
        col = "red",
        lty = 2,
        lwd = 2
      )
      lines(
        Intensity[, 3],
        col = "orange",
        lty = 3,
        lwd = 2
      )
      legend(
        "topright",
        legend = c("Tr-2", "Tr-50", "Tr-100"),
        col = c("black", "red", "orange"),
        lty = 1:3,
        cex = 0.8
      )
    }

    dev.off() # Close the device
  }

  cat("\014")
}
