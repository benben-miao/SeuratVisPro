#' @title Publication-ready ggplot theme
#' @description Consistent theme used across SeuratVisPro plots with optional grid removal.
#' @author benben-miao
#'
#' @return A `ggplot2` theme object.
#' @param grid Logical; if `FALSE`, removes panel grid.
#'
#' @keywords internal
#'
svpp_theme <- function() {
  family <- "Arial"
  psfonts <- names(grDevices::postscriptFonts())
  if (!is.null(psfonts) && !(family %in% psfonts))
    family <- "sans"
  if (!interactive())
    family <- "sans"

  ggplot2::theme_light() +
    ggplot2::theme(
      text = ggplot2::element_text(family = family),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      axis.title = ggplot2::element_text(
        size = 10,
        face = "bold",
        color = "black"
      ),
      legend.title = ggplot2::element_text(size = 8, color = "black"),
      legend.text = ggplot2::element_text(size = 8, color = "black"),
      legend.key.size = ggplot2::unit(0.5, 'cm')
    )
}
