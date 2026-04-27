##' @keywords internal
##' @importFrom utils flush.console
##' @importFrom graphics points
##' @importFrom sf st_nearest_feature
##' @importFrom Deriv Deriv
##' @importFrom numDeriv grad hessian
##' @importFrom ggplot2 ggplot aes geom_line geom_point theme labs
##' @importFrom ggplot2 facet_wrap as_labeller
##' @importFrom ggplot2 geom_histogram geom_vline
##' @importFrom ggplot2 scale_x_continuous scale_color_manual scale_fill_manual scale_linetype_manual
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 theme_bw
##' @importFrom rlang sym
##' @importFrom gridExtra grid.arrange
##' @importFrom stats plogis
##' @importFrom dplyr mutate first
##' @importFrom tibble as_tibble
##' @importFrom sf st_as_text
NULL

utils::globalVariables(c("lo", "hi", "series", "obs", "n_obs"))
