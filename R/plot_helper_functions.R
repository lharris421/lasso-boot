# Function to calculate model results
#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
calculate_model_results <- function(data) {
  data %>%
    mutate(
      covered = lower <= truth & upper >= truth,
      mag_truth = abs(truth),
      covered = as.numeric(covered)
    )
}

# Function to perform fitting and prediction
#' Title
#'
#' @param data
#' @param x_values
#' @param method
#'
#' @return
#' @export
#'
#' @examples
predict_covered <- function(data, x_values, method) {
  # fit <- gam(covered ~ s(mag_truth) + s(group, bs = "re"), data = data, family = binomial)
  fit <- gam(covered ~ s(truth), data = data, family = binomial)
  y_values <- predict(fit, data.frame(truth = x_values, group = 101), type = "response")
  data.frame(x = x_values, y = y_values, method = method)
}

#' Title
#'
#' @param lower
#' @param upper
#' @param truth
#'
#' @return
#' @export
#'
#' @examples
miss_low <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    sign(truth) == 0 ~ 0,
    sign(truth) == 1 & sign(upper) == 1 & upper < truth  ~ 1,
    sign(truth) == -1 & sign(lower) == -1 & lower > truth ~ 1,
    sign(truth) != 0 & sign(lower) == 0 & sign(upper) == 0 ~ 1,
    TRUE ~ 0
  )
}
#' Title
#'
#' @param lower
#' @param upper
#' @param truth
#'
#' @return
#' @export
#'
#' @examples
miss_high <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    sign(truth) == 0 ~ 1,
    sign(truth) == 1 & lower > truth ~ 1,
    sign(truth) == -1 & upper < truth ~ 1,
    TRUE ~ 0
  )
}
#' Title
#'
#' @param lower
#' @param upper
#' @param truth
#'
#' @return
#' @export
#'
#' @examples
sign_inversion <- function(lower, upper, truth) {
  case_when(
    lower <= truth & upper >= truth ~ 0,
    sign(truth) == 0 ~ 0,
    sign(lower) == 0 & sign(upper) == 0 ~ 0,
    (sign(upper) != sign(truth)) & (sign(lower) != sign(truth))  ~ 1,
    # (sign(truth) == 1 & sign(upper) == -1) | (sign(truth) == 1 & sign(upper) == 0 & sign(lower) == -1) ~ 1,
    # (sign(truth) == -1 & sign(lower) == 1) | (sign(truth) == -1 & sign(lower) == 0 & sign(upper) == 1) ~ 1,
    TRUE ~ 0
  )
}
