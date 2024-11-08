#' Title
#'
#' @param res
#' @param params
#'
#' @return
#' @export
#'
#' @examples
bias_decomp_plots <- function(res, params) {

  # Generate the original data to get true beta values
  true_betas <- res$true_betas
  true_a1 <- true_betas["A1"]

  # Initialize data frame with original biases
  original_data <- data.frame(
    original_bias = true_a1 - res$orig_est,
    original_n_bias = res$orig_n_bias,
    original_err_bias = res$orig_err_bias,
    bias = true_a1 - apply(res$ests, 1, mean),
    n_bias = apply(res$n_bias, 1, mean),
    err_bias = apply(res$err_bias, 1, mean)
  )

  # Handle other A variables (excluding A1)
  if (!is.null(res$orig_a_bias)) {
    other_a_vars <- colnames(res$orig_a_bias)
    for (var in other_a_vars) {
      original_data[[paste0("original_a_bias_", var)]] <- res$orig_a_bias[, var]
      original_data[[paste0("a_bias_", var)]] <- apply(res$a_bias[, , var], 1, mean)
    }
  }

  # Handle B variables
  if (!is.null(res$orig_b_bias)) {
    b_vars <- colnames(res$orig_b_bias)
    for (var in b_vars) {
      original_data[[paste0("original_b_bias_", var)]] <- res$orig_b_bias[, var]
      original_data[[paste0("b_bias_", var)]] <- apply(res$b_bias[, , var], 1, mean)
    }
  }

  # Compute sum of original biases
  original_bias_vars <- grep("^original_.*_bias", names(original_data), value = TRUE)
  boot_bias_vars <- names(original_data)[!(names(original_data) %in% c(original_bias_vars, "original_bias", "bias"))]
  original_data$sum_original_bias <- rowSums(original_data[, original_bias_vars], na.rm = TRUE)

  # Compute sum of bootstrap biases
  original_data$sum_bias <- rowSums(original_data[, boot_bias_vars], na.rm = TRUE)

  # Compute additional biases by subtracting bootstrap biases from original biases
  plotting_data <- original_data %>%
    mutate(
      add_bias = bias - original_bias,
      add_n_bias = n_bias - original_n_bias,
      add_err_bias = err_bias - original_err_bias,
      .keep = "all"
    )

  # Additional biases for other A variables
  if (!is.null(res$orig_a_bias)) {
    for (var in other_a_vars) {
      plotting_data[[paste0("add_a_bias_", var)]] <- plotting_data[[paste0("a_bias_", var)]] - plotting_data[[paste0("original_a_bias_", var)]]
    }
  }

  # Additional biases for B variables
  if (!is.null(res$orig_b_bias)) {
    for (var in b_vars) {
      plotting_data[[paste0("add_b_bias_", var)]] <- plotting_data[[paste0("b_bias_", var)]] - plotting_data[[paste0("original_b_bias_", var)]]
    }
  }

  # Compute sum of additional biases
  add_bias_vars <- c("add_n_bias", "add_err_bias")

  if (!is.null(res$orig_a_bias)) {
    add_a_bias_vars <- paste0("add_a_bias_", other_a_vars)
    add_bias_vars <- c(add_bias_vars, add_a_bias_vars)
  }

  if (!is.null(res$orig_b_bias)) {
    add_b_bias_vars <- paste0("add_b_bias_", b_vars)
    add_bias_vars <- c(add_bias_vars, add_b_bias_vars)
  }

  # Sum of all additional biases
  plotting_data$sum_add_bias <- rowSums(plotting_data[, add_bias_vars], na.rm = TRUE)

  # Prepare variables for plotting
  add_bias_vars <- c("add_bias", add_bias_vars, "sum_add_bias")

  # Reshape data into long format for ggplot2
  long_data <- pivot_longer(
    plotting_data %>% select(all_of(add_bias_vars)),
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

  # Create labels for variables
  variable_labels_add <- c(
    "add_bias" = "Additional Bias",
    "add_n_bias" = "Attributable to N",
    "add_err_bias" = "Attr. to Err Corr (Irreducible)"
  )

  # Labels for other A variables
  if (!is.null(res$orig_a_bias)) {
    for (var in other_a_vars) {
      var_name <- paste0("add_a_bias_", var)
      variable_labels_add[var_name] <- paste("Attributable to", var)
    }
  }

  # Labels for B variables
  if (!is.null(res$orig_b_bias)) {
    for (var in b_vars) {
      var_name <- paste0("add_b_bias_", var)
      variable_labels_add[var_name] <- paste("Attributable to", var)
    }
  }

  variable_labels_add["sum_add_bias"] <- "Sum Attributable"

  # Map variables to labels
  long_data$variable <- factor(
    long_data$variable,
    levels = names(variable_labels_add),
    labels = variable_labels_add[names(variable_labels_add) %in% long_data$variable]
  )

  # Similarly reshape original_data for the original biases plot
  original_bias_vars_plot <- c("original_bias", original_bias_vars, "sum_original_bias")
  long_data_orig <- pivot_longer(
    original_data %>% select(all_of(original_bias_vars_plot)),
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

  # Create labels for original biases
  variable_labels_orig <- c(
    "original_bias" = "Bias",
    "original_n_bias" = "Attributable to N",
    "original_err_bias" = "Attr. to Err Corr (Irreducible)"
  )

  # Labels for other A variables
  if (!is.null(res$orig_a_bias)) {
    for (var in other_a_vars) {
      var_name <- paste0("original_a_bias_", var)
      variable_labels_orig[var_name] <- paste("Attributable to", var)
    }
  }

  # Labels for B variables
  if (!is.null(res$orig_b_bias)) {
    for (var in b_vars) {
      var_name <- paste0("original_b_bias_", var)
      variable_labels_orig[var_name] <- paste("Attributable to", var)
    }
  }

  variable_labels_orig["sum_original_bias"] <- "Sum Attributable"

  # Map variables to labels
  long_data_orig$variable <- factor(
    long_data_orig$variable,
    levels = names(variable_labels_orig),
    labels = variable_labels_orig[names(variable_labels_orig) %in% long_data_orig$variable]
  )

  # Reshape data for bootstrap biases plot
  bias_vars_plot <- c("bias", boot_bias_vars, "sum_bias")
  long_data_boot <- pivot_longer(
    original_data %>% select(all_of(bias_vars_plot)),
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

  # Create labels for bootstrap biases
  variable_labels_boot <- c(
    "bias" = "Bias",
    "n_bias" = "Attributable to N",
    "err_bias" = "Attr. to Err Corr (Irreducible)"
  )

  # Labels for other A variables
  if (!is.null(res$a_bias)) {
    for (var in other_a_vars) {
      var_name <- paste0("a_bias_", var)
      variable_labels_boot[var_name] <- paste("Attributable to", var)
    }
  }

  # Labels for B variables
  if (!is.null(res$b_bias)) {
    for (var in b_vars) {
      var_name <- paste0("b_bias_", var)
      variable_labels_boot[var_name] <- paste("Attributable to", var)
    }
  }

  variable_labels_boot["sum_bias"] <- "Sum Attributable"

  # Map variables to labels
  long_data_boot$variable <- factor(
    long_data_boot$variable,
    levels = names(variable_labels_boot),
    labels = variable_labels_boot[names(variable_labels_boot) %in% long_data_boot$variable]
  )

  # Create density plots using ggplot2
  add_bias <- ggplot(long_data, aes(x = value, fill = variable, color = variable)) +
    geom_density(alpha = 0.2) +
    theme_bw() +
    labs(
      title = "Additional Bootstrap Bias Compared to Original Debiased Est (Decomp)",
      x = "Bias", y = "Density", color = "Bias Type", fill = "Bias Type"
    )

  orig_bias <- ggplot(long_data_orig, aes(x = value, fill = variable, color = variable)) +
    geom_density(alpha = 0.2) +
    theme_bw() +
    labs(
      title = "Original Debiased Est Bias Decomp",
      x = "Bias", y = "Density", color = "Bias Type", fill = "Bias Type"
    ) +
    geom_vline(xintercept = median(res$lambdas), color = "red")

  boot_bias <- ggplot(long_data_boot, aes(x = value, fill = variable, color = variable)) +
    geom_density(alpha = 0.2) +
    theme_bw() +
    labs(
      title = "Bootstrap Bias Decomp",
      x = "Bias", y = "Density", color = "Bias Type", fill = "Bias Type"
    ) +
    geom_vline(xintercept = median(res$lambdas), color = "red")

  return(list(add_bias, orig_bias, boot_bias))

}
