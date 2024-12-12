## Setup
source("./scripts/setup.R")

# Assuming `pipe` is already defined and works as intended.

# User-defined arguments
posterior <- TRUE  # TRUE uses 'coef', FALSE uses 'estimate'
n <- 250
p <- 50
p1 <- 10
vars_to_plot <- paste0("V", formatC(1:20, width = 2, flag = "0"))  # Default: V01 to V20

penalties <- c("lasso", "SCAD", "MCP")
families <- c("gaussian", "binomial", "poisson")

# Generate one dataset per family
datasets <- lapply(families, function(fam) {
  hdrm::gen_data(family = fam, n = n, p = p, p1 = p1)
})
names(datasets) <- families

# Helper function to process a single penalty on a given dataset
process_data <- function(data, family, penalty, posterior) {
  # Run pipe function with given arguments
  results <- pipe(X = data$X, y = data$y, family = family, penalty = penalty, posterior = posterior)

  # Choose which estimate column to use
  estimate_col <- if (posterior) "coef" else "estimate"

  # Mutate the results to add truth and coverage
  results <- results %>%
    mutate(
      truth = data$beta,
      covered = lower <= truth & truth <= upper,
      cov_est = lower <= .data[[estimate_col]] & .data[[estimate_col]] <= upper,
      family = family,
      penalty = penalty,
      estimate = .data[[estimate_col]]
    )

  return(results)
}

# Apply all penalties to each of the three datasets (one per family)
all_data <- do.call(rbind, lapply(names(datasets), function(fam) {
  data_fam <- datasets[[fam]]
  do.call(rbind, lapply(penalties, function(pen) {
    process_data(data_fam, family = fam, penalty = pen, posterior = posterior)
  }))
}))

# Filter by the desired variables to plot
all_data <- all_data %>% filter(variable %in% vars_to_plot)

# Create the plot
ggplot(all_data, aes(y = variable)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  geom_point(aes(x = estimate), color = "black") +
  geom_point(aes(x = truth), color = "red", shape = 21) +
  facet_grid(rows = vars(family), cols = vars(penalty)) +
  theme_bw() +
  labs(
    x = "Coefficient Value",
    y = "Variable",
    title = paste(
      "Comparison of Intervals for", n, "obs,", p, "vars,", p1, "nonzero,",
      if (posterior) "Posterior=TRUE" else "Posterior=FALSE"
    )
  )
