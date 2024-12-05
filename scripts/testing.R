rm(list = ls())
source("./scripts/setup.R")
library(ncvreg)

# hdrm::attachData("lister")
# res1 <- pipe(X, y, penalty = "MCP", lambda = 0.059)
# res2 <- pipe_ncvreg(X, y, lambda = 0.059, penalty = "MCP")
#
# res1 %>%
#   filter(p.adjust <= 0.1)
#
# res2 %>%
#   filter(significance <= 0.1)

# data <- hdrm::gen_data(n = 1000, p = 50, p1 = 10)
# pipe(data$X, data$y, penalty = "lasso") %>%
#   mutate(
#     truth = data$beta,
#     covered = lower <= truth & truth <= upper
#   ) %>%
#   pull(covered) %>%
#   mean()
#
data <- hdrm::gen_data(family = "binomial", n = 400, p = 50, p1 = 10)
pipe(data$X, data$y, family = "binomial", penalty = "lasso") %>%
  mutate(
    truth = data$beta,
    covered = lower <= truth & truth <= upper
  ) %>%
  pull(covered) %>%
  mean()


data <- hdrm::gen_data(family = "poisson", n = 1000, p = 50, p1 = 10)
pipe(data$X, data$y, family = "poisson", penalty = "lasso") %>%
  mutate(
    truth = data$beta,
    covered = lower <= truth & truth <= upper
  ) %>%
  pull(covered) %>%
  mean()

tmp <- cv.ncvreg(data$X, data$y, family = "poisson", penalty ="lasso")
hist(tmp$fit$y - predict(tmp$fit, data$X, lambda = tmp$lambda.min, type = "response"))


coef(tmp)
