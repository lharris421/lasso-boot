set.seed(189807771)
X1 <- rnorm(50)
X2 <- sqrt(0.75)*rnorm(50) + 0.5*X1
Y <- X1 - X2 + rnorm(50)
hist(Y)

## 0th Order
confint(lm(Y ~ X1))[2,]

## 1st Order
resid <- Y - fitted(lm(Y ~ X2))
confint(lm(resid ~ 0 + X1))[1,]
## OR (more like what we are doing for a single iteration)

joint_mod <- lm(Y ~ X1 + X2)
partial_resid <- Y - fitted(joint_mod) + joint_mod$coefficients[2]*X1
confint(lm(partial_resid ~ 0 + X1))[1,]

X <- cbind(1, X1, X2)
Xzj <- cbind(1, 0, X2)
proj <- (diag(length(Y)) - Xzj%*%solve(t(X) %*% X)%*%t(X))
R <- proj %*% Y
confint(lm(R ~ 0 + X1))[1,]

## 2nd Order
confint(lm(Y ~ X1 + X2))[2,]

###
X <- X2
proj <- (diag(length(Y)) - X%*%solve(t(X) %*% X)%*%t(X))
R <- as.numeric(proj %*% Y)
confint(lm(R ~ as.numeric(proj %*% X1)))[2,]

# Projection (2nd order...although maybe 1st order versions exist in model selection context?)
X <- model.matrix(~X2)
Q <- diag(length(Y)) - X %*% solve(crossprod(X)) %*% t(X)
yy <- Q %*% Y
XX <- Q %*% X1
lm(yy ~ 0 + XX) |> confint()
