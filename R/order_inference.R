set.seed(189807771)
X1 <- rnorm(50)
X2 <- sqrt(0.75)*rnorm(50) + 0.5*X1
Y <- X1 - X2 + rnorm(50)
hist(Y)

## 0th Order
confint(lm(Y ~ X1))[2,]

## 1st Order
resid <- Y - fitted(lm(Y ~ X2))
confint(lm(resid ~ X1))[2,]
## OR (more like what we are doing for a single iteration)
joint_mod <- lm(Y ~ X1 + X2)
partial_resid <- Y - fitted(joint_mod) + joint_mod$coefficients[2]*X1
confint(lm(partial_resid ~ X1))[2,]

## 2nd Order
confint(lm(Y ~ X1 + X2))[2,]
