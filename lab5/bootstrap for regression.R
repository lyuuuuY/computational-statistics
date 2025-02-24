data <- read.table("kresseertrag.dat", header = TRUE, sep = "", stringsAsFactors = FALSE)
x = data$X0.000
y = data$X215.5
regression_model <- lm(y~poly(x,3,raw = TRUE),data = data)
summary(regression_model)
model_reduced <- lm(y ~ I(x^2) + I(x^3) , data = data)
summary(model_reduced)
final_model <- model_reduced
pred <- predict(final_model) 
confint(final_model, level = 0.95)

plot(x,y, xlab = "Concentration", 
     ylab = "Yield",
     main = "Yield vs. Concentration with Fitted Curve")
lines(x, pred, col = "blue", lwd = 2)

#  95%-Bootstrap Confidence Interval (manual)

set.seed(12345) 
B <- 10000
n <- nrow(data)

beta2_boot <- numeric(B)

for (i in 1:B) {
  sample_id <- sample(1:n, size = n, replace = TRUE)
  data_boot <- data[sample_id, ]
  x = data_boot$X0.000
  y = data_boot$X215.5
  model_boot <- lm(y ~ I(x^2) + I(x^3), data = data_boot)
  
  beta2_boot[i] <- coef(model_boot)[2]
}

hist(beta2_boot, breaks = 30, main = "Bootstrap distribution of beta_2",
     xlab = "beta_2 estimates")
beta2_boot_sort <- sort(beta2_boot)
ci95 <- c(beta2_boot_sort[round(B*0.025)],beta2_boot_sort[round(B*0.975)])

#  95%-Bootstrap Confidence Interval (boot)
library(boot)
boot_func <- function(data, i){
  data_boot <- data[i, ]
  model <- lm(data_boot$X215.5 ~ I(data_boot$X0.000^2) + I(data_boot$X0.000^3))
  return(coef(model)[2])
}
get_boot <- boot(data = data, statistic = boot_func, R = 10000)
boot.ci(get_boot, type = "perc")
boot.ci(get_boot, type = "bca")
