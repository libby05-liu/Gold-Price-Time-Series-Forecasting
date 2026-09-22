# Gold Price Time-Series Analysis and Forecasting
# Portfolio project: unit-root testing, structural change, and 2022 forecast evaluation

# Data setup
library(forecast)
library(urca)

# Read the data file. If GoldPrice.csv is in the RStudio working directory,
# it is read automatically. Otherwise, select the supplied CSV when prompted.
if (file.exists("data/GoldPrice.csv")) {
  data_file <- "data/GoldPrice.csv"
} else if (file.exists("GoldPrice.csv")) {
  data_file <- "GoldPrice.csv"
} else {
  stop("Add GoldPrice.csv to the data/ folder before running this script.")
}

dt <- read.csv(data_file)

# Confirm that the required columns are present
required_columns <- c("Year", "Month", "GoldPrice")
if (!all(required_columns %in% names(dt))) {
  stop("The CSV must contain the columns Year, Month and GoldPrice.")
}

EndSample <- c(2021, 12)
EndForecast <- c(2022, 12)

GoldPrice <- ts(dt$GoldPrice,
                frequency = 12,
                start = c(2000, 1),
                end = EndForecast)
Y <- log(GoldPrice)
Time <- time(Y)

ts <- which((dt$Year < EndSample[1]) |
              (dt$Year == EndSample[1] & dt$Month <= EndSample[2]))
n <- length(ts)
tf <- (n + 1):length(Y)
Ys <- as.numeric(Y[ts])
Times <- as.numeric(Time[ts])


# 1. ADF test: ADF test with a linear trend
# The selected AR(2) levels model implies one augmented difference lag.
ADF1 <- ur.df(Ys, type = "trend", lags = 1)
print(summary(ADF1))

tADF1 <- as.numeric(ADF1@teststat[1])
pADF1 <- punitroot(tADF1, trend = "ct", statistic = "t")
critical5_ADF1 <- qunitroot(0.05,
                            N = length(Ys) - 2,
                            trend = "ct",
                            statistic = "t")
print(c(ADF_statistic = tADF1,
        critical_value_5pct = critical5_ADF1,
        p_value = pADF1))


# 2. Replicate the ADF regression: reproduce the ADF regression using lm()
t1 <- 3:n
DY1dep <- Ys[t1] - Ys[t1 - 1]
Ylag1 <- Ys[t1 - 1]
DYlag1 <- Ys[t1 - 1] - Ys[t1 - 2]

ADF1_lm <- lm(DY1dep ~ Times[t1] + Ylag1 + DYlag1)
print(round(summary(ADF1_lm)$coefficients, 6))
print(c(ur_df = tADF1,
        lm = summary(ADF1_lm)$coefficients["Ylag1", "t value"]))


# 3. Baseline differenced model: estimate the differenced model with AR(1) errors
D <- diff(Ys)
Model1 <- Arima(D,
                order = c(1, 0, 0),
                include.mean = TRUE,
                method = "ML")
print(Model1)

Model1_table <- cbind(
  Estimate = coef(Model1),
  Standard_error = sqrt(diag(Model1$var.coef)),
  t_statistic = coef(Model1) / sqrt(diag(Model1$var.coef)))
print(round(Model1_table, 6))


# 4. ADF test with a June 2011 trend break: ADF test with June 2011 trend break
b <- 2011 + 5/12
TB <- pmax(Times - b, 0)
t2 <- 4:n

DY2dep <- Ys[t2] - Ys[t2 - 1]
Ylag2 <- Ys[t2 - 1]
DYlag1_2 <- Ys[t2 - 1] - Ys[t2 - 2]
DYlag2_2 <- Ys[t2 - 2] - Ys[t2 - 3]

ADF2_lm <- lm(DY2dep ~ Times[t2] + TB[t2] + Ylag2 +
                DYlag1_2 + DYlag2_2)
ADF2_table <- summary(ADF2_lm)$coefficients
print(round(ADF2_table, 6))

tADF2 <- ADF2_table["Ylag2", "t value"]
print(c(ADF_statistic = tADF2))


# 5. Simulated critical value and p-value: simulated critical value and p-value
# Each replication is a random walk under the unit-root null and uses the
# same sample size, deterministic regressors and two augmented difference lags.
set.seed(42)
reps <- 10000
tADF_sim <- numeric(reps)

for (r in seq_len(reps)) {
  simY <- c(0, cumsum(rnorm(n - 1)))
  simDY <- simY[t2] - simY[t2 - 1]
  simYlag <- simY[t2 - 1]
  simDYlag1 <- simY[t2 - 1] - simY[t2 - 2]
  simDYlag2 <- simY[t2 - 2] - simY[t2 - 3]

  simFit <- lm(simDY ~ Times[t2] + TB[t2] + simYlag +
                 simDYlag1 + simDYlag2)
  tADF_sim[r] <- summary(simFit)$coefficients["simYlag", "t value"]
}

simulated_critical_values <- quantile(tADF_sim,
                                      probs = c(0.01, 0.05, 0.10))
simulated_p_value <- mean(tADF_sim <= tADF2)
print(simulated_critical_values)
print(c(ADF_statistic = tADF2,
        simulated_p_value = simulated_p_value))


# 6. Differenced model with a post-break mean: estimate the differenced model with a post-break mean
TimeD <- Times[-1]
Post <- as.numeric(TimeD > b)
Post_matrix <- matrix(Post, ncol = 1,
                      dimnames = list(NULL, "Post"))

Model2 <- Arima(D,
                order = c(2, 0, 0),
                xreg = Post_matrix,
                include.mean = TRUE,
                method = "ML")
print(Model2)

Model2_table <- cbind(
  Estimate = coef(Model2),
  Standard_error = sqrt(diag(Model2$var.coef)),
  t_statistic = coef(Model2) / sqrt(diag(Model2$var.coef)))
print(round(Model2_table, 6))

# Pre-break and post-break monthly mean log changes
pre_break_mean <- coef(Model2)["intercept"]
post_break_mean <- coef(Model2)["intercept"] + coef(Model2)["Post"]
print(c(Pre_break_mean = pre_break_mean,
        Post_break_mean = post_break_mean))


# 7. Forecast evaluation: fixed-origin 1- to 12-step-ahead forecasts for 2022
DYf1 <- as.numeric(forecast(Model1, h = 12)$mean)

Post_future <- matrix(1, nrow = 12, ncol = 1,
                      dimnames = list(NULL, "Post"))
DYf2 <- as.numeric(forecast(Model2, h = 12,
                            xreg = Post_future)$mean)

Yf1 <- Ys[n] + cumsum(DYf1)
Yf2 <- Ys[n] + cumsum(DYf2)
Yactual <- as.numeric(Y[tf])

Error1 <- Yactual - Yf1
Error2 <- Yactual - Yf2
RMSE1 <- sqrt(mean(Error1^2))
RMSE2 <- sqrt(mean(Error2^2))

ForecastTable <- data.frame(
  Month = month.abb,
  Actual = Yactual,
  No_break_forecast = Yf1,
  No_break_error = Error1,
  Break_forecast = Yf2,
  Break_error = Error2)

ForecastTable_print <- ForecastTable
ForecastTable_print[-1] <- round(ForecastTable_print[-1], 6)
print(ForecastTable_print)
print(c(No_break_RMSE = RMSE1,
        Trend_break_RMSE = RMSE2))

RMSE_reduction <- 100 * (RMSE1 - RMSE2) / RMSE1
print(c(RMSE_reduction_percent = RMSE_reduction))


# 7. Forecast evaluation: export and display the forecast comparison plot
draw_forecast_plot <- function() {
  plot(1:12, Yactual,
       type = "n",
       ylim = range(c(Yactual, Yf1, Yf2)),
       xaxt = "n",
       xlab = "Month in 2022",
       ylab = "Log gold price",
       main = "Fixed-origin multi-step forecasts for 2022")

  axis(1, at = 1:12, labels = month.abb)
  grid(nx = NA, ny = NULL, col = "grey85", lty = 1)

  lines(1:12, Yactual,
        type = "o", pch = 16, lwd = 2, col = "black")
  lines(1:12, Yf1,
        type = "o", pch = 15, lwd = 2, col = "#1F4E79")
  lines(1:12, Yf2,
        type = "o", pch = 17, lwd = 2, col = "#B36B00")

  legend("bottomleft",
         legend = c("Actual", "No trend break",
                    "June 2011 trend break"),
         col = c("black", "#1F4E79", "#B36B00"),
         pch = c(16, 15, 17),
         lty = 1,
         lwd = 2,
         bty = "n")
}

png(filename = "figures/forecast_comparison.png",
    width = 2400,
    height = 1600,
    res = 300)
draw_forecast_plot()
dev.off()

# Display the same plot in the RStudio Plots pane
draw_forecast_plot()
