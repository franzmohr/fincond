rm(list = ls())

library(dplyr)
library(readxl)
library(tidyr)


#### Macroeconomic variables ####

# Inflation
p <- read_xlsx("data-raw/koop.xlsx", sheet = "MACRO_P", na = "NaN") %>%
  rename(time = DATE)

# Unemployment
u <- read_xlsx("data-raw/koop.xlsx", sheet = "MACRO_RUC", na = "NaN") %>%
  rename(time = DATE)

# Output
y <- read_xlsx("data-raw/koop.xlsx", sheet = "MACRO_ROUTPUT", na = "NaN") %>%
  rename(time = DATE)

macro <- NULL
for (i in 2:ncol(p)) {
  temp <- ts(cbind(p[, i], u[, i], y[, i]),
             start = c(1970, 1), frequency = 4)
  temp <- na.omit(temp)
  dimnames(temp)[[2]] <- c("p", "u", "y")
  macro[[i - 1]] <- temp
}

#### Financial data ####

fin <- read_xlsx("data-raw/koop.xlsx", sheet = "Financial Data", skip = 3,
                 col_names = c("time", "sp500", "twexmmth", "oil", "spread_ted", "spread_10_2",
                               "spread_2_3m", "spread_comm", "loanhpi", "spread_30_mort",
                               "cmdebt", "cci", "move", "vix", "totalsl", "stdscom",
                               "mich1", "mich2", "mich3"),
                 na = "NaN")

ts_start_yr <- as.numeric(substring(fin[1, "time"], 1, 4))
ts_start_qtr <- as.numeric(substring(fin[1, "time"], 6, 6))

fin <- ts(fin[, -1], frequency = 4, start = c(ts_start_yr, ts_start_qtr))

#### Combine datasets ####
koop <- list("financial" = fin,
             "macro" = macro)

usethis::use_data(koop, overwrite = TRUE)
