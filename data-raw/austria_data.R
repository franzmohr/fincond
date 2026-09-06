rm(list = ls())

library(dplyr)
library(OECD)
library(seasonal)
library(tidyr)
library(zoo)


# Real GDP growth
temp <- get_dataset("QNA", filter = "AUT.B1_GE.VOBARSA.Q") %>%
  arrange(obsTime) %>%
  mutate(value = log(obsValue / 4)) %>%
  select(obsTime, value)

ts_start_yr <- as.numeric(substring(temp[1, "obsTime"], 1, 4))
ts_start_qtr <- as.numeric(substring(temp[1, "obsTime"], 7, 7))

y <- ts(temp[, "value"], frequency = 4, start = c(ts_start_yr, ts_start_qtr))
rm(list = c("ts_start_yr", "ts_start_qtr", "temp"))

# Inflation
download.file("https://www.bis.org/statistics/full_webstats_long_cpi_dataflow_csv.zip", destfile = "data-raw/cpi.zip")
unzip("data-raw/cpi.zip", exdir = "data-raw")

temp <- read.csv("data-raw/WEBSTATS_LONG_CPI_DATAFLOW_csv_col.csv", stringsAsFactors = FALSE) %>%
  filter(FREQ == "M",
         UNIT_MEASURE == 628,
         REF_AREA == "AT") %>%
  gather(key = "date", value = "value", -(FREQ:Time.Period)) %>%
  filter(!is.na(value)) %>%
  mutate(date = as.yearmon(date, "X%Y.%m"),
         obsTime = as.yearqtr(date)) %>%
  group_by(obsTime) %>%
  filter(n() == 3) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  mutate(obsTime = as.character(obsTime)) %>%
  filter(!is.na(value))

ts_start_yr <- as.numeric(substring(temp[1, "obsTime"], 1, 4))
ts_start_qtr <- as.numeric(substring(temp[1, "obsTime"], 7, 7))

Dp <- ts(temp[, "value"], frequency = 4, start = c(ts_start_yr, ts_start_qtr))
Dp <- log(Dp)
rm(list = c("ts_start_yr", "ts_start_qtr", "temp"))

Dp <- seas(diff(Dp))$data[, "final"]

unlink("data-raw/cpi.zip")
unlink("data-raw/WEBSTATS_LONG_CPI_DATAFLOW_csv_col.csv")

# Unemployment rate
temp <- get_dataset("MEI", filter = "AUT.LMUNRRTT.STSA.Q") %>%
  arrange(obsTime) %>%
  mutate(value = obsValue) %>%
  select(obsTime, value)

ts_start_yr <- as.numeric(substring(temp[1, "obsTime"], 1, 4))
ts_start_qtr <- as.numeric(substring(temp[1, "obsTime"], 7, 7))

u <- ts(temp[, "value"], frequency = 4, start = c(ts_start_yr, ts_start_qtr))
rm(list = c("ts_start_yr", "ts_start_qtr", "temp"))


# Make ts object
AT <- cbind(diff(y, 1), u, Dp) %>%
  na.omit()

dimnames(AT)[[2]] <- c("y", "u", "Dp")

plot(AT)

# Save in data
usethis::use_data(AT, overwrite = TRUE)

