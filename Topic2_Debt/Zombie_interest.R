library(tidyverse)
library(tibbletime)
library(lubridate)
library(readxl)

CMMPI <- read_excel("data/CMMPI_Irfil.xlsx")

CMMPI <- CMMPI %>% 
  mutate(
    begin_date = Datesgn %>% ymd(),
    short_rate = if_else(is.na(Irfil02), Irfil28, Irfil02),
    long_rate = if_else(is.na(Irfil04), Irfil29, Irfil04)
  ) %>% 
  select(begin_date, short_rate, long_rate)


CMMPI_year <- CMMPI %>% 
  mutate(
    end_date = lead(begin_date, default = ymd('2017-12-31')),
    date = map2(begin_date, end_date-1, ~seq(.x, .y, by = "day"))
  ) %>% 
  unnest(date) %>% 
  mutate(
    year = year(date)
  ) %>% 
  group_by(year) %>% 
  summarise_at(vars(contains("rate")), mean) %>% 
  write_csv("Zombie_data/Loan_interest.csv")
  
#####################################

BND <- read_excel("data/BND_Ccbdinfo.xlsx")

BND <- BND %>% 
  select(Bndyer, Intrrate)

get_intrrate <- function(year){
  BND %>% 
    filter(Bndyer <= year, Bndyer > year - 5) %>% 
    summarise(Intrrate = min(Intrrate)) %>% 
    pull(Intrrate)
}

tibble(
  year  = 2000:2017
) %>% 
  mutate(
    Bond_rate = year %>% map_dbl(get_intrrate)
  ) %>% 
  write_csv("Zombie_data/Bond_interest.csv")

