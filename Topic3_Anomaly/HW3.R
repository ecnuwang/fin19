# 加载需要的R包
library(tidyverse)
library(readxl)
library(lubridate)
library(lmtest)
library(sandwich)
library(broom)
library(stargazer)
library(kableExtra)
library(tibbletime)

# 基本信息
cominfo <- read_excel("Anomaly_data/TRD_Co.xlsx")

# 财务指标
comfin <- read_excel("Anomaly_data/FI_T2.xlsx", col_names = c("Stkcd", "Accper", "Indcd", "nr", "npexnr", "ROE", "ROE_exnr"), skip = 3)

# 换手率数据
turnover <- read_excel("data/Liq_Tover_M.xlsx")

# 个股交易数据
trdmnth <- read_excel("Anomaly_data/TRD_Mnth.xlsx")

# 无风险利率
nrrate <- read_excel("Anomaly_data/TRD_Nrrate.xlsx")

# 五因子
facmnth <- read_excel("Anomaly_data/STK_MKT_FivefacMonth.xlsx")

comfin <- comfin %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper),
    month = month(Accper),
    npexnr = npexnr %>% accumulate(~if_else(is.na(.y), .x, .y)),
    ROE = ROE %>% accumulate(~if_else(is.na(.y), .x, .y))
  ) %>% 
  filter(year <= 2018)

trdmnth <- trdmnth %>% 
  separate(Trdmnt, into = c("year", "month"), sep = "-") %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month)
  )

nrrate <- nrrate %>% 
  mutate(
    Clsdt = ymd(Clsdt),
    year = year(Clsdt),
    month = month(Clsdt)
  ) %>% 
  group_by(year, month) %>% 
  summarise(rf = mean(Nrrmtdt))

# 挑选综合A股和创业板的因子
facmnth <- facmnth %>% 
  filter(MarkettypeID == "P9709", Portfolios == 1) %>% 
  separate(TradingMonth, into = c("year", "month"), sep = "-") %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month)
  ) %>% 
  mutate_at(vars(RiskPremium1:CMA2), ~.*100)

# 计算超额换手率
mean_roll <- rollify(~mean(.x, na.rm = T), window = 12)

turnover <- turnover %>% 
  separate(Trdmnt, into = c("year", "month"), sep = "-") %>% 
  mutate_at(vars(year:MarketType), ~as.numeric(.)) %>% 
  arrange(Stkcd, year, month) %>% 
  group_by(Stkcd) %>% 
  filter(n() > 12) %>% 
  mutate(
    avgtover = mean_roll(ToverTlMAvg),
    abtover = lag(ToverTlMAvg / avgtover)
  ) %>% 
  filter(year <= 2018)

# 筛选企业，合并数据
comA <- cominfo %>% 
  filter(
    Indcd != "0001" & Indcd != "0002",
    Markettype == 1 | Markettype == 4 | Markettype == 16
  )

comdata <- comfin %>% 
  inner_join(comA, by = c("Stkcd" = "Stkcd"))

trddata <- trdmnth %>% 
  inner_join(comA, by = c("Stkcd" = "Stkcd")) %>% 
  inner_join(turnover, by = c("Stkcd" = "Stkcd", "year" = "year", "month" = "month")) %>%
  mutate(
    year1 = if_else(month == 1 | month == 2 | month == 3, year - 1, year),
    month1 = case_when(
      month == 1 | month == 2 | month == 3 ~ 12,
      month == 4 | month == 5 | month == 6 ~ 3,
      month == 7 | month == 8 | month == 9 ~ 6,
      TRUE ~ 9
    )
  ) %>% 
  inner_join(comdata, by = c("Stkcd" = "Stkcd", "year1" = "year", "month1" = "month")) %>% 
  complete(Stkcd, year, month) %>% 
  arrange(Stkcd, year, month) %>% 
  group_by(Stkcd) %>% 
  mutate(
    size = lag(Mclsprc) * (Msmvttl/Mclsprc),
    EP = npexnr / size
  ) %>% 
  group_by(year, month) %>% 
  filter(size >= quantile(size, 0.3, na.rm = T)) %>% 
  left_join(nrrate) %>% 
  mutate(
    exret = Mretwd*100 - rf
  ) 

# 因子分组
trddata_grp <- trddata %>% 
  group_by(year, month) %>% 
  mutate(
    grp_size = ntile(size, 2),
    grp_value = case_when(
      is.na(EP) ~ NA_real_,
      EP >= quantile(EP, 0.7, na.rm = T) ~ 1,
      EP <= quantile(EP, 0.7, na.rm = T) & EP >= quantile(EP, 0.3, na.rm = T) ~ 2,
      EP <= quantile(EP, 0.3, na.rm = T) ~ 3,
      TRUE ~ NA_real_
    ),
    grp_abto = ntile(abtover, 10),
    grp_size1 = ntile(size, 10)
  ) %>% 
  group_by(year, month, grp_size1) %>% 
  mutate(grp_abto_con = ntile(abtover, 10)) %>% 
  ungroup()

# 构造因子
# Value-weighted return
factor1 <- trddata_grp %>% 
  group_by(year, month) %>% 
  summarise(
    MKT = sum(exret*size, na.rm = TRUE) / sum(size, na.rm = TRUE)
  ) %>% 
  ungroup()

factor2 <- trddata_grp %>% 
  group_by(year, month, grp_size, grp_value) %>% 
  summarise(
    vwret = sum(Mretwd*100*size, na.rm = TRUE) / sum(size, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(!is.na(grp_size), !is.na(grp_value)) %>% 
  mutate(grp = grp_size*10 + grp_value) %>% 
  select(year, month, grp, vwret) %>% 
  spread(grp, vwret) %>% 
  mutate(
    SMB = 1/3*(`11` + `12` + `13`) - 1/3*(`21` + `22` + `23`),
    VMG = 1/2*(`11` + `21`) - 1/2*(`13` + `23`)
  )

# 与Fama-French因子数据合并
factor <- factor2 %>%
  left_join(factor1) %>% 
  left_join(facmnth)

# 换手率因子
factor3 <- trddata_grp %>% 
  group_by(year, month, grp_abto) %>% 
  summarise(
    vwret = sum(Mretwd*100*size, na.rm = TRUE) / sum(size, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(!is.na(grp_abto)) %>% 
  spread(grp_abto, vwret) %>% 
  mutate(
    liq = `1` - `10`
  ) %>% 
  left_join(factor)

factor3 <- trddata_grp %>% 
  group_by(year, month, grp_abto_con) %>% 
  summarise(
    vwret = sum(Mretwd*100*size, na.rm = TRUE) / sum(size, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(!is.na(grp_abto_con)) %>% 
  spread(grp_abto_con, vwret) %>% 
  mutate(
    liq = `1` - `10`
  ) %>% 
  left_join(factor)

# Factor model
## CAPM
capmtest1 <- lm(liq ~ MKT, factor3)
capmsd1 <- coeftest(capmtest1, vcov = NeweyWest(capmtest1, lag = 4))

## 3-factor model
thftest1 <- lm(liq ~ MKT + SMB + VMG, factor3)
thfsd1 <- coeftest(thftest1, vcov = NeweyWest(thftest1, lag = 4))

## Fama-French 3-factor model
thftest1_FF <- lm(liq ~ MKT + SMB2 + HML2, factor3)
thfsd1_FF <- coeftest(thftest1_FF, vcov = NeweyWest(thftest1_FF, lag = 4))

stargazer(capmsd1, thfsd1, thfsd1_FF, 
          intercept.bottom = FALSE, intercept.top = TRUE, 
          style = "aer", omit.table.layout = "n")

# Fama-Macbeth
stk_beta <- read_csv("Anomaly_data/beta.csv")

trddata_fm <- trddata %>% 
  left_join(stk_beta, by = c("Stkcd" = "Stkcd", "year" = "year", "month" = "month")) %>% 
  rename(beta = estimate) %>% 
  mutate(
    lnMe = log(size),
    lnEP_pos = if_else(EP > 0, log(EP), 0),
    EP_dum = if_else(EP <= 0, 1, 0)
  )

fm_model <- function(model){
  fmcoef <- trddata_fm %>%
    nest(-year, -month) %>%
    mutate(
      coefs = data %>% map(~lm(model, .)) %>% map(tidy)
    ) %>%
    select(year, month, coefs) %>%
    unnest(coefs) %>% 
    nest(-term) %>%
    mutate(
      fmtest = data %>% map(~lm(estimate ~ 1, .)),
      fmsd = fmtest %>% map(~coeftest(., vcov = NeweyWest(., lag = 4))) %>% map(tidy)
    ) %>%
    select(term, fmsd) %>%
    unnest() %>% 
    mutate(varnm = str_c(round(estimate, 3), "\n", "(", round(statistic, 2), ")")) %>% 
    select(term, varnm) 
  
  fmrsq <- trddata_fm %>%
    nest(-year, -month) %>%
    mutate(
      rsq = data %>% map(~lm(model, .)) %>% map(glance),
      N = data %>% map(~count(.))
    ) %>%
    select(year, month, rsq, N) %>%
    unnest(rsq, N) %>% 
    summarise_at(vars(r.squared, n), mean) %>% 
    mutate(r.squared = round(r.squared, 3),
           n = round(n, 0)) %>% 
    gather(term, varnm, r.squared:n) 
  
  rbind(fmcoef, fmrsq)
}

fmcoef <- fm_model(exret ~ beta + lnMe + lnEP_pos + EP_dum + abtover)

fmcoef %>% 
  rename(`(1)` = varnm) %>% 
  kable() %>% 
  kable_styling("striped")
