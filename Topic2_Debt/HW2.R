library(tidyverse)
library(readxl)
library(lubridate)
library(stargazer)

# 基本信息
cominfo <- read_excel("Zombie_data/TRD_Co.xlsx")

# 资产负债表
combas <- read_excel("Zombie_data/FS_Combas.xlsx", col_names = c("Stkcd", "Accper", "Typrep", "asset", "shortloan", "CP", "longloan", "bond"), skip = 3)

# 利润表
comins <- read_excel("Zombie_data/FS_Comins.xlsx", col_names = c("Stkcd", "Accper", "Typrep", "sales", "profit", "npr"), skip = 3)

# 财务报表附注-财务费用
finexp <- read_excel("Zombie_data/FN_Fn051.xlsx")

# 财务报表附注-政府补助
subsidy <- read_excel("Zombie_data/FN_Fn056.xlsx")

# 利率信息
loan_intrate <- read_csv("Zombie_data/Loan_interest.csv")
bond_intrate <- read_csv("Zombie_data/Bond_interest.csv")


# 初步筛选
combas <- combas %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper)
  ) %>% 
  filter(month(Accper) == 12, Typrep == "A")

comins <- comins %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper)
  ) %>% 
  filter(month(Accper) == 12, Typrep == "A")

# A股非金融企业
com_A <- cominfo %>% 
  filter(
    (Markettype == 1 | Markettype == 4 | Markettype == 16), 
    Indcd != "0001"
  )

# 利息支出
intexp <- finexp %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper)
  ) %>% 
  filter(
    month(Accper) == 12, 
    Typrep == 1,
    Sgnyea == 1,
    Fn05101 == 1,
    !is.na(Fn05102)
  ) %>% 
  rename(intexp = Fn05102) %>% 
  select(Stkcd, year, intexp)

# 计入营业外收入的政府补助
subsidy <- subsidy %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper)
  ) %>% 
  filter(
    month(Accper) == 12, 
    DataSources == "02", 
    Typrep == 1, 
    Fn05601 == "合计",
    !is.na(Fn05602)
  ) %>% 
  rename(subsidy = Fn05602) %>% 
  select(Stkcd, year, subsidy)

# 合并数据集
coms <- combas %>% 
  full_join(comins) %>% 
  left_join(intexp) %>% 
  left_join(subsidy)

comdata <- com_A %>% 
  left_join(coms) %>% 
  filter(asset > 0, sales >= 0) %>% 
  left_join(loan_intrate) %>% 
  left_join(bond_intrate)

# 处理确实值
comdata <- comdata %>%
  mutate_at(
    c("shortloan", "CP", "longloan", "bond", "intexp", "subsidy"), 
    ~replace_na(., 0)
  )

# 构造变量
comdata_all <- comdata %>% 
  mutate(
    interest = 0.9*short_rate/100*shortloan + 0.9*roll_long_rate/100*longloan + bond_rate/100*bond,
    debt = shortloan + CP + longloan + bond
  ) %>% 
  complete(Stkcd, year) %>% 
  arrange(Stkcd, year) %>% 
  group_by(Stkcd) %>% 
  mutate(
    minint = lag(interest),
    int_sub = if_else(intexp < minint, minint - intexp, 0),
    exprofit = profit - subsidy - int_sub,
    exprofit1 = npr,
    debt_lag = lag(debt),
    exprofit_lag = lag(exprofit),
    exprofit_lead = lead(exprofit),
    sumprofit = exprofit_lead + exprofit + exprofit_lag,
    sumprofit1 = lead(exprofit1) + exprofit1 + lag(exprofit1)
  )

# 识别僵尸企业
## CHK方法
com_CHK <- comdata_all %>% 
  filter(year >= 2003, year <= 2015) %>% 
  mutate(
    zombie = case_when(
      is.na(minint) ~ NA_real_,
      intexp < minint ~ 1,
      TRUE ~ 0
    )
  )

## FN方法
com_FN <- comdata_all %>% 
  filter(year >= 2003, year <= 2015) %>% 
  mutate(
    zombie = case_when(
      is.na(minint) ~ NA_real_,
      profit >= minint ~ 0,
      (profit < minint) & (debt_lag > 0.5*asset) & (debt > debt_lag) ~ 1,
      TRUE ~ 0
    )
  )

## 改进方法
com_zombie <- comdata_all %>% 
  mutate(
    sumprofit_lead = lead(sumprofit),
    sumprofit_lag = lag(sumprofit)
  ) %>% 
  filter(year >= 2003, year <= 2015) %>% 
  mutate(
    zombie = case_when(
      is.na(sumprofit_lead) | is.na(sumprofit) | is.na(sumprofit_lag) ~ NA_real_,
      (sumprofit_lead < 0) | (sumprofit < 0) | (sumprofit_lag < 0) ~ 1,
      TRUE ~ 0
    )
  )

## HW2-2 仅考虑净利润
com_zombie1 <- comdata_all %>% 
  mutate(
    sumprofit1_lead = lead(sumprofit1),
    sumprofit1_lag = lag(sumprofit1)
  ) %>% 
  filter(year >= 2003, year <= 2015) %>% 
  mutate(
    zombie = case_when(
      is.na(sumprofit1_lead) | is.na(sumprofit1) | is.na(sumprofit1_lag) ~ NA_real_,
      (sumprofit1_lead < 0) | (sumprofit1 < 0) | (sumprofit1_lag < 0) ~ 1,
      TRUE ~ 0
    )
  )


# 描述性统计
## CHK方法
com_CHK %>%
  filter(!is.na(zombie), !is.na(asset)) %>% 
  group_by(Nnindcd) %>% 
  summarise(
    n = n(),
    z = sum(zombie), 
    pct = mean(zombie), 
    sumat = sum(zombie*asset), 
    pctat = weighted.mean(zombie, asset)
  )

## FN方法
com_FN %>%
  filter(!is.na(zombie), !is.na(asset)) %>% 
  group_by(Nnindcd) %>% 
  summarise(
    n = n(),
    z = sum(zombie), 
    pct = mean(zombie), 
    sumat = sum(zombie*asset), 
    pctat = weighted.mean(zombie, asset)
  )

## 改进方法
com_zombie %>% 
  filter(!is.na(zombie), !is.na(asset)) %>% 
  group_by(Nnindcd) %>% 
  summarise(
    n = n(),
    z = sum(zombie), 
    pct = mean(zombie), 
    sumat = sum(zombie*asset), 
    pctat = weighted.mean(zombie, asset)
  )

## HW2-2 仅考虑政府补贴
### 分年份
com_zombie1 %>% 
  filter(!is.na(zombie), !is.na(asset)) %>% 
  group_by(year) %>% 
  summarise(
    n = n(),
    z = sum(zombie), 
    pct = mean(zombie), 
    sumat = sum(zombie*asset), 
    pctat = weighted.mean(zombie, asset)
  )

### 分行业
com_zombie1 %>% 
  filter(!is.na(zombie), !is.na(asset)) %>% 
  group_by(Nnindcd) %>% 
  summarise(
    n = n(),
    z = sum(zombie), 
    pct = mean(zombie), 
    sumat = sum(zombie*asset), 
    pctat = weighted.mean(zombie, asset)
  )
