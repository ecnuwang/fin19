library(tidyverse)
library(readxl)
library(lubridate)
library(broom)
library(kableExtra)
library(stargazer)
library(lmtest) #coeftest
library(sandwich) #NeweyWest

cominfo <- read_excel("Event_data/TRD_Co.xlsx")

location <- read_excel("data/location.xlsx") %>%
  mutate(
    Stkcd = str_sub(`证券代码`, 1, 6),
    Prvcnm = if_else(str_count(`省份`) == 2, str_c(`省份`, "市"), `省份`)
  ) %>% 
  select(Stkcd, Prvcnm)

comA <- cominfo %>%
  left_join(location) %>% 
  filter(
    Indcd != "0001",
    Markettype == 1 | Markettype == 4 | Markettype == 16
  )

# nrrate_d <- read_excel("Event_data/TRD_Nrrate.xlsx") %>%
#   mutate(rf = Nrrdaydt/100,
#          Clsdt = ymd(Clsdt)) %>%
#   select(Clsdt, rf)
# 
# fac_d <- read_excel("Event_data/STK_MKT_ThrfacDay.xlsx") %>%
#   filter(MarkettypeID == "P9709") %>%
#   mutate(
#     TradingDate = ymd(TradingDate)
#   ) %>%
#   select(TradingDate, RiskPremium2)
# 
# 
# # read stock daily return
# read_file <- function(df){
#   read_excel(df) %>%
#     filter(Markettype == 1 | Markettype == 4 | Markettype == 16) %>%
#     select(Stkcd, Trddt, Dretwd, Dsmvtll) %>%
#     mutate(
#       Trddt = ymd(Trddt)
#     )
# }
# 
# # read and merge with factors and risk-free rate
# stkdaily <- tibble(
#   name = dir(path = "D:/Qsync/Courses/金融研究方法/2019本科/上机课/T4_Event/Event_data", pattern = "TRD_Dalyr0")
# ) %>%
#   mutate(
#     f = str_c("D:/Qsync/Courses/金融研究方法/2019本科/上机课/T4_Event/Event_data", name, sep = "/"),
#     data = f %>% map(read_file)
#   ) %>%
#   unnest(data) %>%
#   select(-name, -f) %>%
#   left_join(nrrate_d, by = c("Trddt" = "Clsdt")) %>%
#   left_join(fac_d, by = c("Trddt" = "TradingDate")) %>%
#   mutate(exret = Dretwd - rf) %>% 
#   semi_join(comA)
  
stkdaily <- read_rds("D:/Qsync/Courses/金融研究方法/2019本科/上机课/T4_Event/Event_data/stkdaily.rds")

stk_event <- stkdaily %>% 
  filter(!is.na(exret),
         Trddt >= as.Date("2012-3-13"), 
         Trddt <= as.Date("2012-3-15")) %>% 
  mutate(
    event = case_when(Trddt == as.Date("2012-3-14") ~ 0,
                      Trddt == as.Date("2012-3-13") ~ -1,
                      TRUE ~ 1)
  ) %>% 
  group_by(Stkcd) %>% 
  filter(n() >= 1)

stk_est <- stkdaily %>% 
  filter(!is.na(exret),
         Trddt >= as.Date("2011-8-7"), 
         Trddt <= as.Date("2012-2-6")) %>% 
  mutate(event = -99) %>% 
  group_by(Stkcd) %>% 
  filter(n() >= 100) %>% 
  semi_join(stk_event, by = c("Stkcd" = "Stkcd"))

stk_event <- stk_event %>% 
  semi_join(stk_est, by = c("Stkcd" = "Stkcd"))

eventdata <- bind_rows(stk_event, stk_est) %>% 
  arrange(Stkcd, Trddt)

# Market model
AR <- eventdata %>%
  select(Stkcd, event, exret, RiskPremium2) %>% 
  nest(-Stkcd) %>% 
  mutate(
    model = data %>% map(filter, event == -99) %>% map(~lm(exret ~ RiskPremium2, .)),
    pred_data = data %>% map(filter, event >= -1, event <= 1),
    pred = map2(model, pred_data, predict),
    event = data %>% map(filter, event >= -1, event <= 1) %>% map("event"),
    ret = data %>% map(filter, event >= -1, event <= 1) %>% map("exret")
  ) %>% 
  select(Stkcd, event, pred, ret) %>% 
  unnest() %>% 
  mutate(
    ret = as.numeric(ret),
    AR = ret - pred
  )

CAR <- AR %>% 
  group_by(Stkcd) %>% 
  mutate(
    CAR = cumsum(AR)
  ) 

## Fixed investment
inv <- read_excel("data/CRE_Ifa02.xlsx") %>% 
  filter(Sgnyea == 2010 | Sgnyea == 2011 | Sgnyea == 2012) %>% 
  mutate(invr = Ifa0202 / Ifa0201) %>% 
  group_by(Prvcnm) %>% 
  summarise(invr = mean(invr)) %>% 
  mutate(grp = ntile(invr, 3)) 

# aggregate into portfolio
event_grp <- eventdata %>%
  left_join(comA) %>% 
  left_join(inv) %>% 
  filter(!is.na(grp)) %>% 
  group_by(grp, Trddt, event) %>% 
  summarise_at(vars(exret, RiskPremium2), mean, na.rm = T) %>% 
  select(grp, event, exret, RiskPremium2) 

AR_grp <- event_grp %>% 
  group_by(grp) %>% 
  nest() %>% 
  mutate(
    model = data %>% map(filter, event == -99) %>% map(~lm(exret ~ RiskPremium2, .)),
    pred_data = data %>% map(filter, event >= -1, event <= 1),
    pred = map2(model, pred_data, predict),
    event = data %>% map(filter, event >= -1, event <= 1) %>% map("event"),
    ret = data %>% map(filter, event >= -1, event <= 1) %>% map("exret")
  ) %>% 
  select(grp, event, pred, ret) %>% 
  unnest() %>% 
  mutate(
    ret = as.numeric(ret),
    AR = ret - pred
  ) 

CAR_grp <- AR_grp %>% 
  group_by(grp) %>% 
  mutate(
    CAR = cumsum(AR)
  ) %>% 
  filter(event == 1)

# portfolio method
CAR_Var1 <- event_grp %>% 
  group_by(grp) %>% 
  nest() %>% 
  mutate(
    model_data = data %>% map(filter, event == -99),
    pred_data = data %>% map(filter, event >= -1, event <= 1),
    model = model_data %>% map(~lm(exret ~ RiskPremium2, .)),
    sigma = model %>% map(summary) %>% map_dbl("sigma"),
    X = model_data %>% map(~model.matrix(exret ~ RiskPremium2, .)),
    X_star = pred_data %>%  map(~model.matrix(exret ~ RiskPremium2, .)),
    dim_star = X_star %>% map(dim) %>% map_int(1),
    V = pmap(list(dim_star, X, X_star, sigma),
             function(dim_star, X, X_star, sigma)
               (diag(dim_star) + X_star %*% solve(t(X) %*% X) %*% t(X_star)) * sigma^2
    ),
    Var_CAR = V %>% map_dbl(sum)
  ) %>%
  select(grp, Var_CAR) %>%
  unnest()

# method used in the paper
CAR_Var2 <- eventdata %>%
  select(Stkcd, Trddt, event, exret, RiskPremium2) %>%
  nest(-Stkcd) %>%
  mutate(
    model_data = data %>% map(filter, event == -99),
    resid = model_data %>% map(~lm(exret ~ RiskPremium2, .)) %>% map(augment),
    Trddt = data %>% map(filter, event == -99) %>% map("Trddt")
  ) %>%
  select(Stkcd, Trddt, resid) %>%
  unnest() %>% 
  left_join(comA) %>% 
  left_join(inv) %>% 
  filter(!is.na(grp)) %>% 
  group_by(grp, Trddt) %>% 
  summarise(AR = mean(`.resid`, na.rm = T)) %>% 
  group_by(grp) %>% 
  summarise(std_CAR = sqrt(3)*sd(AR, na.rm = T))

CAR_est <- CAR_grp %>%
  left_join(CAR_Var1) %>%
  left_join(CAR_Var2) %>% 
  mutate(
    J1_1 = CAR/sqrt(Var_CAR),
    J1_2 = CAR/std_CAR
  )

CAR_est %>%
  select(grp, CAR, Var_CAR, J1_1, std_CAR, J1_2) %>% 
  kable(digits = 5) %>% 
  kable_styling("striped")


# regression
# calculate firm-level controls: size, bm, leverage
combas <- read_excel("Event_data/FS_Combas.xlsx", col_names = c("Stkcd", "Accper", "Typrep", "Asset", "Liability", "Equity"), skip = 3) %>% 
  mutate(Accper = ymd(Accper)) %>% 
  filter(year(Accper) == 2011, month(Accper) == 12, Typrep == "A")

CAR_control <- stkdaily %>% 
  filter(Trddt == as.Date("2012-3-7")) %>% 
  left_join(combas) %>% 
  mutate(
    lnSZ = log(Dsmvtll),
    BM = Equity / Dsmvtll,
    leverage = Liability / Asset
  ) %>% 
  select(Stkcd, lnSZ, BM, leverage)

# winsorize at 0.5% and 99.5%
CAR_control <- CAR_control %>%
  mutate_if(is.numeric,
            ~ case_when(
              . < quantile(., 0.005, na.rm = T) ~ quantile(., 0.005, na.rm = T),
              . > quantile(., 0.995, na.rm = T) ~ quantile(., 0.995, na.rm = T),
              TRUE ~ as.numeric(.)),
            na.rm = TRUE)

# merge with CAR and sensitivity data
CAR_reg <- CAR %>%
  filter(event == 1) %>% 
  left_join(comA) %>% 
  left_join(inv) %>% 
  left_join(CAR_control) %>% 
  mutate(
    lnBM = if_else(BM < 0, NA_real_, log(BM))
  )

# run regression
CARreg <- lm(CAR ~ invr + lnSZ + lnBM + leverage, CAR_reg)

# show results
stargazer(CARreg, type = "text", 
          style = "aer", omit.table.layout = "n", omit.stat = c("rsq", "ser"))
