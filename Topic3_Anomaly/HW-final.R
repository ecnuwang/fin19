# ??????Ҫ??R??
library(tidyverse)
library(readxl)
library(lubridate)
library(TTR)
library(lmtest)
library(sandwich) #NeweyWest
library(broom)
library(stargazer)
library(kableExtra)
library(tibbletime)

#??ȡ??????????
turnover <- read_excel("data/Liq_Tover_M.xlsx")

#?ڶ?????ĸ
turnover1 <- turnover[c(1,2,8)] %>% 
  mutate_if(is.numeric, ~replace_na(ToverTlMAvg, 0)) %>% 
  mutate(
    turnover_avg = SMA (ToverTlMAvg, n = 13) *13/12 - ToverTlMAvg / 12
  )

turnover2 <- turnover %>% 
  left_join(., turnover1)

#?㳬?????
turnover2 <- turnover2 %>% 
  mutate(
    abnormal_turnover = lag(ToverTlMAvg) / turnover_avg
  )
##?ڶ???????

#?ָ?ʱ????��
turnover2 <- turnover2 %>% 
  separate(Trdmnt, into = c("year", "month"), sep = "-") %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month)
  )

#????????
# ??ȡ???ɽ??????ݣ??Ա?ʹ?ø???????ֵ????
trdmnth <- read_excel("Anomaly_data/TRD_Mnth.xlsx")

trdmnth1 <- trdmnth[c(1,2,10)] %>% 
  separate(Trdmnt, into = c("year", "month"), sep = "-") %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month)
  )

turnover2 <- turnover2 %>% 
  left_join(., trdmnth1) 

#????????
turnover3 <- turnover2 %>% 
  group_by(year, month) %>% 
  mutate(
    ab_grp = ntile(abnormal_turnover, 10)
  )  %>% 
  ungroup()

turnover4 <- turnover3 %>% 
  group_by(year, month, ab_grp) %>% 
  summarise(
    vwr = sum(Msmvttl*turnover_avg, na.rm = TRUE) / sum(Msmvttl, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  filter(!is.na(ab_grp)) %>% 
  spread(ab_grp, vwr) %>% 
  mutate(
    rmw = `1` - `10`
  ) 

#?ϲ?????????
factor <- read_excel("Anomaly_data/factor.xlsx")

factor1 <- factor[c(1,2,9,10,11,16,17,18,19,20,21,22,23)]

factor2 <- turnover4 %>% 
  left_join(., factor1)

##??????????????ģ??
#CAPM
capmtest1 <- lm(rmw ~ MKT, factor2)
capmsd1 <- coeftest(capmtest1, vcov = NeweyWest(capmtest1, lag = 4))

# 3-factor model
thftest1 <- lm(rmw ~ MKT + SMB + VMG, factor2)
thfsd1 <- coeftest(thftest1, vcov = NeweyWest(thftest1, lag = 4))

## Fama-French 3-factor model
factor3 <- factor2 %>% 
  filter(!is.na(MKT)) %>% 
  filter(!is.na(HML2)) %>% 
  filter(!is.na(SMB2)) 

thftest1_FF1 <- lm(rmw ~ MKT + SMB2 + HML2, factor3)
thfsd1_FF <- coeftest(thftest1_FF1, vcov = NeweyWest(thftest1_FF1, lag = 4))

stargazer(capmsd1, thfsd1, thfsd1_FF, type = "text",
          intercept.bottom = FALSE, intercept.top = TRUE, 
          column.labels = c("Unconditional Sorts"), 
          column.separate = c(3, 3),
          style = "aer", omit.table.layout = "n")


##??????
#??ȡ?޷????????ʣ????㳬???ر????ã?
# ?޷???????
nrrate <- read_excel("Anomaly_data/TRD_Nrrate.xlsx")

nrrate <- nrrate %>% 
  mutate(
    Clsdt = ymd(Clsdt),
    year = year(Clsdt),
    month = month(Clsdt)
  ) %>% 
  group_by(year, month) %>% 
  summarise(rf = mean(Nrrmtdt))


#?ر???
trdmnth2 <- trdmnth[c(1,2,13)] %>% 
  separate(Trdmnt, into = c("year", "month"), sep = "-") %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month)
  )

data1 <- turnover3 %>% 
  left_join(., nrrate) %>% 
  left_join(., trdmnth2) %>% 
  mutate(
    excess_return = Mretnd - rf
  ) 

#??ȡbeta
stk_beta <- read_csv("Anomaly_data/beta.csv")

data_fm <- data1 %>% 
  left_join(stk_beta, by = c("Stkcd" = "Stkcd", "year" = "year", "month" = "month")) %>% 
  rename(beta = estimate) 
 
#??EP
#??ȡ???????ݣ????þ?????
comfin <- read_excel("Anomaly_data/FI_T2.xlsx", col_names = c("Stkcd", "Accper", "Indcd", "nr", "npexnr", "ROE", "ROE_exnr"), skip = 3)

comfin <- comfin %>% 
  mutate(
    Accper = ymd(Accper),
    year = year(Accper),
    month = month(Accper)
  ) %>% 
  group_by(Stkcd) %>% 
  mutate(
    npexnr = npexnr %>% accumulate(~if_else(is.na(.y), .x, .y))
  )

comfin1 <- comfin[c(1,5,8,9)]

data_fm <- data_fm %>% 
  left_join(., comfin1) %>% 
  mutate(EP = npexnr / Msmvttl)  
  

##fm_model
data_fm3 <- data_fm %>% 
  filter(!is.na( beta)) %>% 
  filter(!is.na( Msmvttl)) %>% 
  filter(!is.na( EP)) %>% 
  filter(!is.na( abnormal_turnover)) 

fm_model <- function(model){
  fmcoef <- data_fm3 %>%
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
}

fmcoef <- fm_model(excess_return ~ beta + Msmvttl + EP + abnormal_turnover)

fmcoef %>% 
  kable() %>% 
  kable_styling("striped")


