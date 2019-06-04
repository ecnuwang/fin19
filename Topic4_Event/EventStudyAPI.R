library(tidyverse)
library(readxl)
library(lubridate)
library(EventStudy)

cominfo <- read_excel("Event_data/TRD_Co.xlsx")

comA <- cominfo %>% 
  filter(
    Indcd != "0001",
    Markettype == 1 | Markettype == 4 | Markettype == 16
  )

read_file <- function(df){
  read_excel(df) %>%
    filter(Markettype == 1 | Markettype == 4 | Markettype == 16) %>%
    select(Stkcd, Trddt, Clsprc, Adjprcwd, Dsmvtll) %>%
    mutate(
      Trddt = ymd(Trddt)
    )
}

# read and merge with factors and risk-free rate
stkdaily <- tibble(
  name = dir(path = "D:/Qsync/Courses/金融研究方法/2019本科/上机课/T4_Event/Event_data", pattern = "TRD_Dalyr0")
) %>%
  mutate(
    f = str_c("D:/Qsync/Courses/金融研究方法/2019本科/上机课/T4_Event/Event_data", name, sep = "/"),
    data = f %>% map(read_file)
  ) %>%
  unnest(data) %>%
  select(-name, -f) %>%
  semi_join(comA)

fac_d <- stkdaily %>% 
  group_by(Trddt) %>% 
  summarise(`Closing Price` = weighted.mean(Clsprc*Adjprcwd, Dsmvtll, na.rm = T)) %>% 
  rename(Date = Trddt) %>% 
  mutate(`Market ID` = "Index300") %>% 
  write_excel_csv("data/03_marketData.csv")
  
stk_d <- stkdaily %>% 
  select(Stkcd, Trddt, Clsprc) %>% 
  rename(`Firm ID` = Stkcd,
         Date = Trddt,
         `Closing Price` = Clsprc) %>% 
  write_excel_csv("data/02_firmData.csv")

RRR <- read_excel("Event_data/RRR.xlsx") %>% 
  mutate(AnnDate = ymd(AnnDate),
         `Event ID` = row_number()) %>% 
  left_join(stkdaily, by = c("AnnDate" = "Trddt")) %>% 
  select(`Event ID`, Stkcd, AnnDate) %>% 
  rename(`Firm ID` = Stkcd,
         `Event Date` = AnnDate) %>% 
  mutate(`Market ID` = "Index300",
         `Grouping Variable` = "all",
         `Start Event Window` = -1,
         `End Event Window` = 1,
         `End of Estimation Window` = -5,
         `Estimation Window Length` = 120) %>% 
  write_excel_csv("data/01_RequestFile.csv")

#####################

apiUrl <- "http://api.eventstudytools.com"
apiKey <- "573e58c665fcc08cc6e5a660beaad0cb"

# Setup API Connection
estSetup <- EventStudyAPI$new(apiUrl)
estSetup$authentication(apiKey)

# Perform Abnormal Return Event Study
returnEstParams <- ARCApplicationInput$new()

returnEstParams$setTestStatistics(c("caarptlz"))

returnEstParams$setResultFileType("xlsx")

returnEstParams$setReturnType("simple")

returnEstParams$setNonTradingDays("skip")

returnEstParams$setBenchmarkModel("mm")

estResults <- estSetup$performEventStudy(estParams     = returnEstParams, 
                                         dataFiles     = c("request_file" = "data/01_RequestFile.csv",
                                                           "firm_data"    = "data/02_firmData.csv",
                                                           "market_data"  = "data/03_marketData.csv"),
                                         downloadFiles = T, # download result files
                                         checkFiles    = T) # check input files 




