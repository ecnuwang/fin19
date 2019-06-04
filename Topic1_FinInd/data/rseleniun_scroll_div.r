library(tidyverse)
library(stringr)
library(rvest)
library(RSelenium)
library(lubridate)
library(plotly)

########################

remDr <- remoteDriver(port = 4444L)
remDr$open()

remDr$navigate("http://data.stats.gov.cn/easyquery.htm?cn=E0103")

remDr$findElement("css", "body")$sendKeysToElement(list(key = "end"))

remDr$findElement("css", "[title=国民经济核算]")$clickElement()
remDr$findElement("css", "[title=地区生产总值]")$clickElement()


remDr$findElement("css", "#mySelect_sj")$clickElement()
year <- remDr$findElements("css", "#mySelect_sj .dtList li")
year[[3]]$clickElement()


remDr$findElement("css", "#mySelect_reg")$clickElement()
prov <- remDr$findElements("css", "#mySelect_reg .dtList li")


get_data <- function(prov){
  
  remDr$findElement("css", "#mySelect_reg")$clickElement()
  remDr$executeScript("arguments[0].scrollIntoView(true);", list(prov))
  prov$clickElement()
  
  data <- remDr$getPageSource()[[1]] %>% 
    read_html() %>% 
    html_node(".public_table.table_column") %>% 
    html_table() %>% 
    t() %>% 
    as_tibble(rownames = "year")
  
  colnames(data) <- data[1,]
  
  data[-1,] %>% 
    mutate(
      `指标` = str_extract(`指标`, "\\d{4}")
    ) %>% 
    mutate_all(as.numeric)
  
}


fin_data <- tibble(
  prov_elem = prov[-1],
  prov_name = prov_elem %>% map( ~.$getElementAttribute("title")) %>% simplify_all() %>% flatten_chr(),
  data = prov_elem %>% map(get_data)
)


fin_data %>% 
  unnest(data) %>%
  write_excel_csv("fin_data.csv")
  
##########

popu_data <- read_csv("popu_data.csv")
wage_data <- read_csv("wage_data.csv")
gdp_data <- read_csv("fin_data.csv")

wage_date <- wage_data %>% 
  gather(hangye, wage, -prov_name, -`指标`) %>% 
  rename(year = `指标`)

fin_gdp_pect <- gdp_data %>% 
  mutate(
    fin_gdp_pect = `金融业增加值(亿元)`/`地区生产总值(亿元)`
  ) %>% 
  select(prov_name, year = `指标`, fin_gdp_pect) %>% 
  mutate(
    prov_name = prov_name %>% str_extract("^.{2}"),
    prov_name = if_else(prov_name == "黑龙", "黑龙江", prov_name),
    prov_name = if_else(prov_name == "内蒙", "内蒙古", prov_name)
  )
  
fin_wage <- wage_date %>% 
  spread(hangye, wage) %>% 
  mutate(
    pect = `金融业城镇单位就业人员平均工资(元)`/`城镇单位就业人员平均工资(元)`
  ) %>% 
  select(prov_name, year, pect)
  

p <- ggplot(fin_wage) + geom_line(aes(x = year, y = pect, color = prov_name)) +
  geom_line(aes(x = year, y = pect), size = 2, color = "darkorange", stat = "summary", fun.y = median)

ggplotly(p)

#######################

popu_data <- popu_data %>% 
  gather(hangye, popu, -prov_name, -`指标`)
  
full_data <- wage_date %>% 
  mutate(
    hangye = hangye %>% str_extract(".*(?=平均工资)")
    ) %>% 
  full_join(
    popu_data %>% 
      rename(year = `指标`) %>% 
      mutate(
        hangye = hangye %>% str_extract(".*(?=\\()")
      )
  )


full_data <- full_data %>% 
  filter(hangye %in% c("城镇单位就业人员", "金融业城镇单位就业人员")) %>% 
  drop_na() %>% 
  unite("wage_popu", wage, popu) %>% 
  spread(hangye, wage_popu) %>% 
  separate(`城镇单位就业人员`, c("城镇工资", "城镇人口"), sep = "_") %>% 
  separate(`金融业城镇单位就业人员`, c("金融工资", "金融人口"), sep = "_") %>% 
  type_convert() %>% 
  mutate(
    popu_pect = `金融人口`/`城镇人口`,
    wage_pect = `金融工资`/`城镇工资`
  )

############

prov_location <- tibble(
  `东部沿海地区` = "北京、天津、河北、辽宁、上海、江苏、浙江、福建、山东、广东、广西、海南、重庆、大连、宁波、厦门、青岛、深圳",
  `中部内陆地区` = "山西、内蒙古、吉林、黑龙江、安徽、江西、河南、湖北、湖南",
  `西部边远地区` = "四川、贵州、云南、西藏、陕西、甘肃、青海、宁夏、新疆"
) %>% 
  gather(`经济带划分`, prov_name, everything()) %>% 
  mutate(
    prov_name = prov_name %>% str_split("、")
  ) %>% 
  unnest()

#############

p <- full_data %>% 
  mutate(
    prov_name = prov_name %>% str_extract("^.{2}"),
    prov_name = if_else(prov_name == "黑龙", "黑龙江", prov_name),
    prov_name = if_else(prov_name == "内蒙", "内蒙古", prov_name)
  ) %>% 
  # select(prov_name, year, popu_pect, wage_pect) %>% 
  left_join(fin_gdp_pect) %>% 
  filter(prov_name != "宁夏", prov_name != "西藏") %>% 
  left_join(prov_location) %>% 
  ggplot(aes(x = fin_gdp_pect, y = wage_pect, frame = year)) +  
    geom_point(aes(size = `金融人口`, fill = `经济带划分`), alpha = 0.6, color = "white") +
    geom_text(aes(label = prov_name), check_overlap = TRUE, vjust = 0, nudge_y = 0, size = 3.5) +
    scale_size(range = c(1, 20))

ggplotly(p) %>% animation_opts(frame = 1000)


############




