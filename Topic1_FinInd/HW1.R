# 作图1：A股上市公司行业利润总额

comsum_ind <- comdata_all %>% 
  group_by(year, fin) %>% 
  summarise(tot_profit = sum(profit))

comsum_ind %>% 
  ggplot(aes(x = year, y = tot_profit/1e12)) +
  geom_area(aes(fill = factor(fin, levels = c("1", "0"), labels = c("金融", "其他")))) +
  labs(title = "A股上市企业行业利润总额", subtitle = "1992-2017年", x = "", y = "行业利润总额(万亿元)", caption = "数据来源: CSMAR", fill = "行业") +
  scale_x_continuous(breaks = seq(1992, 2017, 5))

# 作图2：银行业上市企业利息收入占比

comsum_bank <- comdata_all %>% 
  filter(Nindcd == "I01") %>% 
  group_by(year) %>% 
  summarise(`利息收入` = sum(interest), `中间业务收入` = sum(interim)) %>% 
  mutate(
    interest_ratio = `利息收入` / (`利息收入`+`中间业务收入`) * 100
  )

comsum_bank%>% 
  ggplot() +
  geom_line(aes(x = year, y = interest_ratio)) +
  labs(title = "A股银行上市企业利息收入占比", subtitle = "1992-2017年", x = "", y = "利息收入占比(%)", caption = "数据来源: CSMAR") +
  scale_x_continuous(breaks = seq(1992, 2017, 5))
