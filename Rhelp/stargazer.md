---
title: "A Stargazer Cheatsheet"
date: "Updated 27 三月, 2019"
output: 
  html_document:
    theme: spacelab
    keep_md: true
    toc: yes
    highlight: haddock
---

## Dataset: dplyr and nycflights13

Setting up a dataset for this cheatsheet allows me to spotlight two recent
R packages created by [Hadley Wickham](http://had.co.nz/). The first, 
[dplyr](https://github.com/hadley/dplyr), is a set of new tools for data 
manipulation. Using `dplyr`, I will extract flights and weather data from 
another new package called 
[nycflights13](https://github.com/hadley/nycflights13). With this data I 
will show how to estimate a couple of regression models and nicely format the 
results into tables with `stargazer`.

Note: `stargazer` v. 5.1 does not play nicely with `dplyr`'s tbl_df class.
As a temporary work-around I pipe the merged dataset to `data.frame`.


```r
library("dplyr")
library("nycflights13")
library("AER") # Applied Econometrics with R
library("stargazer")
daily <- flights %>%
  filter(origin == "EWR") %>%
  group_by(year, month, day) %>%
  summarise(delay = mean(dep_delay, na.rm = TRUE))
daily_weather <- weather %>%
  filter(origin == "EWR") %>%
  group_by(year, month, day) %>%
  summarise(temp   = mean(temp, na.rm = TRUE),
            wind   = mean(wind_speed, na.rm = TRUE),
            precip = sum(precip, na.rm = TRUE))
# Merge flights with weather data frames
both <- inner_join(daily, daily_weather) %>% 
  data.frame()  # Temporary fix
# Create an indicator for quarter
both$quarter <- cut(both$month, breaks = c(0, 3, 6, 9, 12), 
                                labels = c("1", "2", "3", "4"))
# Create a vector of class logical
both$hot <- as.logical(both$temp > 85)
head(both)
```

```
##   year month day     delay     temp      wind precip quarter   hot
## 1 2013     1   1 17.483553 36.81909 13.233970      0       1 FALSE
## 2 2013     1   2 25.322674 28.70000 10.884461      0       1 FALSE
## 3 2013     1   3  8.450450 29.57750  8.582901      0       1 FALSE
## 4 2013     1   4 12.103858 34.33250 14.001157      0       1 FALSE
## 5 2013     1   5  5.696203 36.56000  9.398037      0       1 FALSE
## 6 2013     1   6 12.383333 39.92000  9.110342      0       1 FALSE
```

\
We can use the `both` data frame to estimate a couple of linear models 
predicting the average delay out of Newark controlling for the weather. The 
first model will use only the weather variables and in the second I'll add 
dummy variables indicating the quarter. I also estimate a third model, using 
using the `ivreg` command from package 
[AER](http://cran.r-project.org/web/packages/AER/index.html) to demonstrate 
output with mixed models. The raw R output:


```r
output  <- lm(delay ~ temp + wind + precip, data = both)
output2 <- lm(delay ~ temp + wind + precip + quarter, data = both)
# Instrumental variables model 
output3 <- ivreg(delay ~ temp + wind + precip | . - temp + hot, data = both)
summary(output)
```

```
## 
## Call:
## lm(formula = delay ~ temp + wind + precip, data = both)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -47.164  -8.150  -3.148   4.667  70.493 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  6.15166    2.98897   2.058   0.0403 *  
## temp         0.07905    0.03933   2.010   0.0452 *  
## wind         0.28523    0.15741   1.812   0.0708 .  
## precip      15.25495    2.01280   7.579 2.98e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.82 on 360 degrees of freedom
## Multiple R-squared:  0.1541,	Adjusted R-squared:  0.147 
## F-statistic: 21.86 on 3 and 360 DF,  p-value: 5.003e-13
```

```r
summary(output2)
```

```
## 
## Call:
## lm(formula = delay ~ temp + wind + precip + quarter, data = both)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -47.458  -8.086  -3.294   4.542  69.426 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  5.24570    3.41233   1.537   0.1251    
## temp         0.16433    0.06628   2.479   0.0136 *  
## wind         0.23710    0.15733   1.507   0.1327    
## precip      14.84493    2.00409   7.407 9.38e-13 ***
## quarter2    -2.06050    2.57919  -0.799   0.4249    
## quarter3    -6.75538    3.12100  -2.164   0.0311 *  
## quarter4    -4.40842    2.03707  -2.164   0.0311 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.72 on 357 degrees of freedom
## Multiple R-squared:  0.1747,	Adjusted R-squared:  0.1608 
## F-statistic: 12.59 on 6 and 357 DF,  p-value: 6.823e-13
```

```r
summary(output3)
```

```
## 
## Call:
## ivreg(formula = delay ~ temp + wind + precip | . - temp + hot, 
##     data = both)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -47.100  -8.386  -3.451   4.613  69.801 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   8.5765     8.7968   0.975    0.330    
## temp          0.0412     0.1350   0.305    0.760    
## wind          0.2501     0.1980   1.263    0.207    
## precip       15.3351     2.0338   7.540 3.86e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.84 on 360 degrees of freedom
## Multiple R-Squared: 0.1519,	Adjusted R-squared: 0.1449 
## Wald test: 20.49 on 3 and 360 DF,  p-value: 2.788e-12
```

\
[Back to table of contents](#TOC)

## Quick notes

Since I'm using [knitr](http://cran.rstudio.com/web/packages/knitr/index.html) 
and [R markdown](http://rmarkdown.rstudio.com/) to create this webpage, in the 
code that follows I will include the `stargazer` option `type = "html"`. 
`stargazer` is set to produce **LaTeX** output by default. If you desire 
**LaTeX** output, just remove the type option from the code below.

\
Also, while I have added an example for many of the available `stargazer` 
options, I have not included all of them. So while you're likely to find a 
relevant example don't assume if it's not listed that `stargazer` can't do 
it. Check the documentation for additional features and updates to the package.
It is often the case that an omitted argument is specific for **LaTeX** output 
and I can't demonstrate it here.

### HTML formatting

It is possible to change the formatting of html tables generated with `stargazer` 
via an html style sheet. See the [R Markdown documentation](http://rmarkdown.rstudio.com/html_document_format.html#custom_css) 
about incorporating a custom CSS. 

\
[Back to table of contents](#TOC)

## The default summary statistics table

```r
stargazer(both, type = "html")
```


<table style="text-align:center"><tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>N</td><td>Mean</td><td>St. Dev.</td><td>Min</td><td>Pctl(25)</td><td>Pctl(75)</td><td>Max</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">year</td><td>364</td><td>2,013.000</td><td>0.000</td><td>2,013</td><td>2,013</td><td>2,013</td><td>2,013</td></tr>
<tr><td style="text-align:left">month</td><td>364</td><td>6.511</td><td>3.445</td><td>1</td><td>4</td><td>9.2</td><td>12</td></tr>
<tr><td style="text-align:left">day</td><td>364</td><td>15.679</td><td>8.784</td><td>1</td><td>8</td><td>23</td><td>31</td></tr>
<tr><td style="text-align:left">delay</td><td>364</td><td>15.080</td><td>13.883</td><td>-1.349</td><td>5.446</td><td>20.007</td><td>97.771</td></tr>
<tr><td style="text-align:left">temp</td><td>364</td><td>55.531</td><td>17.598</td><td>15.455</td><td>40.114</td><td>71.647</td><td>91.085</td></tr>
<tr><td style="text-align:left">wind</td><td>364</td><td>9.464</td><td>4.398</td><td>2.637</td><td>6.797</td><td>11.508</td><td>56.388</td></tr>
<tr><td style="text-align:left">precip</td><td>364</td><td>0.121</td><td>0.335</td><td>0</td><td>0</td><td>0.03</td><td>4</td></tr>
<tr><td style="text-align:left">hot</td><td>364</td><td>0.022</td><td>0.147</td><td>0</td><td>0</td><td>0</td><td>1</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>

\
When supplied a data frame, by default `stargazer` creates a table with summary 
statistics. If the `summary` option is set to `FALSE` then stargazer 
will instead print the contents of the data frame.


```r
# Use only a few rows
both2 <- both %>% slice(1:6)
stargazer(both2, type = "html", summary = FALSE, rownames = FALSE)
```


<table style="text-align:center"><tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">year</td><td>month</td><td>day</td><td>delay</td><td>temp</td><td>wind</td><td>precip</td><td>quarter</td><td>hot</td></tr>
<tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">2,013</td><td>1</td><td>1</td><td>17.484</td><td>36.819</td><td>13.234</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>2</td><td>25.323</td><td>28.700</td><td>10.884</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>3</td><td>8.450</td><td>29.578</td><td>8.583</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>4</td><td>12.104</td><td>34.333</td><td>14.001</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>5</td><td>5.696</td><td>36.560</td><td>9.398</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>6</td><td>12.383</td><td>39.920</td><td>9.110</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr></table>

### Remove row and column names


```r
stargazer(both2, type = "html", summary = FALSE,
          rownames = FALSE,
          colnames = FALSE)
```


<table style="text-align:center"><tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">2,013</td><td>1</td><td>1</td><td>17.484</td><td>36.819</td><td>13.234</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>2</td><td>25.323</td><td>28.700</td><td>10.884</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>3</td><td>8.450</td><td>29.578</td><td>8.583</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>4</td><td>12.104</td><td>34.333</td><td>14.001</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>5</td><td>5.696</td><td>36.560</td><td>9.398</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td style="text-align:left">2,013</td><td>1</td><td>6</td><td>12.383</td><td>39.920</td><td>9.110</td><td>0</td><td>1</td><td>FALSE</td></tr>
<tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr></table>

### Change which statistics are displayed

In order to customize which summary statistics are displayed, change any of the
the following (logical) parameters, `nobs`, `mean.sd`, `min.max`, `median`, and
`iqr`. 


```r
stargazer(both, type = "html", nobs = FALSE, mean.sd = TRUE, median = TRUE,
          iqr = TRUE)
```


<table style="text-align:center"><tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>Mean</td><td>St. Dev.</td><td>Min</td><td>Pctl(25)</td><td>Median</td><td>Pctl(75)</td><td>Max</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">year</td><td>2,013.000</td><td>0.000</td><td>2,013</td><td>2,013</td><td>2,013</td><td>2,013</td><td>2,013</td></tr>
<tr><td style="text-align:left">month</td><td>6.511</td><td>3.445</td><td>1</td><td>4</td><td>7</td><td>9.2</td><td>12</td></tr>
<tr><td style="text-align:left">day</td><td>15.679</td><td>8.784</td><td>1</td><td>8</td><td>16</td><td>23</td><td>31</td></tr>
<tr><td style="text-align:left">delay</td><td>15.080</td><td>13.883</td><td>-1.349</td><td>5.446</td><td>10.501</td><td>20.007</td><td>97.771</td></tr>
<tr><td style="text-align:left">temp</td><td>55.531</td><td>17.598</td><td>15.455</td><td>40.114</td><td>56.356</td><td>71.647</td><td>91.085</td></tr>
<tr><td style="text-align:left">wind</td><td>9.464</td><td>4.398</td><td>2.637</td><td>6.797</td><td>8.895</td><td>11.508</td><td>56.388</td></tr>
<tr><td style="text-align:left">precip</td><td>0.121</td><td>0.335</td><td>0</td><td>0</td><td>0</td><td>0.03</td><td>4</td></tr>
<tr><td style="text-align:left">hot</td><td>0.022</td><td>0.147</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>

### Change which statistics are displayed (a second way)


```r
stargazer(both, type = "html", summary.stat = c("n", "p75", "sd"))
```


<table style="text-align:center"><tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>N</td><td>Pctl(75)</td><td>St. Dev.</td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">year</td><td>364</td><td>2,013</td><td>0.000</td></tr>
<tr><td style="text-align:left">month</td><td>364</td><td>9.2</td><td>3.445</td></tr>
<tr><td style="text-align:left">day</td><td>364</td><td>23</td><td>8.784</td></tr>
<tr><td style="text-align:left">delay</td><td>364</td><td>20.007</td><td>13.883</td></tr>
<tr><td style="text-align:left">temp</td><td>364</td><td>71.647</td><td>17.598</td></tr>
<tr><td style="text-align:left">wind</td><td>364</td><td>11.508</td><td>4.398</td></tr>
<tr><td style="text-align:left">precip</td><td>364</td><td>0.03</td><td>0.335</td></tr>
<tr><td style="text-align:left">hot</td><td>364</td><td>0</td><td>0.147</td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr></table>

### Remove logical variables in the summary statistics

`stargazer` reports summary statistics for logical variables by default 
(0 = FALSE and 1 = TRUE). To supress the reporting of logical vectors change 
`summary.logical` to `FALSE`. Note the stats for our vector `hot` are gone. 


```r
stargazer(both, type = "html", summary.logical = FALSE)
```


<table style="text-align:center"><tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>N</td><td>Mean</td><td>St. Dev.</td><td>Min</td><td>Pctl(25)</td><td>Pctl(75)</td><td>Max</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">year</td><td>364</td><td>2,013.000</td><td>0.000</td><td>2,013</td><td>2,013</td><td>2,013</td><td>2,013</td></tr>
<tr><td style="text-align:left">month</td><td>364</td><td>6.511</td><td>3.445</td><td>1</td><td>4</td><td>9.2</td><td>12</td></tr>
<tr><td style="text-align:left">day</td><td>364</td><td>15.679</td><td>8.784</td><td>1</td><td>8</td><td>23</td><td>31</td></tr>
<tr><td style="text-align:left">delay</td><td>364</td><td>15.080</td><td>13.883</td><td>-1.349</td><td>5.446</td><td>20.007</td><td>97.771</td></tr>
<tr><td style="text-align:left">temp</td><td>364</td><td>55.531</td><td>17.598</td><td>15.455</td><td>40.114</td><td>71.647</td><td>91.085</td></tr>
<tr><td style="text-align:left">wind</td><td>364</td><td>9.464</td><td>4.398</td><td>2.637</td><td>6.797</td><td>11.508</td><td>56.388</td></tr>
<tr><td style="text-align:left">precip</td><td>364</td><td>0.121</td><td>0.335</td><td>0</td><td>0</td><td>0.03</td><td>4</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>

### Flip the table axes


```r
stargazer(both, type = "html", flip = TRUE)
```


<table style="text-align:center"><tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>year</td><td>month</td><td>day</td><td>delay</td><td>temp</td><td>wind</td><td>precip</td><td>hot</td></tr>
<tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">N</td><td>364</td><td>364</td><td>364</td><td>364</td><td>364</td><td>364</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">Mean</td><td>2,013.000</td><td>6.511</td><td>15.679</td><td>15.080</td><td>55.531</td><td>9.464</td><td>0.121</td><td>0.022</td></tr>
<tr><td style="text-align:left">St. Dev.</td><td>0.000</td><td>3.445</td><td>8.784</td><td>13.883</td><td>17.598</td><td>4.398</td><td>0.335</td><td>0.147</td></tr>
<tr><td style="text-align:left">Min</td><td>2,013</td><td>1</td><td>1</td><td>-1.349</td><td>15.455</td><td>2.637</td><td>0</td><td>0</td></tr>
<tr><td style="text-align:left">Pctl(25)</td><td>2,013</td><td>4</td><td>8</td><td>5.446</td><td>40.114</td><td>6.797</td><td>0</td><td>0</td></tr>
<tr><td style="text-align:left">Pctl(75)</td><td>2,013</td><td>9.2</td><td>23</td><td>20.007</td><td>71.647</td><td>11.508</td><td>0.03</td><td>0</td></tr>
<tr><td style="text-align:left">Max</td><td>2,013</td><td>12</td><td>31</td><td>97.771</td><td>91.085</td><td>56.388</td><td>4</td><td>1</td></tr>
<tr><td colspan="9" style="border-bottom: 1px solid black"></td></tr></table>

\
[Back to table of contents](#TOC)

## The default regression table


```r
stargazer(output, output2, type = "html")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Change the style

`stargazer` includes several pre-formatted styles that imitate popular
academic journals. Use the `style` argument.


```r
stargazer(output, output2, type = "html", style = "qje")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left"><em>N</em></td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Notes:</em></td><td colspan="2" style="text-align:right"><sup>***</sup>Significant at the 1 percent level.</td></tr>
<tr><td style="text-align:left"></td><td colspan="2" style="text-align:right"><sup>**</sup>Significant at the 5 percent level.</td></tr>
<tr><td style="text-align:left"></td><td colspan="2" style="text-align:right"><sup>*</sup>Significant at the 10 percent level.</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Labelling the table

### Add a title; change the variable labels 


```r
stargazer(output, output2, type = "html", 
          title            = "These are awesome results!",
          covariate.labels = c("Temperature", "Wind speed", "Rain (inches)",
                               "2nd quarter", "3rd quarter", "Fourth quarter"),
          dep.var.caption  = "A better caption",
          dep.var.labels   = "Flight delay (in minutes)")
```


<table style="text-align:center"><caption><strong>These are awesome results!</strong></caption>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2">A better caption</td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">Flight delay (in minutes)</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Temperature</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Wind speed</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Rain (inches)</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">2nd quarter</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">3rd quarter</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Fourth quarter</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Exclude the dependent variable label or the model numbers

Note the dependent variable caption stays. To additionally remove the caption 
add `dep.var.caption = ""`.


```r
stargazer(output, output2, type = "html", 
          dep.var.labels.include = FALSE,
          model.numbers          = FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Change the column names

To change the column names just supply a character vector with the new labels, 
as shown below. 


```r
stargazer(output, output2, type = "html", column.labels = c("Good", "Better"))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>Good</td><td>Better</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Apply a label to more than one column

The option `column.separate` allows for assigning a label to more than one 
column. In this example I told `stargazer` to report each regression twice, for a
total of four columns. Using `column.separate`, `stargazer` now applies the first
label to the first two columns and the second label to the next two columns.


```r
stargazer(output, output, output2, output2, type = "html", 
          column.labels   = c("Good", "Better"),
          column.separate = c(2, 2))
```


<table style="text-align:center"><tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="4"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="4" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="4">delay</td></tr>
<tr><td style="text-align:left"></td><td colspan="2">Good</td><td colspan="2">Better</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td><td>(3)</td><td>(4)</td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.039)</td><td>(0.066)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.285<sup>*</sup></td><td>0.237</td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.013)</td><td>(2.004)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td></td><td>-2.060</td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td>(2.579)</td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td></td><td>-6.755<sup>**</sup></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td>(3.121)</td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td></td><td>-4.408<sup>**</sup></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td>(2.037)</td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>6.152<sup>**</sup></td><td>5.246</td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(2.989)</td><td>(3.412)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.154</td><td>0.175</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.147</td><td>0.161</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="4" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Model names

When the results from different types of regression models (e.g., "OLS", 
"probit") are displayed in the same table `stargazer` adds a row indicating 
model type. Remove these labels by including `model.names = FALSE` (not shown).


```r
stargazer(output, output2, output3, type = "html")
```


<table style="text-align:center"><tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="3"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="3" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="3">delay</td></tr>
<tr><td style="text-align:left"></td><td colspan="2"><em>OLS</em></td><td><em>instrumental</em></td></tr>
<tr><td style="text-align:left"></td><td colspan="2"><em></em></td><td><em>variable</em></td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td><td>(3)</td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td><td>0.041</td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td><td>(0.135)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td><td>0.250</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td><td>(0.198)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td><td>15.335<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td><td>(2.034)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td><td>8.577</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td><td>(8.797)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td><td>0.152</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td><td>0.145</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td><td>12.839 (df = 360)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td><td></td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="3" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Model names (again)

The example above shows the default behavior of `stargazer` is to display only 
one model name (and dependent variable caption) for adjacent columns with the 
same model type. To repeat these labels for all of the columns, do the 
following:


```r
stargazer(output, output2, output3, type = "html",
          multicolumn = FALSE)
```


<table style="text-align:center"><tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="3"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="3" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td>delay</td><td>delay</td><td>delay</td></tr>
<tr><td style="text-align:left"></td><td><em>OLS</em></td><td><em>OLS</em></td><td><em>instrumental</em></td></tr>
<tr><td style="text-align:left"></td><td><em></em></td><td><em></em></td><td><em>variable</em></td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td><td>(3)</td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td><td>0.041</td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td><td>(0.135)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td><td>0.250</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td><td>(0.198)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td><td>15.335<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td><td>(2.034)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td><td>8.577</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td><td>(8.797)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td><td>0.152</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td><td>0.145</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td><td>12.839 (df = 360)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td><td></td></tr>
<tr><td colspan="4" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="3" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Add a custom row to the reported statistics

I use this example to show how to add a row(s), such as reporting fixed effects.


```r
stargazer(output, output2, type = "html",
          add.lines = list(c("Fixed effects?", "No", "No"),
                           c("Results believable?", "Maybe", "Try again later")))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Fixed effects?</td><td>No</td><td>No</td></tr>
<tr><td style="text-align:left">Results believable?</td><td>Maybe</td><td>Try again later</td></tr>
<tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Include R object names


```r
stargazer(output, output2, type = "html", 
          object.names = TRUE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td style="text-align:left"></td><td>output</td><td>output2</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Change the default output

### Report t-statistics or p-values instead of standard errors

Standard errors are reported by default. To report the t-statistics or
p-values instead, see the `report` argument. Notice that I've used this option 
to move the star characters to the t-statistics instead of being next to the
coefficients.


```r
stargazer(output, output2, type = "html",
          report = "vct*")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079</td><td>0.164</td></tr>
<tr><td style="text-align:left"></td><td>t = 2.010<sup>**</sup></td><td>t = 2.479<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285</td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>t = 1.812<sup>*</sup></td><td>t = 1.507</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255</td><td>14.845</td></tr>
<tr><td style="text-align:left"></td><td>t = 7.579<sup>***</sup></td><td>t = 7.407<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>t = -0.799</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755</td></tr>
<tr><td style="text-align:left"></td><td></td><td>t = -2.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408</td></tr>
<tr><td style="text-align:left"></td><td></td><td>t = -2.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152</td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>t = 2.058<sup>**</sup></td><td>t = 1.537</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Report confidence intervals


```r
stargazer(output, output2, type = "html",
          ci = TRUE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.002, 0.156)</td><td>(0.034, 0.294)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(-0.023, 0.594)</td><td>(-0.071, 0.545)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(11.310, 19.200)</td><td>(10.917, 18.773)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-7.116, 2.995)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-12.872, -0.638)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-8.401, -0.416)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(0.293, 12.010)</td><td>(-1.442, 11.934)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Adjust the confidence intervals
 
By default `ci.level = 0.95`. You may also change the character that separates 
the intervals with the `ci.separator` argument.


```r
stargazer(output, output2, type = "html",
          ci = TRUE, ci.level = 0.90, ci.separator = " @@ ")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.014 @@ 0.144)</td><td>(0.055 @@ 0.273)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.026 @@ 0.544)</td><td>(-0.022 @@ 0.496)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(11.944 @@ 18.566)</td><td>(11.548 @@ 18.141)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-6.303 @@ 2.182)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-11.889 @@ -1.622)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-7.759 @@ -1.058)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(1.235 @@ 11.068)</td><td>(-0.367 @@ 10.858)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Robust standard errors (replicating Stata's robust option)

If you want to use robust standard errors (or clustered), `stargazer` allows 
for replacing the default output by supplying a new vector of values to the 
option `se`. For this example I will display the same model twice and adjust the 
standard errors in the second column with the `HC1` correction from the `sandwich`
package (i.e. the same correction Stata uses). 

I also need to adjust the F statistic with the corrected variance-covariance 
matrix (matching Stata's results). Currently, this must be done manually 
(via `add.lines`) as `stargazer` does not (yet) have an option for directly 
replacing the F statistic.

Similar options exist to supply adjusted values to the coefficients, 
t-statistics, confidence intervals, and p-values. See `coef`, `t`, `ci.custom`, 
or `p`.


```r
library(sandwich)
library(lmtest)   # waldtest; see also coeftest.
# Adjust standard errors
cov1         <- vcovHC(output, type = "HC1")
robust_se    <- sqrt(diag(cov1))
# Adjust F statistic 
wald_results <- waldtest(output, vcov = cov1)
stargazer(output, output, type = "html",
          se        = list(NULL, robust_se),
          omit.stat = "f",
          add.lines = list(c("F Statistic (df = 3; 360)", "12.879***", "7.73***")))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.079<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.041)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.285</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.177)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>15.255<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(4.827)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>6.152<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.074)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">F Statistic (df = 3; 360)</td><td>12.879***</td><td>7.73***</td></tr>
<tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.154</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.147</td></tr>
<tr><td style="text-align:left">Residual Std. Error (df = 360)</td><td>12.822</td><td>12.822</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Move the intercept term to the top of the table


```r
stargazer(output, output2, type = "html",
          intercept.bottom = FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Compress the table output

We can condense the table output by placing all the the output on the same row.
When `single.row` is set to `TRUE`, the argument `no.space` is automatically 
set to `TRUE`.


```r
stargazer(output, output2, type = "html",
          single.row = TRUE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup> (0.039)</td><td>0.164<sup>**</sup> (0.066)</td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup> (0.157)</td><td>0.237 (0.157)</td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup> (2.013)</td><td>14.845<sup>***</sup> (2.004)</td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060 (2.579)</td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup> (3.121)</td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup> (2.037)</td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup> (2.989)</td><td>5.246 (3.412)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Omit parts of the default output

In fixed effect model specifications it is often undesirable to report the 
fixed effect coefficients. To omit any coefficient, supply a regular 
expression to `omit`.


```r
stargazer(output, output2, type = "html", omit = "quarter")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Reporting omitted variables

Add the `omit.labels` parameter to report which variables have been omitted. 
Must be the same length as the number of regular expressions supplied to `omit`.
By default `omit.labels` reports "Yes" or "No". To change this supply 
a new vector of length 2 to `omit.yes.no = c("Yes", "No")`.


```r
stargazer(output, output2, type = "html", 
          omit        = "quarter",
          omit.labels = "Quarter dummies?")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Quarter dummies?</td><td>No</td><td>No</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Omit summary statistics


```r
## Remove r-square and f-statistic
stargazer(output, output2, type = "html", 
          omit.stat = c("rsq", "f"))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
See also `keep.stat` a related argument with the opposite behavior.

### Omit whole parts of the table

If you just want to remove parts of the table it is easier to use 
`omit.table.layout` to explicitly specify table elements. See 
`table layout chracters` for a list of codes.


```r
# Remove statistics and notes sections completely
stargazer(output, output2, type = "html", 
          omit.table.layout = "sn")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr></table>

### Omit whole parts of the table (a second way)

Another way to achieve the result above is through the argument `table.layout`. It 
also accepts a character string that tells `stargazer` which table elements to 
**include**.


```r
# Include everything except the statistics and notes sections
stargazer(output, output2, type = "html", 
          table.layout = "-ld#-t-")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr></table>

### Omit the degrees of freedom


```r
stargazer(output, output2, type = "html", 
          df = FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822</td><td>12.718</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup></td><td>12.593<sup>***</sup></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Statistical significance options

By default `stargazer` uses \*\*\*, \*\*, and \* to denote statistical 
significance at the one, five, and ten percent levels. This behavior can be 
changed by altering the `star.char` option. 


```r
stargazer(output, output2, type = "html", 
          star.char = c("@", "@@", "@@@"))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>@@</sup></td><td>0.164<sup>@@</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>@</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>@@@</sup></td><td>14.845<sup>@@@</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>@@</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>@@</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>@@</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>@@@</sup> (df = 3; 360)</td><td>12.593<sup>@@@</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Change the cutoffs for significance

Notice that temperature, quarter3, and quarter4 have each lost a gold star 
because we made it tougher to earn them.


```r
stargazer(output, output2, type = "html", 
          star.cutoffs = c(0.05, 0.01, 0.001)) # star.cutoffs = NULL by default
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>*</sup></td><td>0.164<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285</td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>*</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.05; <sup>**</sup>p<0.01; <sup>***</sup>p<0.001</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Modifying table notes

### Make an addition to the existing note section


```r
stargazer(output, output2, type = "html", 
          notes = "I make this look good!")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
<tr><td style="text-align:left"></td><td colspan="2" style="text-align:right">I make this look good!</td></tr>
</table>

### Replace the note section


```r
stargazer(output, output2, type = "html", 
          notes        = "Sometimes you just have to start over.", 
          notes.append = FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right">Sometimes you just have to start over.</td></tr>
</table>

### Change note alignment


```r
stargazer(output, output2, type = "html", 
          notes.align = "l")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:left"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Change the note section label


```r
stargazer(output, output2, type = "html", 
          notes.label = "New note label")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">New note label</td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Table aesthetics

### Use html tags to modify table elements

```r
# For LaTeX output you can also wrap table text in commands like \textbf{...}, 
# just remember to escape the first backslash, e.g., "A \\textbf{better} caption"
stargazer(output, output2, type = "html", 
          title = "These are <em> awesome </em> results!",  # Italics
          dep.var.caption  = "A <b> better </b> caption")   # Bold
```


<table style="text-align:center"><caption><strong>These are <em> awesome </em> results!</strong></caption>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2">A <b> better </b> caption</td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Change decimal character


```r
stargazer(output, output2, type = "html", 
          decimal.mark = ",")
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0,079<sup>**</sup></td><td>0,164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0,039)</td><td>(0,066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0,285<sup>*</sup></td><td>0,237</td></tr>
<tr><td style="text-align:left"></td><td>(0,157)</td><td>(0,157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15,255<sup>***</sup></td><td>14,845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2,013)</td><td>(2,004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2,060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2,579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6,755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3,121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4,408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2,037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6,152<sup>**</sup></td><td>5,246</td></tr>
<tr><td style="text-align:left"></td><td>(2,989)</td><td>(3,412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0,154</td><td>0,175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0,147</td><td>0,161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12,822 (df = 360)</td><td>12,718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21,860<sup>***</sup> (df = 3; 360)</td><td>12,593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0,1; <sup>**</sup>p<0,05; <sup>***</sup>p<0,01</td></tr>
</table>

### Control the number of decimal places


```r
stargazer(output, output2, type = "html", 
          digits = 1)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.1<sup>**</sup></td><td>0.2<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.04)</td><td>(0.1)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.3<sup>*</sup></td><td>0.2</td></tr>
<tr><td style="text-align:left"></td><td>(0.2)</td><td>(0.2)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.3<sup>***</sup></td><td>14.8<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.0)</td><td>(2.0)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.1</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.6)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.8<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.1)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.4<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.0)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.2<sup>**</sup></td><td>5.2</td></tr>
<tr><td style="text-align:left"></td><td>(3.0)</td><td>(3.4)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.2</td><td>0.2</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.1</td><td>0.2</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.8 (df = 360)</td><td>12.7 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.9<sup>***</sup> (df = 3; 360)</td><td>12.6<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Additional decimal controls

You may also specify the number of additional decimal places to be used if a 
number, when rounded to `digits` decimal places, is equal to zero (Use argument 
`digits.extra`).

My example models do not have any numbers in the thousands, so I won't show 
them, but `digit.separate` and `digits.separator` are also available for 
customizing the output of those characters.


```r
stargazer(output, output2, type = "html", 
          digits       = 1,
          digits.extra = 1)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>0.1<sup>**</sup></td><td>0.2<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.04)</td><td>(0.1)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.3<sup>*</sup></td><td>0.2</td></tr>
<tr><td style="text-align:left"></td><td>(0.2)</td><td>(0.2)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.3<sup>***</sup></td><td>14.8<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.0)</td><td>(2.0)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.1</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.6)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.8<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.1)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.4<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.0)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.2<sup>**</sup></td><td>5.2</td></tr>
<tr><td style="text-align:left"></td><td>(3.0)</td><td>(3.4)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.2</td><td>0.2</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.1</td><td>0.2</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.8 (df = 360)</td><td>12.7 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.9<sup>***</sup> (df = 3; 360)</td><td>12.6<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Drop leading zeros from decimals


```r
stargazer(output, output2, type = "html", 
          initial.zero = FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">temp</td><td>.079<sup>**</sup></td><td>.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(.039)</td><td>(.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>.285<sup>*</sup></td><td>.237</td></tr>
<tr><td style="text-align:left"></td><td>(.157)</td><td>(.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>.154</td><td>.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>.147</td><td>.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Change the order of the variables

The `order` argument will also accept a vector of regular expressions.


```r
stargazer(output, output2, type = "html", 
          order = c(4, 5, 6, 3, 2, 1))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">quarter2</td><td></td><td>-2.060</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.579)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter3</td><td></td><td>-6.755<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(3.121)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">quarter4</td><td></td><td>-4.408<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(2.037)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">wind</td><td>0.285<sup>*</sup></td><td>0.237</td></tr>
<tr><td style="text-align:left"></td><td>(0.157)</td><td>(0.157)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">temp</td><td>0.079<sup>**</sup></td><td>0.164<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.066)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>6.152<sup>**</sup></td><td>5.246</td></tr>
<tr><td style="text-align:left"></td><td>(2.989)</td><td>(3.412)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

### Select which variables to keep in the table

By default `keep = NULL` meaning all variables are included. `keep` accepts a
vector of regular expressions.


```r
# Regex for keep "precip" but not "precipitation"
stargazer(output, output2, type = "html", 
          keep = c("\\bprecip\\b"))
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">delay</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">precip</td><td>15.255<sup>***</sup></td><td>14.845<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(2.013)</td><td>(2.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>364</td><td>364</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.154</td><td>0.175</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.147</td><td>0.161</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>12.822 (df = 360)</td><td>12.718 (df = 357)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>21.860<sup>***</sup> (df = 3; 360)</td><td>12.593<sup>***</sup> (df = 6; 357)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

\
[Back to table of contents](#TOC)

## Feedback

The .Rmd file for this cheatsheet [is on GitHub](https://github.com/JakeRuss/cheatsheets) 
and I welcome suggestions or pull requests.

\
[Back to table of contents](#TOC)
