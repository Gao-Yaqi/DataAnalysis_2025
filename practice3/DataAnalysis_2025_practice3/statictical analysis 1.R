#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
example_df<-read.csv("distribution.csv", header = TRUE,dec = ',', sep = ";")
factor_df <- read.csv("factor_data.csv")
imputed_df <- read.csv("imputed_data.csv")
# Display structure with variable types
str(example_df)
str(factor_df)
str(imputed_df)
#---------------merge two files-------------------------
data_for_analysis <- merge(
  factor_df, 
  imputed_df, 
  by = "record_id",        # column for merge
  all = FALSE       # FALSE = INNER JOIN (only coincidences), TRUE = FULL JOIN
)
str(data_for_analysis)

# save data_for_analysis in CSV
write.csv(data_for_analysis, "data_for_analysis.csv", row.names = FALSE)  


#----------------Check for missing values----------------
# Check the number of missing values
sum(is.na(data_for_analysis$lipids5)) 

# View the distribution of missing values (by group)
table(data_for_analysis$outcome, is.na(data_for_analysis$lipids5))

#----------------Handling missing values----------------
library(dplyr)
data_for_analysis <- data_for_analysis %>%
  group_by(outcome) %>%
  mutate(lipids5 = ifelse(is.na(lipids5), 
                          median(lipids5, na.rm = TRUE),
                          lipids5)) %>%
  ungroup()


#------------------Probability Distributions----------------------- 
install.packages("MASS", dependencies=T)
library(MASS)
#----------------example-----------------------------------------
summary (example_df)
example_df$value <- as.numeric(example_df$value)
summary(example_df)
#building histograms for example
# normal distribution
val<-example_df[example_df$distribution=="norm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="normal")

fit
#lognormal distribution
val<-example_df[example_df$distribution=="lognorm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="lognormal")

fit

unname(fit$estimate[1])

unname(fit$estimate[2])

m_log<-exp(unname(fit$estimate[1]))*sqrt(exp(unname(fit$estimate[2])^2))
m_log
sd_log<-sqrt(exp(2*unname(fit$estimate[1]))*(exp(unname(fit$estimate[2])^2)-1)*sqrt(exp(unname(fit$estimate[2])^2)))
sd_log
#exponential distribution

val<-example_df[example_df$distribution=="exp",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="exponential")

fit

unname(fit$estimate[1])

m_exp<-1/unname(fit$estimate[1])
m_exp
#Poisson distribution
val<-example_df[example_df$distribution=="pois",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="Poisson")

fit

unname(fit$estimate[1])

sd_pois<-sqrt(unname(fit$estimate[1]))
sd_pois
#Selecting a Distribution Model

val<-example_df[example_df$distribution=="lognorm",]$value

fit_1<-fitdistr(val, densfun="normal")
fit_2<-fitdistr(val, densfun="lognormal")
fit_3<-fitdistr(val, densfun="exponential")

#Bayesian Information Criterion calculation
BIC(fit_3)

#calculation of the Bayesian information criterion for all models
BIC_value<-c(BIC(fit_1), BIC(fit_2), BIC(fit_3))

#forming a vector with the name of the models
distribution<-c("normal", "lognormal", "exponential")

#combining the results into a final table
rez<-data.frame(BIC_value=BIC_value, distribution=distribution)

#sort table in ascending order of Bayesian Information Criterion value
rez<-rez[order(rez$BIC_value, decreasing=F),]

rez


#calculation of absolute values of the confidence interval for the mean of a lognormal distribution
error_min<-unname(fit_2$estimate[1])-unname(fit_2$sd[1])
error_max<-unname(fit_2$estimate[1])+unname(fit_2$sd[1])

error_min
error_max

m<-exp(unname(fit_2$estimate[1]))*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_min<-exp(error_min)*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_max<-exp(error_max)*sqrt(exp(unname(fit_2$estimate[2])^2))

value_error_min
m
value_error_max

#--------------data for analysis--------------------------
#building histograms
value_d1<-data_for_analysis$lipids1
hist(value_d1)
value_d2<-data_for_analysis$lipids2
hist(value_d2)
value_d3<-data_for_analysis$lipids3
hist(value_d3)
value_d4<-data_for_analysis$lipids4
hist(value_d4)


# d1 distribution estimate


fit_d1_1<-fitdistr(value_d1,densfun="normal")
fit_d1_2<-fitdistr(value_d1,densfun="lognormal")
fit_d1_3<-fitdistr(value_d1,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d1

BIC_value_d1 <- c(BIC(fit_d1_1),BIC(fit_d1_2),BIC(fit_d1_3))
distribution <-c("normal","lognormal","exponential")
result_d1<-data.frame(BIC_value_d1=BIC_value_d1, distribution=distribution)
result_d1
min(result_d1$BIC_value_d1)
distribution_d1<-result_d1[result_d1$BIC_value_d1==min(result_d1$BIC_value_d1),]$distribution
distribution_d1
# Finding parameters for d1
fit_d1_1$estimate[1:2]

# d2 distribution estimate


fit_d2_1<-fitdistr(value_d2,densfun="normal")
fit_d2_2<-fitdistr(value_d2,densfun="lognormal")
fit_d2_3<-fitdistr(value_d2,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d2

BIC_value_d2 <- c(BIC(fit_d2_1),BIC(fit_d2_2),BIC(fit_d2_3))
distribution <-c("normal","lognormal","exponential")
result_d2<-data.frame(BIC_value_d2=BIC_value_d2, distribution=distribution)
result_d2
min(result_d2$BIC_value_d2)
distribution_d2<-result_d2[result_d2$BIC_value_d2==min(result_d2$BIC_value_d2),]$distribution
distribution_d2
# Finding parameters for d2
fit_d2_1$estimate[1:2]

#------------------ d3 distribution estimate --------------------
fit_d3_1 <- fitdistr(value_d3, densfun = "normal")
fit_d3_2 <- fitdistr(value_d3, densfun = "lognormal")
fit_d3_3 <- fitdistr(value_d3, densfun = "exponential")

# BIC计算与最佳分布选择
BIC_value_d3 <- c(BIC(fit_d3_1), BIC(fit_d3_2), BIC(fit_d3_3))
distribution <- c("normal", "lognormal", "exponential")
result_d3 <- data.frame(BIC_value_d3 = BIC_value_d3, distribution = distribution)

# 显示结果
result_d3
distribution_d3 <- result_d3[result_d3$BIC_value_d3 == min(result_d3$BIC_value_d3), ]$distribution
distribution_d3

# 最佳分布参数估计
if(distribution_d3 == "normal") fit_d3_1$estimate[1:2]
if(distribution_d3 == "lognormal") fit_d3_2$estimate[1:2]
if(distribution_d3 == "exponential") fit_d3_3$estimate[1]

#------------------ d4 distribution estimate --------------------
fit_d4_1 <- fitdistr(value_d4, densfun = "normal")
fit_d4_2 <- fitdistr(value_d4, densfun = "lognormal")
fit_d4_3 <- fitdistr(value_d4, densfun = "exponential")

# BIC计算与最佳分布选择
BIC_value_d4 <- c(BIC(fit_d4_1), BIC(fit_d4_2), BIC(fit_d4_3))
distribution <- c("normal", "lognormal", "exponential")
result_d4 <- data.frame(BIC_value_d4 = BIC_value_d4, distribution = distribution)

# 显示结果
result_d4
distribution_d4 <- result_d4[result_d4$BIC_value_d4 == min(result_d4$BIC_value_d4), ]$distribution
distribution_d4

# 最佳分布参数估计
if(distribution_d4 == "normal") fit_d4_1$estimate[1:2]
if(distribution_d4 == "lognormal") fit_d4_2$estimate[1:2]
if(distribution_d4 == "exponential") fit_d4_3$estimate[1]
#-----------descriptive statistics------------------
# -----------------按分布类型生成描述统计表----------------
# 自定义统计函数
custom_stats <- function(variable, distribution) {
  switch(distribution,
         "normal" = list("{mean} ({sd})"),
         "lognormal" = list("exp({mean})*sqrt(exp({sd}^2)) ({exp(mean)*sqrt(exp(sd^2)*(sqrt(exp(sd^2)-1)})"),
         "exponential" = list("{1/rate} ({1/rate})")
  )
}

dist_types <- c(
  lipids1 = "normal",
  lipids2 = "lognormal",
  lipids3 = "exponential",
  lipids4 = "normal"
)

# 创建自定义统计表格
custom_table <- tbl_summary(
  data = data_for_analysis,
  by = outcome,
  include = c(lipids1, lipids2, lipids3, lipids4),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})"  # 默认显示均值和标准差
  )
)

# 添加Brunner-Munzel检验p值
# -----------------添加统计检验结果----------------
add_brunner_p <- function(data, variable) {
  group0 <- data[[variable]][data$outcome == "0"]
  group1 <- data[[variable]][data$outcome == "1"]
  test <- brunner.munzel.test(group0, group1)
  return(test$p.value)
}

p_values <- sapply(c("lipids1", "lipids2", "lipids3", "lipids4"), 
                   function(x) add_brunner_p(data_for_analysis, x))

# 合并到表格
final_table <- custom_table %>%
  add_p() %>%
  modify_table_body(
    ~ .x %>%
      dplyr::mutate(
        p.value = ifelse(variable %in% names(p_values),
                         p_values[variable],
                         p.value)
      )
  )

# 显示最终表格
final_table

#-----------for publication tables-----------------
install.packages("gtsummary")
library(gtsummary)

tbl_summary(data_for_analysis)  # Automatic table
tbl_summary(data_for_analysis, by = outcome)  # By groups

#---------------Creating a custom table--------------
# Homework: Creating a custom table with descriptive statistics results
# 生成分组统计表
tbl_summary(
  data = data_for_analysis,
  by = outcome,  # 按 outcome 分组
  include = c(lipids1, lipids2, lipids3, lipids4),  # 指定分析的变量
  type = list(all_continuous() ~ "continuous"),  # 确保所有变量作为连续型处理
  statistic = all_continuous() ~ "{mean} ({sd})"  # 显示均值±标准差
) %>%
  add_p()  # 添加组间差异检验的 p 值

library(ggplot2)
library(tidyr)

# 将数据转换为长格式（便于分面）
data_long <- data_for_analysis %>%
  pivot_longer(cols = c(lipids1, lipids2, lipids3, lipids4), 
               names_to = "variable", 
               values_to = "value")

# 绘制分面直方图
ggplot(data_long, aes(x = value)) +
  geom_histogram(fill = "lightblue", bins = 20) +
  facet_grid(variable ~ outcome, scales = "free") +  # 按变量和 outcome 分面
  labs(title = "Distribution of Lipids by Outcome Group")

#--------------Statistical Tests---------------------
value_outcome1<-data_for_analysis[data_for_analysis$outcome=="1",]$lipids1
hist(value_outcome1)
value_outcome0<-data_for_analysis[data_for_analysis$outcome=="0",]$lipids1
hist(value_outcome0)

#-------Levene's Test for Homogeneity of Variance--------------
install.packages("car")
library(car)
str(data_for_analysis)
data_for_analysis$outcome<- as.factor(data_for_analysis$outcome)
car::leveneTest(lipids1 ~ outcome, data = data_for_analysis)
#---------------Application of the Brunner-Munzel test----------
install.packages("lawstat")
library(lawstat)
group1 <- data_for_analysis$lipids1[data_for_analysis$outcome == "0"]
group2 <- data_for_analysis$lipids1[data_for_analysis$outcome == "1"]

brunner.munzel.test(group1, group2)
#-------------comparison of results with other tests--------------
t.test(group1, group2)
wilcox.test(group1, group2)

#----------------------------EDA----------------------------------
install.packages("DataExplorer")
library(DataExplorer)
create_report(data_for_analysis)  # Generates HTML report with graphs and statistics
create_report(
  data = data_for_analysis,
  output_file = "EDA_Report.html",  
  output_dir = getwd(),                
  report_title = "EDA Report"          
)