#--------------------start-------------------------------
# Get current working directory
getwd()
# create path relative to project root
#data_path <- "/cloud/project/data/DataSet_No_Details.csv"
data_path <- "DataSet_No_Details.csv"
#----------------read dataset--------------------------
df <- read.csv(data_path)
# Display structure with variable types
str(df)
# Beautiful summary with histograms for numeric variables
install.packages("skimr")
library(skimr)
skim(df) 
#---------------data set preparation------------------
library(dplyr)
# Delete a few columns 
cols_to_remove <- c("h_index_34", "h_index_56", "hormone10_1", "hormone10_2","an_index_23","outcome","factor_eth","factor_h","factor_pcos","factor_prl")
MD_df <- df %>% select(-any_of(cols_to_remove))
factor_df <- df %>% select (record_id, outcome, factor_eth, factor_h,factor_pcos,factor_prl)
str(MD_df)
summary(factor_df)

#--------------Identify Missing Values-----------------
sum(is.na(MD_df))               # Total NAs in entire dataset
colSums(is.na(MD_df))           # NA counts per column
skim(MD_df)
na_stats <- colMeans(is.na(MD_df)) * 100 # % missing data
na_stats
na_stats_filtered <- na_stats[na_stats <= 35] #  missing data <=35 %
# result as a table
data.frame(
  Column = names(na_stats_filtered),
  NA_Percent = na_stats_filtered,
  row.names = NULL
)

na_stats_filtered_1 <- na_stats[na_stats > 35] # missing data >35 %
# result as a table
data.frame(
  Column = names(na_stats_filtered_1),
  NA_Percent = na_stats_filtered_1,
  row.names = NULL
)

#-------------------Visualizing Missing Data Patterns------------------
install.packages(visdat)
library(visdat)
vis_miss(MD_df)  # Visualizes NA patterns

install.packages(naniar)
library(naniar)
gg_miss_var(MD_df)  # Barplot of missingness per variable
#------------------ Analyzing the Impact of Missing Data--------------
# Delete a few columns 
library(dplyr)
cols_to_remove1 <- c("hormone9", "hormone11", "hormone12", "hormone13","hormone14")
handle_MD_df <- MD_df %>% select(-any_of(cols_to_remove1))
str(handle_MD_df)

#------------------Performing Little's MCAR Test----------------------
#Homework!!!
#Hypotheses:
#  H₀ (Null Hypothesis): Data is MCAR.
library(naniar)

handle_MD_df_numeric <- handle_MD_df %>% select(-record_id)
mcar_test_result <- mcar_test(handle_MD_df_numeric)
print(mcar_test_result)

#H₁ (Alternative Hypothesis): Data is not MCAR (either MAR or MNAR).

if (mcar_test_result$p.value > 0.05) {
  cat("Fail to reject H0: Data is likely MCAR.\n")
} else {
  cat("Reject H0: Data is likely not MCAR (either MAR or MNAR).\n")
}
#If p-value > 0.05, we fail to reject H₀ (data is likely MCAR).
#If p-value ≤ 0.05, we reject H₀ (data is likely not MCAR).
#------------------Imputation with MICE-------------------------------
# Install packages if they are not already installed
install.packages(c("mice", "ggplot2", "naniar"))
# Load the packages
library(mice)
library(ggplot2)
library(naniar)

# Perform Multiple Imputation
imputed_handle_MD_df <- mice(handle_MD_df, m=5, method='pmm', print=FALSE)
# Perform Multiple Imputation 
#Random Forest method 
#------For complex nonlinear relationships between variables------------

imputed_handle_MD_df <- mice(handle_MD_df[, !names(handle_MD_df) %in% "New"], method="rf")  
imputed_handle_MD_df_final <- complete(imputed_handle_MD_df)  # generate full data
# Density plots 
ggplot(handle_MD_df, aes(x=hormone10_generated, fill="Original")) +
  geom_density(alpha=0.5) +
  geom_density(data=imputed_handle_MD_df_final, aes(x=hormone10_generated, fill="Imputed"), alpha=0.5) +
  labs(title="Density Plot of hormone10_generated: Original vs. Imputed")+
  scale_x_continuous(limits = c(0, 2))

#Predictive Mean Matching 
#------default for numeric data------------
#Homework!!!
imputed_handle_MD_df1 <- mice(handle_MD_df[, !names(handle_MD_df) %in% "New"], method="pmm")  
imputed_handle_MD_df_final1 <- complete(imputed_handle_MD_df)  # generate full data
# Density plots 
ggplot(handle_MD_df, aes(x=hormone10_generated, fill="Original")) +
  geom_density(alpha=0.5) +
  geom_density(data=imputed_handle_MD_df_final1, aes(x=hormone10_generated, fill="Imputed"), alpha=0.5) +
  labs(title="Density Plot of hormone10_generated: Original vs. Imputed")+
  scale_x_continuous(limits = c(0, 2))

write.csv(imputed_handle_MD_df_final,
          "imputed_handle_MD_df_final.csv",
          row.names = FALSE)
                     
#----------------Outlier Detection Methods------------------------
library(ggplot2)
library(tidyr)
outliers_data <- imputed_handle_MD_df_final %>%
  select(lipids1, lipids2, lipids3, lipids4, lipids5) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# build a graph 
ggplot(outliers_data, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Outlier Detection",
       x = "variables",
       y = "value") +
  theme_minimal()
#build a graph for all dataset
imputed_handle_MD_df_final %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Boxplots for Outlier Detection")

# 加载必要的包
library(dbscan)
library(ggplot2)

# 假设你的数据框是 imputed_handle_MD_df_final，包含 hormone10_generated 列
# 提取需要进行 LOF 分析的列，假设是 'hormone10_generated'
lof_data <- imputed_handle_MD_df_final[, c("hormone10_generated")]

# 将数据转换为矩阵
lof_data_matrix <- as.matrix(lof_data)

# 进行局部离群因子（LOF）分析，设置 minPts=11（邻域大小，minPts 值可以根据需要调整）
lof_result <- lof(lof_data_matrix, minPts = 11)

# 将 LOF 分数添加到数据框中
imputed_handle_MD_df_final$lof_score <- lof_result

# 创建 LOF 分数的可视化，散点图显示 hormone10_generated 和 LOF 分数的关系
ggplot(imputed_handle_MD_df_final, aes(x = hormone10_generated, y = lof_score)) +
  geom_point() + 
  theme_minimal() +
  labs(title = "Local Outlier Factor (LOF) Analysis", 
       x = "Hormone 10 Generated", 
       y = "LOF Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
