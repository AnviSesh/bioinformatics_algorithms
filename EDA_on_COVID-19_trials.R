library(tidyverse)
library(lubridate)
library(janitor)
#loading the dataset from kaggle 
df <- read.csv("C:/Users/ANVITHA/Downloads/archive (1)/COVID clinical trials.csv",
               stringsAsFactors = FALSE)
#clean_names helps you standardise column names
df <- clean_names(df)
#Helps you view the data types of each column
str(df)
#summary helps with type of data present in each column, if its chr data type, its gonna depict class,mode and length. If it is a int dtype, MIN,1st Quartile, median, mean, 3rd Quartile and Max are shown.
summary(df)
#To count the number of rows and columns
dim(df)
#Check whether the dataset has missing values or not
colSums(is.na(df)) %>% sort(decreasing = TRUE)
#You're changing the format into an R in-built object
df$start_date <- as.Date(df$start_date, format="%B %d, %Y")
df$completion_date <- as.Date(df$completion_date, format="%B %d, %Y")
df$start_year <- year(df$start_date)

df %>%
  count(status) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x = reorder(status, n), y = n)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Trial Status Distribution")


df %>%
  count(phases) %>%
  ggplot(aes(x = phases, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df %>%
  count(locations) %>%
  top_n(10, n) %>%
  ggplot(aes(x = reorder(locations, n), y = n)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Top Countries Conducting Trials")


df %>%
  count(sponsor_collaborators) %>%
  top_n(10, n)

intervention_data <- df %>%
  count(interventions)
ggplot(intervention_data, aes(x = "", y = n, fill = interventions)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            size = 3) +
  labs(title = "Distribution of Intervention Types") +
  theme_void()
head(interventions)

df %>%
  count(start_year) %>%
  ggplot(aes(x = start_year, y = n)) +
  geom_line() +
  geom_point() +
  labs(title="Trials Over Time")

df$duration_days <- as.numeric(df$completion_date - df$start_date)

summary(df$duration_days)

ggplot(df, aes(x = duration_days)) +
  geom_histogram(bins = 50)
