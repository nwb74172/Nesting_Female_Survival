library(tidyverse)

df <- read_csv("D:/Nest_Structure_Female_Mortality/Died_On_Nest_Veg/wla_nuds_roble_readings.csv")

df$overall_avg_nuds <- rowMeans(df[, 15:18], na.rm = TRUE)

df$overall_avg_robel <- rowMeans(df[, c(4,7,10,13)], na.rm = TRUE)

df$overall_min_robel <- rowMeans(df[, c(5,8,11,14)], na.rm = TRUE)

df$overall_max_robel <- rowMeans(df[, c(3,6,9,12)], na.rm = TRUE)

df$overall_avg_nuds_height <- rowMeans(df[, c(19:22)], na.rm = TRUE)

df$new_column <- ifelse(df[[2]] == "n", 1,
                        ifelse(df[[2]] == "r", 0, NA))


cor(df[, c(23:27)], method = "pearson")
