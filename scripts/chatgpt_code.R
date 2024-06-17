# Load the necessary library
install.packages("tidyverse")
library(tidyverse)

# Check the current working directory
getwd()

#Set working directory
setwd("~/Users/kitchel/Dropbox")


# Set a new working directory if needed
# setwd("path/to/your/directory")

# Create a data frame
data <- tibble(
  Name = c("Alice", "Bob", "Charlie", "David", "Eve"),
  Age = c(23, 35, 45, 25, 30),
  Gender = c("F", "M", "M", "M", "F"),
  Score = c(85, 90, 78, 88, 95)
)

# Print the data frame
print("Original Data Frame:")
print(data)

# Perform some basic operations
# Filter rows where Age is greater than 30
filtered_data <- data %>%
  filter(Age > 30)

print("Filtered Data Frame (Age > 30):")
print(filtered_data)

# Add a new column with Score squared
mutated_data <- data %>%
  mutate(Score_Squared = Score^2)

print("Data Frame with Score Squared:")
print(mutated_data)

# Select specific columns
selected_data <- data %>%
  select(Name, Score)

print("Data Frame with Selected Columns (Name, Score):")
print(selected_data)

# Save the original data frame to a CSV file
write_csv(data, "data_frame.csv")
print("Data frame saved to 'data_frame.csv'")

# Save the filtered data frame to an RDS file
saveRDS(filtered_data, "filtered_data_frame.rds")
print("Filtered data frame saved to 'filtered_data_frame.rds'")

# Load the filtered data frame from the RDS file
loaded_filtered_data <- readRDS("filtered_data_frame.rds")
print("Loaded Filtered Data Frame from RDS file:")
print(loaded_filtered_data)

#practice test practice test