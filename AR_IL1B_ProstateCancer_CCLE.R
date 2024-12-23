CCLE AR activity # Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
x <- read.csv("C:/Users/u/Downloads/CCLE_prostateCancer.csv", header = TRUE)
names(x)[7] <- 'Type'
names(x)[8] <- 'IL1B'
names(x)[9] <- 'AR'
names(x)[10] <- 'CX3CR1'

# View the structure of the data (optional)
# str(data)

# Assuming gene expression columns start from the 8th column
gene_columns <- colnames(x)[8:ncol(x)]

# Convert data from wide to long format
long_data <- x %>%
  pivot_longer(cols = all_of(gene_columns), 
               names_to = "Gene", 
               values_to = "Expression") %>%
  filter(!is.na(Expression))  # Remove rows with NA values in Expression

# Calculate average AR expression for each displayName
ar_order <- long_data %>%
  filter(Gene == "AR") %>%
  group_by(displayName) %>%
  summarise(AR_Expression = mean(Expression, na.rm = TRUE)) %>%
  arrange(desc(AR_Expression))

long_data$displayName <- factor(long_data$displayName, 
                                levels = ar_order$displayName)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
x <- read.csv("C:/Users/u/Downloads/CCLE_prostateCancer.csv", header = TRUE)
names(x)[7] <- 'Type'
names(x)[8] <- 'IL1B'
names(x)[9] <- 'AR'
names(x)[10] <- 'CX3CR1'

# Convert data from wide to long format
gene_columns <- colnames(x)[8:ncol(x)]

long_data <- x %>%
  pivot_longer(cols = all_of(gene_columns), 
               names_to = "Gene", 
               values_to = "Expression") %>%
  filter(!is.na(Expression))  # Remove rows with NA values in Expression

# Calculate average AR expression for each displayName
ar_order <- long_data %>%
  filter(Gene == "AR") %>%
  group_by(displayName) %>%
  summarise(AR_Expression = mean(Expression, na.rm = TRUE)) %>%
  arrange(desc(AR_Expression))

# Set displayName as a factor ordered by average AR expression
long_data$displayName <- factor(long_data$displayName, 
                                levels = ar_order$displayName)

# Create the bar plot with adjusted spacing and white background
p <- ggplot(long_data, aes(x = displayName, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Expression per Tumor Type", 
       x = "Tumor Type", 
       y = "Expression Level") +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +  # Adjust gap
  theme_light() +  # Use theme_light() for a white background
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

correlation_data <- long_data %>%
  filter(Gene %in% c("IL1B", "AR")) %>%
  pivot_wider(names_from = Gene, values_from = Expression) %>%
  select(displayName, IL1B, AR) %>%
  drop_na()  # Remove rows with NA values

# Calculate correlation coefficient
correlation_result <- cor(correlation_data$IL1B, correlation_data$AR)


# Save the plot as a PNG file
ggsave("C:/Users/u/Downloads/gene_expression_plot.png", plot = p, width = 10, height = 6)
scatter_plot <- ggplot(correlation_data, aes(x = IL1B, y = AR)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of IL1B vs AR Expression",
       x = "IL1B Expression",
       y = "AR Expression") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add regression line
  annotate("text", x = max(correlation_data$IL1B, na.rm = TRUE) * 0.8, 
           y = max(correlation_data$AR, na.rm = TRUE) * 0.8,
           label = paste("Correlation: ", round(correlation_result, 2)),
           size = 5, color = "black") +
  theme_light()

# Save the scatter plot as a PNG file
ggsave("C:/Users/u/Downloads/il1b_ar_scatter_plot.png", plot = scatter_plot, width = 10, height = 6)

# Print the correlation result to console
print(paste("Correlation between IL1B and AR: ", round(correlation_result, 2)))

