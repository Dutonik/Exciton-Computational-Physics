library(fs)
library(ggplot2)
library(dplyr)

#VARIANCE CALCULATIONS
# Read the variance data from the csv file
variance_data <- read.csv("/Users/david/Desktop/C++/Exciton Research/Heat Map R/2023-04-11 13:32:04 Variance/Variance.csv")

# Load the rsq package
library(rsq)

# Plot the data
attach(variance_data)
plot(Time, Variance, main="Variance vs Time", xlab="Time (cycles)", ylab = expression(paste("Variance (", "tile"^2, ")")), pch=19, cex = 0.1)
fit <- lm(Variance ~ Time)
abline(fit, col="red")
rsquared <- round(rsq(fit), digits = 3) # calculate R^2 and round to 3 decimal places
legend("topleft", legend=bquote(hat(y) == .(format(coef(fit)[2], digits=3))*x + .(format(coef(fit)[1], digits=3)) ~ "; " ~ R^2 == .(rsquared)),
       col="red", lty=1)

# Extract the variance value at time = 0
variance_value <- variance_data[variance_data$Time == 0, "Variance"]
print(variance_value)
# Plot the Gaussian
x <- seq(-50, 50, length.out = 1000)
y <- exp(-(x^2/variance_value))
plot(x, y, type = "l", xlab = expression(paste("Position (", "Unit Tile Lengths", ")")), ylab = "exp(-x^2/variance)", main = "Gaussian Plot (t=0)")

# Extract the variance value at time = 59 cycles
variance_value2 <- variance_data[variance_data$Time == 58, "Variance"]
print(variance_value2)
# Plot the Gaussian
x <- seq(-50, 50, length.out = 1000)
y <- exp(-(x^2/variance_value2))
plot(x, y, col = 'red', type = "l", xlab = expression(paste("Position (", "Unit Tile Lengths", ")")), ylab = "exp(-x^2/variance)", main = "Gaussian Plot (t=59 cycles)")

# Extract the variance value at time = 590 cycles
variance_value3 <- variance_data[variance_data$Time == 589, "Variance"]
print(variance_value3)
# Plot the Gaussian
x <- seq(-50, 50, length.out = 1000)
y <- exp(-(x^2/variance_value3))
plot(x, y, col = 'blue', type = "l", xlab = expression(paste("Position (", "Unit Tile Lengths", ")")), ylab = "exp(-x^2/variance)", main = "Gaussian Plot (t=590 cycles)")





###############





# Read the CSV file
exciton_data <- read.csv("/Users/david/Desktop/C++/Exciton Research/Heat Map R/2023-04-11 13:12:13 eDensity1/Exciton_Regional_Density.csv")

# Rename the columns for better readability
colnames(exciton_data) <- c("concentration", "time")

# Create the plot
density_vs_time_plot <- ggplot(exciton_data, aes(x = time, y = density)) +
  geom_point(size=0.5) +
  geom_line() +
  labs(title = "Exciton Regional Concentration vs Time",
       x = "Time (cycles)",
       y = "Concentration (excitons/unit tile)") +
  theme_minimal()

# Display the plot
print(density_vs_time_plot)




##############





#POSITION HEAT HAP
# Get file paths for all .csv files in the directory
file_paths <- fs::dir_ls("/Users/david/Desktop/C++/Exciton Research/Heat Map R/2023-04-06 21:20:42 Positons")

# Iterate through each file
for (i in seq_along(file_paths)) {
  # Read the file into a data frame
  excitons_data <- read.csv(file_paths[[i]])
  
  # Create a new dataframe that represents the count of excitons at each position on the grid
  excitons_grid <- data.frame(X = rep(1:100, each = 100), Y = rep(1:100, times = 100), count = 0)  

  # Count the number of excitons at each position on the grid
  for (j in 1:nrow(excitons_data)) {
    x <- excitons_data$X_pos[j]
    y <- excitons_data$Y_pos[j]
    excitons_grid[excitons_grid$X == x & excitons_grid$Y == y, "count"] <- excitons_grid[excitons_grid$X == x & excitons_grid$Y == y, "count"] + 1
  }
  

  # Create a heatmap using ggplot
  ggplot(data = excitons_grid, aes(x = X, y = Y, fill = count)) + 
    geom_tile() + 
    scale_fill_gradient(limits=c(0,20),low = "blue", high = "red",breaks=seq(0,20,2)) + 
    ggtitle(paste0("Excitons Heatmap for File ",i)) + 
    xlab("X Position") + 
    ylab("Y Position") + 
    theme_void()
  #Save the output images to a specified folder with a unique label
  ggsave(paste0("/Users/david/Desktop/C++/Exciton Research/Heat map images/heatmap_","file_",i,".png"))
}

