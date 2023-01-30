library(fs)
library(ggplot2)
library(dplyr)

# Get file paths for all .csv files in the directory
file_paths <- fs::dir_ls("/Users/david/Desktop/C++/Exciton Research/Heat Map R/2023-01-28 21:58:41 Data")

# Iterate through each file
for (i in seq_along(file_paths)) {
  # Read the file into a data frame
  excitons_data <- read.csv(file_paths[[i]])
  
  # Create a new dataframe that represents the count of excitons at each position on the grid
  excitons_grid <- data.frame(X = rep(1:100, each = 100), Y = rep(1:100, 100), count = 0)
  
  # Count the number of excitons at each position on the grid
  for (j in 1:nrow(excitons_data)) {
    x <- excitons_data$X_pos[j]
    y <- excitons_data$Y_pos[j]
    excitons_grid[excitons_grid$X == x & excitons_grid$Y == y, "count"] <- excitons_grid[excitons_grid$X == x & excitons_grid$Y == y, "count"] + 1
  }
  
  # Create a heatmap using ggplot
  ggplot(data = excitons_grid, aes(x = X, y = Y, fill = count)) + 
    geom_tile() + 
    scale_fill_gradient(limits=c(0,5),low = "blue", high = "red",breaks=seq(0,5,1)) + 
    ggtitle(paste0("Excitons Heatmap for File ",i)) + 
    xlab("X Position") + 
    ylab("Y Position") + 
    theme_void()
  ggsave(paste0("/Users/david/Desktop/C++/Exciton Research/Heat map images/heatmap_","file_",i,".png"))
}



  








'#
  # Create a heatmap using ggplot
  print(ggplot(data = excitons_grid, aes(x = X, y = Y, fill = count)) + 
          geom_tile() + 
          scale_fill_gradient(low = "blue", high = "red") + 
          ggtitle("Excitons Heatmap") + 
          xlab("X Position") + 
          ylab("Y Position") + 
          theme_void())
  
  # Save the plot to a file
  ggsave(paste0(csv_path, "/heatmap_", basename(file), ".png"))
}

  
  # Create a heatmap using ggplot
  print(ggplot(data = excitons_grid, aes(x = X, y = Y, fill = count)) + 
          geom_tile() + 
          scale_fill_gradient(low = "blue", high = "red") + 
          ggtitle("Excitons Heatmap") + 
          xlab("X Position") + 
          ylab("Y Position") + 
          theme_void())
  
  # Save the plot to a file
  ggsave(paste0(csv_path, "/heatmap_", basename(file), ".png"))
}
#'





