# Setting the working directory to the current directory
setwd(getwd())

# Loading the tidyr library
library(tidyr)

# Reading the data from the DecayTimecourse.txt file
Timecourse_data <- read.delim("DecayTimecourse.txt")

# Converting the data to a data frame
Timecourse_data <- data.frame(Timecourse_data)

# Extracting the relevant columns for each time course
tc1 <- Timecourse_data[,c(1:10)] 
tc2 <- Timecourse_data[,c(1, 11:19)]
tc3 <- Timecourse_data[,c(1, 20:28)]

# Setting the column names for each time course using the first row of data
colnames(tc1) <- tc1[1,]
colnames(tc2) <- tc2[1,]
colnames(tc3) <- tc3[1,]

# Removing the first row from each time course
tc1 <- tc1[-1,]
tc2 <- tc2[-1,]
tc3 <- tc3[-1,]

# Removing rows where all values are NA
tc1 <- tc1[!apply(tc1[,c(2:10)], 1, function(x){all(is.na(x))}),]
tc2 <- tc2[!apply(tc2[,c(2:10)], 1, function(x){all(is.na(x))}),]
tc3 <- tc3[!apply(tc3[,c(2:10)], 1, function(x){all(is.na(x))}),]

# Extracting the first column from each time course
Trans_1 <- data.frame(tc1[,1])
Trans_2 <- data.frame(tc2[,1])
Trans_3 <- data.frame(tc3[,1])

# Applying the natural logarithm to the remaining columns of each time course
tc1 <- log(tc1[,2:10])
tc2 <- log(tc2[,2:10])
tc3 <- log(tc3[,2:10])

# Defining a function to calculate the slope for each row of a time course
reg <- function(a){
  slope <- c()
  for (i in c(1:nrow(a))){
    k <- c()
    for (j in a[i,]) {
      k = c(k,j)
    }
    k <- replace(k, is.infinite(k) & k <0, NA)
    x <- as.integer(colnames(a)[-which(is.na(k) == TRUE)])
    y <- k[-which(is.na(k)==TRUE)]
    if (length(y) != 0){
      z = data.frame(x,y)
      l = lm(z[,2]~z[,1], data = z)
      slope = c(slope, l$coefficients[2])
    }else{slope = c(slope, NA)}
  }
  return(slope)
}

# Calculating the slope for each time course
slope_tc1 <- data.frame(reg(tc1))
slope_tc2 <- data.frame(reg(tc2))
slope_tc3 <- data.frame(reg(tc3))

# Calculating the half life for each time course
half_life1 <- cbind(Trans_1, log(2)/slope_tc1)
half_life2 <- cbind(Trans_2, log(2)/slope_tc2)
half_life3 <- cbind(Trans_3, log(2)/slope_tc3)

# Setting the column names for the half life data frames
colnames(half_life1)[1] <- "YORF"
colnames(half_life2)[1] <- "YORF"
colnames(half_life3)[1] <- "YORF"

# Merging the half life data frames
merged_half_lifes <- merge(half_life1, half_life2, by = "YORF", all = TRUE)
final_merged_half_life <- merge(merged_half_lifes, half_life3, by = "YORF", all = TRUE)

# Calculating the mean half life for each gene
mean <- rowMeans(final_merged_half_life[,c(2:4)])

# Combining the gene names and mean half lives into a single data frame
mean <- cbind(final_merged_half_life[,c(1)], mean)

# Sorting the genes by mean half life
sorted_means <- mean[order(-mean$mean),]

# Removing rows with missing values
sorted_means <- na.omit(sorted_means)

# Extracting the top and bottom genes
top_genes <- head(sorted_means, n = 618)
bottom_genes <- tail(sorted_means, n = 618)

# Writing the top and bottom genes to csv files
write.csv(top_genes, file = "top_genes.csv", row.names = FALSE)
write.csv(bottom_genes, file = "bottom_genes.csv", row.names = FALSE)