# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  dir <- args[1]  
} else {
  quit(status = 0)
}

# Define species name and the list of clusters
species <- sub('.*5_(.*?)_clustering.*', '\\1', dir)
dir_parts <- unlist(strsplit(dir, "5_"))
dir_home <- dir_parts[1]
output_dir <- paste0(dir_home ,'6_', species,'_plots')
clusters <- sort(unique(sub("\\..*", "", list.files(dir))))

# Column names
genomecov_col <- c("chr", "chrPos", "depth")
bed_col <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

# Empty lists, to be filled and transformed as a single data frame
chr <- chrPos <- depth_total <- depth_5 <- depth_3 <- centroid <- vector()

# Read each file and fill lists
for (cluster in clusters) {
  total_genomecov <- read.table(file.path(dir, paste0(cluster, '.genomecov')), sep = '\t', col.names = genomecov_col, header = FALSE)
  total_genomecov_5 <- read.table(file.path(dir, paste0(cluster, '.genomecov.5')), sep = '\t', col.names = genomecov_col, header = FALSE)
  total_genomecov_3 <- read.table(file.path(dir, paste0(cluster, '.genomecov.3')), sep = '\t', col.names = genomecov_col, header = FALSE)
  
  total_genomecov <- total_genomecov[total_genomecov$depth != 0, ]
  total_genomecov_5 <- total_genomecov_5[total_genomecov_5$depth != 0, ]
  total_genomecov_3 <- total_genomecov_3[total_genomecov_3$depth != 0, ]

  
  df_partial <- merge(total_genomecov, total_genomecov_5, by = c('chr', 'chrPos'), all.x = TRUE)
  df_temp <- merge(df_partial, total_genomecov_3, by = c('chr', 'chrPos'), all.x = TRUE)
  names(df_temp) <- c('chr', 'chrPos', 'depth_total', 'depth_5', 'depth_3')
  df_temp$centroid <- as.numeric(cluster)
  
  chr <- c(chr, df_temp$chr)
  chrPos <- c(chrPos, df_temp$chrPos)
  depth_total <- c(depth_total, df_temp$depth_total)
  depth_5 <- c(depth_5, df_temp$depth_5)
  depth_3 <- c(depth_3, df_temp$depth_3)
  centroid <- c(centroid, df_temp$centroid)
  
  bed <- read.table(file.path(dir, paste0(cluster, '.bed')), sep = '\t', col.names = bed_col, header = TRUE)
  bed_reads <- bed[, 'name']
}




# Create the final data frame
data <- data.frame(chr = chr, chrPos = chrPos, depth_total = depth_total, depth_5 = depth_5, depth_3 = depth_3, centroid = as.numeric(centroid))
data <- data[order(data$centroid),]


pdf_file <- paste0(output_dir,"/Plot.pdf")

# Start a new PDF device for saving plots
pdf(pdf_file)

# Create a list of unique centroid values
centroid_values <- unique(data$centroid)

# Set up a loop to create and save plots for each centroid
for (i in seq_along(centroid_values)) {
  # Subset the data for the current centroid value
  centroid_data <- data[data$centroid == centroid_values[i], ]
  
  # Set common plot parameters
  min.depth <- 0
  max.depth <- max(c(centroid_data$depth_total, centroid_data$depth_5, centroid_data$depth_3), na.rm = TRUE)
  
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  # Create the plot
  plot(centroid_data$chrPos, centroid_data$depth_total, type = "l", lty = 2,
       col = 'grey20', xlab = "Position on mtDNA", ylab = "Depth",
       las = 2, main = paste("Cluster n°:", centroid_values[i]),
       ylim = c(min.depth, max.depth), xaxt = 'n')
  rect(centroid_data$chrPos - 0.45, 0, centroid_data$chrPos + 0.45, centroid_data$depth_5, border = "#DC3220", col = "#DC322080", las = 2)
  rect(centroid_data$chrPos - 0.45, 0, centroid_data$chrPos + 0.45, centroid_data$depth_3, border = "#005AB5", col = "#005AB580", las = 2)
  axis(side = 1, at = seq(min(centroid_data$chrPos), max(centroid_data$chrPos), by = 1), las = 2)
  legend('topright', legend = c("5' clusters ", "3' clusters"), inset=c(-0.3,0), 
         fill = c("#DC3220", "#005AB5"), cex = .8, box.lty = 0)
  
}

# Close the PDF device
garbage <- dev.off()
