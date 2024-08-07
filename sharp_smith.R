#!/usr/bin/env Rscript


# # default arguments
args <- list(
  mode = "score",
  t_five_score = 0.50,
  t_three_score = 0.00,
  penalty = 0.10,
  n_thre = 0.50,
  path_bedfiles = NULL,
  path_smith = NULL,
  sex = "NA"
)

# in order to execute from example (script in $PATH) sensible defoults are needed for path_bedfiles and path_smith
# these have to be relative to example
baserunfolder = getwd() # run folder (example)
alldirs = list.dirs(getwd(), recursive = TRUE) # all directories
args$path_bedfiles = alldirs[ grepl('*.results.clusters.bedfiles', alldirs) ] # correct bedfiles dir, defoult
dir_smith = alldirs[ grepl('*_smithRNAs', alldirs) ] 
args$path_smith = paste(dir_smith, 'presumptive_smithRNAs.fa', sep = '/') # correct smith dir/file, defoult


print_help <- function() {
  cat("Usage: Rscript script_name.R [options]\n")
  cat("Options:\n")
  cat("  --mode=<mode>                  Mode of operation: filter, score, list. Default is 'score'.\n")
  cat("  --t_five_score=<value>         T-five score threshold. Default is 0.50.\n")
  cat("  --t_three_score=<value>        T-three score threshold. Default is 0.00.\n")
  cat("  --penalty=<value>              Penalty value. Default is 0.10.\n")
  cat("  --n_thre=<value>               N threshold value. Default is 0.50.\n")
  cat("  --path_bedfiles=<path>         Path to bedfiles. Default is 5_.../5.2_... within the current folder.\n")
  cat("  --path_smith=<path>            Path to Smith file. Default is 7_.../presumptive_smithRNAs.fa within the current folder.\n")
  cat("  --sex=<value>                  Sex information. Default is 'NA'.\n")
  cat("  --help, -h                     Show this help message and exit.\n")
}

# Parse command-line arguments if provided
args_list <- strsplit(commandArgs(trailingOnly = TRUE), "=")
for (arg in args_list) {
  match <- gsub('--', '', arg[1])
  if (match %in% c("help", "h")) {
    print_help()
    quit(save = "no", status = 0)
  }
  if (match %in% names(args)) {
    args[[match]] <- ifelse(length(arg) > 1, arg[2], "")
  }
}

# Check if required arguments are provided
if (is.null(args$path_bedfiles) || args$path_bedfiles == "") {
  stop("Error: 'path_bedfiles' argument is required.")
}

if (args$mode == "filter" && (is.null(args$path_smith) || args$path_smith == "")) {
  stop("Error: 'path_smith' argument is required when mode is 'filter'.")
}

# rename arguments in objects
path_bedfiles = args$path_bedfiles
penalty = args$penalty
t_five_score = args$t_five_score
t_three_score = args$t_three_score
n_thre = args$n_thre
path_smith = args$path_smith
sex = args$sex
mode = args$mode

# restrict the span on observed peaks
restrict_range <- function(input) {
  start_range <- NULL
  end_range <- NULL
  
  # pick up the non 0 coverages indexes (rownames)
  above_zeros = rownames(input[input$cov!=0,])
  # get the start range, aka the rowname from where the coverage raises
  start_range <-  above_zeros[1]
  # get the start range, aka the rowname from where the coverage drops
  end_range <- tail(above_zeros, 1)
  #subset the original dataframe with that range and returns it
  subs_data = subset(input, rownames(input) %in% start_range:end_range)
  return(subs_data)
}

# get "good" coverages
n_anta <- function(input) {
  # Order the input in descending order
  decre_input <- input[order(input$cov, decreasing = TRUE),]
  
  # cumulative sum of coverages ( the first cov will be summed to 0, the second cov will be summed with the previous and so on)
  decre_input$cumulative_sum <- cumsum(decre_input$cov)
  # calculate the threshold to the sum of everythin. It's the half of the total sums
  threshold = n_thre*sum(decre_input$cov)
  # tag those with upper threshold as passed
  decre_input$pass <- ifelse(decre_input$cumulative_sum < threshold, 'not_passed', 'passed') 
  # Find the index of the first occurrence of "passed"
  passed_index <- which(decre_input$pass == "passed")[1]  
  # and subset the dataframe accordingly (from the beginning to the first "peak")
  subset_df <- decre_input[1:passed_index, ]
  # empty rownames
  rownames(subset_df) <- NULL
 
  # Return the count of peaks and the list of peaks
  return(subset_df)
}

# 
find_gappy <- function(input1, input2){
  # merge the two results
  rownames(input1) <- NULL
  merged = merge(input1, input2, by=c('species','pos','cov'))
  duplicates <- duplicated(merged$cov, fromLast = T)
  merged <- merged[!duplicates, ]
  
  # get the lowest position value (before, in Diego script corresponded to the list index)
  maximum = max(merged$pos)
  # get the highest position value (before, in Diego script corresponded to the list index)
  minumum = min(merged$pos)
  # calculate the gap between the beginning and the end of passed coverages and ane it as GAPPY
  gappy = maximum - minumum -1
  return(gappy)
}

# list all coverages files
files = list.files(path=path_bedfiles)
# and make absolute paths
file_paths <- file.path(path_bedfiles, files)


# create an empty list of clusters
cluster_list <- vector("character", length(file_paths))

# Loop through each cluster path
for (i in seq_along(file_paths)) {
  # Extract filename from full path
  file_name <- basename(file_paths[i])
  # Split filename by '.' and select the first part and add it to the list
  cluster_list[i] <- unlist(strsplit(file_name, "\\."))[1]
}

# make unique names
cluster_list = unique(cluster_list)

# create the final printed dataframe (useful later)
final_df=structure(list(cluster = character(), five_score=integer(),three_score=integer()),class = "data.frame")

# Loop through each cluster in cluster_list
for (cluster in cluster_list) {
  # Read .genomecov.5 file and extract second column
  five_end_data <- read.table(paste0(path_bedfiles, '/', cluster, ".genomecov.5"), header = FALSE, sep = "\t", col.names = c('species','pos','cov'))

  # Read .genomecov.3 file and extract second column
  three_end_data <- read.table(paste0(path_bedfiles, '/', cluster, ".genomecov.3"), header = FALSE, sep = "\t",col.names = c('species','pos','cov'))

  # apply the functions to the 5' and 3'
  five_end_res = restrict_range(five_end_data)
  n_anta_out_five = n_anta(five_end_res)
  five_score = 1/(nrow(n_anta_out_five) + (penalty*find_gappy(five_end_res,n_anta_out_five)))
  
  three_end_res = restrict_range(three_end_data)
  n_anta_out_three = n_anta(three_end_res)
  three_score = 1/(nrow(n_anta_out_three) + (penalty*find_gappy(three_end_res,n_anta_out_three)))
  
  # merge everything in the final df
  final_df <- rbind(final_df, data.frame(cluster=cluster, five_score = five_score, three_score = three_score))
}

# add the sex column
final_df$sex <- sex


# if mode score print the scores
if (mode == 'score'){
  cat("\n")
  print(final_df, row.names = F)
} else{ # otherwise tell the user which clusters passed or didn't
  list_df <- final_df
  list_df$passed <- ifelse(list_df$five_score >  t_five_score & three_score > t_three_score, 'Y', 'N')
  list_df <- list_df[,-c(2,3)]
  if (mode == 'list'){ # in case the mode is set lo list, then print it
    cat("\n")
    print(list_df, row.names = F)}
  if (mode == 'filter' ){ # otherwise write a tmp file to filter fasta files
    # make the df
    list_df = list_df[list_df$passed == 'Y',]
    # write.table(list_df, file = paste0(path_bedfiles,'/table.tpm'), quote = F, sep='\t', row.names = F)
    cat("\n")
    print(list_df,row.names = F)
    # pick up the interesting clusters
    ids <- list_df[,1]
    # create a new smith file but rename the old one
    file.rename(path_smith, paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", path_smith), '_old.fa'))
    path_smith_old = paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", path_smith), '_old.fa')
    path_smith_new = path_smith
    output_file <- file(path_smith_new, open = "wt")
    # read the old smith file
    path_smith_old <- readLines(path_smith_old)
    
    for (i in ids) {
      # Construct the pattern to search for
      pattern <- paste0("clusterid", i, "_")
      
      # Find lines that match the pattern
      matching_indices <- grep(pattern, path_smith_old)
      
      # Write the matching lines and the following line (sequence) to the output file
      for (index in matching_indices) {
        # Write the header line
        writeLines(path_smith_old[index], output_file)
        # Write the sequence line (the following line)
        if (index + 1 <= length(path_smith_old)) {
          writeLines(path_smith_old[index + 1], output_file)
        }
      }
    }
    
    # Close the output file
    close(output_file)}
}
