### SETUP ###
# Load all required packages
required_packages <- c("deepG", "ggplot2", "microseq", "data.table", "seqinr", 
                       "dplyr", "caret", "keras", "magrittr", "patchwork", 
                       "ggseqlogo", "openxlsx", "zoo", 
                       "reticulate", "future", "future.apply", "rtracklayer")
invisible(lapply(required_packages, library, character.only = TRUE))

### Hyperparameters ###
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
target_df <- read.csv("sporulation/sporeinfo.csv")
maxlen <- 1e6
window.size <- 0.25 # start with 40% of the sequence
total.sub <- round(window.size * maxlen)
seg.len <- 1000 # 1000 bp ~ avg. length of a gene
seg.num <- total.sub / seg.len
iterations <- 100

# define functions
give_window <- function(i) {
  # Randomly select windows of seg.len from the sequence
  samps <- sample(1:(maxlen - seg.len), seg.num, replace = FALSE)
  idx_matrix <- do.call(rbind, lapply(samps, function(s) s:(s + seg.len)))
  mask <- array(0, dim = c(1, maxlen, 4))
  mask[1, idx_matrix, ] <- 1
  new_instance <- onehot_instance * mask
  pred <- spore.model %>% predict(new_instance, verbose = 0)
  if (length(pred) > 1) {
    return(list(samps = samps, pred = pred[1, 2]))
  } else {
    return(NULL)
  }
}
trim_window <- function(q) {
  # Remove q% of the best_samps randomly
  num_to_keep <- round((1 - q) * length(best_samps))  # Number of samples to keep
  samps <- sample(best_samps, num_to_keep, replace = FALSE)
  idx_matrix <- do.call(rbind, lapply(samps, function(s) s:(s + seg.len)))
  mask <- array(0, dim = c(1, maxlen, 4))
  mask[1, idx_matrix, ] <- 1
  new_instance <- onehot_instance * mask
  pred <- spore.model %>% predict(new_instance, verbose = 0)
  if (length(pred) > 1) {
    return(list(samps = samps, pred = pred[1, 2]))
  } else {
    return(NULL)
  }
}
extract_sequence <- function(instance_sub, samp, samp_end) {
  return(substr(instance_sub, samp, samp_end))
}

set.seed(123)

#### PROCESSING MULTIPLE FILES ####
# Function to process each file
process_file <- function(fasta_path, maxlen, spore.model) {
  # Read FASTA file and create one-hot encoding
  instance <- microseq::readFasta(fasta_path)$Sequence[1]
  instance_sub <- substr(instance, 1, maxlen)
  sequence_length <- nchar(instance)
  if (sequence_length < maxlen) {
    instance_sub <- paste0(instance_sub, paste(rep("N", maxlen - sequence_length), collapse = ""))
  }
  onehot_instance <- seq_encoding_label(char_sequence = instance_sub,
                                        maxlen = maxlen,
                                        start_ind = 1,
                                        vocabulary = c("A", "C", "G", "T"))
  #### GIVING WINDOWS ####
  # Run iterations
  results <- future_lapply(1:iterations, give_window, future.seed = TRUE)
  predictions <- unlist(lapply(results, function(res) if (!is.null(res)) res$pred))
  # Sanity check with skipping condition
  max_prediction <- max(predictions)
  if (max_prediction < 0.75) {
    cat("Skipping file:", fasta_path, "\n")
    return(NULL)  # Exit this function early and skip to the next file
  }
  cat("Max prediction:", max_prediction, "\n")
  
  best_samps <- results[[which.max(predictions)]]$samps
  
  #### TRIMMING WINDOWS ####
  max_pred <- Inf
  q <- 0
  list_results <- list()
  
  # Continue trimming until max_pred drops below 0.75
  while (max_pred > 0.75) {
    q <- q + 0.05
    print(paste0("Removing q: ", q * 100, "%"))
    
    results <- future_lapply(1:iterations, function(x) trim_window(q), future.seed = TRUE)

    predictions <- unlist(lapply(results, function(res) if (!is.null(res)) res$pred))
    max_pred <- max(predictions, na.rm = TRUE)
    list_results[[as.character(q)]] <- list(predictions = predictions, max_pred = max_pred, samps = unlist(lapply(results, function(res) if (!is.null(res)) res$samps)))
    
    print(paste("Max pred:", max_pred))
  }
  
  #### EXTRACT BEST SAMPS BASED ON SECOND LAST Q ####
  q_values <- as.numeric(names(list_results))
  second_last_q <- q_values[length(q_values) - 1]
  # Extract the predictions for the second last q
  predictions <- list_results[[as.character(second_last_q)]]$predictions
  num_to_keep <- round((1 - second_last_q) * length(best_samps))
  index <- which(predictions == max(predictions))
  best_samps <- list_results[[as.character(second_last_q)]]$samps[(num_to_keep * (index - 1) + 1):(num_to_keep * index)]
  #### SYNTHETIC SEQUENCE ####
  windows <- data.table(samp = best_samps, samp_end = best_samps + seg.len)
  setkey(windows, samp, samp_end)
  # Extract sequences and generate synthetic sequence
  windows[, sequence := mapply(extract_sequence, instance_sub, samp, samp_end)]
  # if a sequence is only containing N, remove it
  windows <- windows[!grepl("^N+$", windows$sequence)]
  synthetic <- paste(windows$sequence, collapse = "")
  # Calculate nucleotide percentages
  nucleotide_percentages <- prop.table(table(strsplit(synthetic, "")[[1]])) * 100
  # coverage: length(best_samps) * seg.len / maxlen
  coverage <- nrow(windows) * seg.len / maxlen
  cat("Coverage:", coverage, "\n")
  return(list(synthetic = synthetic, sequence = windows$sequence, predictions = predictions,
              nucleotide_percentages = nucleotide_percentages, coverage = coverage))
}

#### MAIN FUNCTION FOR MULTIPLE FILES ####
process_all_files <- function(target_df, maxlen, spore.model, sample.count=NULL) {
  # Filter the files with ability_TRUE == 1
  files_to_process <- target_df %>% filter(ability_TRUE == 1)
  # if sample.count is not NULL, sample the files
  if (!is.null(sample.count)) {
    files_to_process <- files_to_process %>% sample_n(sample.count)
  }
  results <- lapply(files_to_process$file, function(fasta_file) {
    fasta_path <- file.path("sporulation/genome", fasta_file)
    cat("Processing file:", fasta_file, "\n")
    result <- process_file(fasta_path, maxlen, spore.model)
    return(list(file_name = fasta_file, result = result))
  })
  
  return(results)
}

#### EXECUTION ####
plan(multicore, workers = 15)
all_results <- process_all_files(target_df, maxlen, spore.model, sample.count = 40)
saveRDS(all_results, "permutation/all_results_40.rds")

#### POST-PROCESSING ####
# load rds
# all_results <- readRDS("permutation/all_results.rds")

# filter for results != NULL
filtered_results <- Filter(function(x) !is.null(x$result), all_results)

# Step 2: Extract the necessary information and concatenate them into a data table
concatenated_results <- rbindlist(lapply(filtered_results, function(x) {
  data.table(
    file_name = x$file_name,
    synthetic = x$result$synthetic,
    max_prediction = max(x$result$predictions),
    nucleotide_percentages = as.array(x$result$nucleotide_percentages),
    coverage = x$result$coverage
  )
}), fill = TRUE)

# get all x$result$sequence and save to a vector
segments <- unlist(lapply(filtered_results, function(x) x$result$sequence))
ggseqlogo(segments, namespace=c("A", "C", "G", "T")) +
  scale_x_continuous(breaks = c(1, seq(100, 1000-1, by = 100)))

boxplot(concatenated_results$nucleotide_percentages.N ~ concatenated_results$nucleotide_percentages.V1, 
        main = "Distribution of Nucleotide Frequencies (window seq)",
        xlab = "Nucleotide", 
        ylab = "Frequency (%)")

## check ACGT distribution in the full sequences of the selected samples
files <- unlist(lapply(filtered_results, function(x) x$file_name))
full_sequences <- future_lapply(files, function(fasta_file) {
  fasta_path <- file.path("sporulation/genome", fasta_file)
  instance <- microseq::readFasta(fasta_path)$Sequence[1]
  # Only take the first maxlen characters
  instance <- substr(instance, 1, maxlen)
  return(instance)
}, future.seed = TRUE)  # Ensures reproducibility

# Function to calculate nucleotide frequencies
calculate_freq <- function(sequence) {
  nucleotides <- c("A", "C", "G", "T")  # Only consider A, C, G, T
  seq_split <- unlist(strsplit(sequence, ""))
  
  # Predefine frequency table for ACGT, fill missing with 0
  freq_table <- prop.table(table(factor(seq_split, levels = nucleotides))) * 100
  return(freq_table)
}

# Apply the function to each sequence in the list
freq_list <- lapply(full_sequences, calculate_freq)

# Combine all frequency tables into a data frame
freq_matrix <- do.call(rbind, freq_list)

# Draw the boxplot
boxplot(freq_matrix, 
        main = "Distribution of Nucleotide Frequencies (full seq)",
        xlab = "Nucleotide", 
        ylab = "Frequency (%)")
plan(sequential)

## save segments to a file
saveRDS(segments, "permutation/segments.rds")
