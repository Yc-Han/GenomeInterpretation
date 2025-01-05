# define functions
# Function to calculate nucleotide frequencies
nt_density <- function(gene_seq) {
  nt_bases <- c("A", "C", "G", "T")
  # if "N" is present in the sequence, remove it
  gene_seq <- gsub("N", "", gene_seq)
  gene_chars <- unlist(strsplit(gene_seq, NULL))
  gene_factor <- factor(gene_chars, levels = nt_bases)
  nt_counts <- table(gene_factor)
  pop.dense <- nt_counts / length(gene_chars)
  pop.dense <- as.numeric(pop.dense)
  return(pop.dense)
}
pop.dense <- c(0.25, 0.25, 0.25, 0.25)

give_window <- function(i) {
  # Randomly select windows of seg.len from the sequence
  samps <- sample(1:(maxlen - seg.len), seg.num, replace = FALSE)
  idx_matrix <- do.call(rbind, lapply(samps, function(s) s:(s + seg.len)))
  new_instance <- onehot_instance
  new_instance[1, idx_matrix, ] <- matrix(rep(pop.dense, length(idx_matrix)),
                                          nrow = length(idx_matrix), byrow = TRUE)
  pred <- model %>% predict(new_instance, verbose = 0)
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
  mask <- array(1, dim = c(1, maxlen, 4))
  mask[1, idx_matrix, ] <- 0
  new_instance <- onehot_instance * mask
  new_instance[mask == 0] <- 0.25
  pred <- model %>% predict(new_instance, verbose = 0)
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
process_file <- function(fasta_path, maxlen, model) {
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
  cat("Max false prediction:", max_prediction, "\n")
  
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
    
    print(paste("Max false pred:", max_pred))
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
process_all_files <- function(target_df, maxlen, model, sample.count=NULL) {
  # Filter the files with ability_TRUE == 1
  files_to_process <- target_df %>% filter(ability_TRUE == 1)
  # if sample.count is not NULL, sample the files
  if (!is.null(sample.count)) {
    files_to_process <- files_to_process %>% sample_n(sample.count)
  }
  results <- lapply(files_to_process$file, function(fasta_file) {
    fasta_path <- file.path("sporulation/genome", fasta_file)
    cat("Processing file:", fasta_file, "\n")
    result <- process_file(fasta_path, maxlen, model)
    return(list(file_name = fasta_file, result = result))
  })
  
  return(results)
}