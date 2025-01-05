### Motif Detection

# setup
library(deepG)
library(ggplot2)
library(microseq)
library(seqinr)
library(dplyr)
library(caret)
library(pROC)
library(keras)
library(magrittr)
library(patchwork)
library(ggseqlogo)
set.seed(42)
library(reticulate)
use_python("E:/miniconda/envs/r-reticulate/python.exe", required = TRUE)
library(innsight)
source("genepermutation.R")
source("ig_modified.R")
source("seqsyn.R")

### Experiment Two: Different Entropy levels
original_seq <- "TATAGCGCAGCTGTCAGCATGAGTCCATGA"

# Function to shuffle sequences
shuffle_seq <- function(seq) {
  paste(sample(unlist(strsplit(seq, "")), length(unlist(strsplit(seq, "")))), collapse = "")
}

# 1. High Entropy (Randomized)
high_entropy <- shuffle_seq(original_seq)

# 2. Moderately High Entropy (Shuffled Regions)
mod_high_entropy <- paste(shuffle_seq(substr(original_seq, 1, 10)),
                          shuffle_seq(substr(original_seq, 11, 20)),
                          shuffle_seq(substr(original_seq, 21, 30)),
                          sep = "")

# 3. Medium Entropy (Partial Repeats)
medium_entropy <- paste(substr(original_seq, 1, 10), shuffle_seq(substr(original_seq, 11, 20)),
                        substr(original_seq, 21, 30), sep = "")

# 4. Moderate to Low Entropy (Increasing Repeats)
mod_low_entropy <- paste(substr(original_seq, 1, 15), substr(original_seq, 1, 15), sep = "")

# 5. Low Entropy (Major Repeats)
low_entropy <- paste(rep(substr(original_seq, 1, 5), 6), collapse = "")

# 6. Very Low Entropy (Single Repeat)
very_low_entropy <- strrep("TA", 15)

# 7. Lowest Entropy (Homogeneous Sequence)
lowest_entropy <- strrep("A", 30)

# Collect sequences into a list for further use
motif.vec <- c(
  high_entropy = high_entropy,
  mod_high_entropy = mod_high_entropy,
  medium_entropy = medium_entropy,
  mod_low_entropy = mod_low_entropy,
  low_entropy = low_entropy,
  very_low_entropy = very_low_entropy,
  lowest_entropy = lowest_entropy
)

nsamp <- 10000
lenseq <- 600

for (length_name in names(motif.vec)) {
  motif <- motif.vec[[length_name]]
  synthetic <- SeqSyn(lenseq, nsamp, codon.dict, by.codon=FALSE)
  synthetic <- synthetic[vapply(synthetic, function(seq) {
    is.null(find_motif(seq, motif))
  }, logical(1))]
  n <- length(synthetic)
  selected_indices <- sample(1:n, size = floor(n / 2), replace = FALSE)
  for (i in selected_indices) {
    synthetic[[i]] <- insert_motif(synthetic[[i]], motif)
  }
  labels <- ifelse(1:n %in% selected_indices, "motif", "normal")
  df <- data.frame(sequence = I(synthetic), label = labels, stringsAsFactors = FALSE)
  print(paste0(length_name, ": ", n, " samples"))
  print(table(df$label))
  file_name <- sprintf("motif/ent/synthetic_%s.csv", length_name)
  write.csv(df, file_name, row.names = FALSE)
}

## data processing
library(caret)

# Specify the directory containing your CSV files
data_directory <- "motif/ent"

# Retrieve all CSV files in the directory
file_paths <- list.files(data_directory, pattern = "\\.csv$", full.names = TRUE)

# Loop over the file paths, load, split, train, and predict
for (file_path in file_paths) {
  # Extract file name without extension for naming output files later
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Load the data from the file
  dataset <- read.csv(file_path)
  
  # Split the data into training, testing, and validation sets
  set.seed(123) # Ensure reproducibility
  # Create indices for creating training (60%), testing (20%), and validation (20%) sets
  trainIndex <- createDataPartition(dataset$label, p = 0.6, list = FALSE, times = 1)
  trainSet <- dataset[trainIndex, ]
  testValidSet <- dataset[-trainIndex, ]
  testIndex <- createDataPartition(testValidSet$label, p = 0.5, list = FALSE, times = 1)
  testSet <- testValidSet[testIndex, ]
  validSet <- testValidSet[-testIndex, ]
  
  # Directory Structure Setup
  main_dir <- "motif"
  sub_dirs <- c("normal", "motif")
  data_sets <- c("train", "validation", "test")
  
  # Create directories
  dir.create(file.path(main_dir, file_name), showWarnings = FALSE, recursive = TRUE)
  for (sub in sub_dirs) {
    dir.create(file.path(main_dir, file_name, sub), showWarnings = FALSE)
    for (data_set in data_sets) {
      dir.create(file.path(main_dir, file_name, sub, data_set), showWarnings = FALSE)
    }
  }
  
  # Function to write each row into a FASTA file in the corresponding folder, according to label
  write_fasta <- function(data_set, set_name) {
    for (i in 1:nrow(data_set)) {
      file_path <- file.path(main_dir, file_name, data_set$label[i], set_name, paste0(i, ".fasta"))
      seqinr::write.fasta(sequences = as.list(data_set$sequence[i]), 
                          names = as.character(i), 
                          file.out = file_path)
    }
  }
  
  # Apply the function to each dataset
  write_fasta(trainSet, "train")
  write_fasta(validSet, "validation")
  write_fasta(testSet, "test")
}

## model training

evaluation_results <- list()

for (file_path in file_paths) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  message(paste0("Processing ", file_name))
  model <- create_model_lstm_cnn(
    maxlen = lenseq,
    layer_lstm = NULL,
    layer_dense = c(2L),
    vocabulary_size = 4,
    kernel_size = c(5, 7, 9),  # Smaller and varied filter sizes
    filters = c(32, 64, 128),
    pool_size = c(3, 3, 3),    # Smaller pooling sizes
    learning_rate = 0.0005     # Slightly increased learning rate
  )
  path_checkpoint <- file.path("motif/checkpoints")
  dir_path <- file.path("motif/outputs")
  
  path_normal_train <- file.path("motif", file_name, "normal/train")
  path_normal_validation <- file.path("motif", file_name, "normal/validation")
  path_motif_train <- file.path("motif", file_name, "motif/train")
  path_motif_validation <- file.path("motif", file_name, "motif/validation")
  path_normal_test <- file.path("motif", file_name, "normal/test")
  path_motif_test <- file.path("motif", file_name, "motif/test")
  
  hist <- train_model(
    train_type = "label_folder",
    model = model,
    path = c(path_normal_train, path_motif_train),
    path_val = c(path_normal_validation, path_motif_validation),
    vocabulary_label = c("normal", "motif"),
    path_checkpoint = path_checkpoint,
    train_val_ratio = 0.2,
    run_name = file_name,
    batch_size = 64,          # Increased batch size
    steps_per_epoch = 30,
    epochs = 15                # Increased number of epochs
  )
  plot(hist)
  # save hist
  saveRDS(hist, file = paste0("motif/outputs/", file_name, ".RDS"))
  
  # Evaluate the model on the test data
  eval_model <- evaluate_model(
    path_input = c(path_normal_test, path_motif_test),
    model = model,
    batch_size = 64,
    step = 10,
    vocabulary_label = list(c("normal", "motif")),
    number_batches = 10,
    mode = "label_folder",
    verbose = FALSE
  )
  print(file_name)
  print("Model Evaluation: \n")
  print(eval_model)
  evaluation_results[[file_name]] <- eval_model
}

saveRDS(evaluation_results, file = "motif/outputs/evaluation_results_ent.RDS")

accuracies <- numeric(length(evaluation_results))
matches <- c()
# Extract accuracy for each model and store in the vector
for (i in seq_along(evaluation_results)) {
  accuracies[i] <- evaluation_results[[i]][[1]]$accuracy
}
acc_df <- data.frame(Motif = rev(names(motif.vec)), Accuracy = accuracies)

### Integrated Gradients
# Define a function to perform analysis
perform_analysis <- function(entropy_level) {
  set.seed(123)
  # Load the model
  model_path <- paste0("motif/checkpoints/synthetic_", entropy_level)
  model <- load_cp(model_path, cp_filter = "last_ep")
  
  # Load the motif
  motif <- motif.vec[[entropy_level]]
  
  # List and sample files
  folder <- paste0("motif/synthetic_", entropy_level, "/motif/test")
  files_in_folder <- list.files(folder, full.names = TRUE)
  sample_file <- sample(files_in_folder, 1)
  instance <- microseq::readFasta(sample_file)$Sequence[1]
  
  # Encode the sequence
  onehot_instance <- seq_encoding_label(char_sequence = instance,
                                        maxlen = 600,
                                        start_ind = 1,
                                        vocabulary = c("A", "C", "G", "T"))
  
  # Find motif position
  motif_pos <- find_motif(instance, motif)
  motif_end <- motif_pos + nchar(motif) - 1
  
  # Create a baseline
  onehot_baseline_25 <- onehot_instance * 0 + 0.25
  
  # Predict
  pred <- predict(model, onehot_instance, verbose = 0)
  print(pred[2])
  # Compute Integrated Gradients
  ig <- ig_modified(m_steps = 400,
                    input_seq = onehot_instance,
                    baseline_type = "modify",
                    baseline_onehot = onehot_baseline_25,
                    target_class_idx = 2,
                    model = model,
                    num_baseline_repeats = 1)
  sum <- rowSums(as.array(ig))
  abs_sum <- rowSums(abs(as.array(ig)))
  df25 <- data.frame(abs_sum = abs_sum, sum = sum, position = 1:600)
  print(motif_pos)
  print(motif_end)
  # Plot and save
  p <- ggplot(df25, aes(x = position, y = sum)) +
    geom_rect(aes(xmin = motif_pos, xmax = motif_end, ymin = -Inf, ymax = Inf),
              fill = "lightblue", alpha = 0.2) + geom_point() + theme_bw() +
    labs(subtitle = paste("Baseline 0.25 for", entropy_level))
  
  ggsave(filename = paste0("motif/outputs/", entropy_level, "_plot.svg"),
         plot = p, device = "svg", width = 10, height = 6)
  
  # Convert IG to matrix and prepare for ggseqlogo
  igmat <- as.data.frame(t(as.matrix(ig)))
  rownames(igmat) <- c("A", "C", "G", "T")
  igmat <- as.matrix(igmat)
  
  # Generate ggseqlogo plot
  p_ig <- ggseqlogo(igmat, method='custom', seq_type='dna') +
    xlim(motif_pos - 20, motif_end + 20) +
    labs(x = "bp", y = "IG", subtitle = paste0(entropy_level, ": ", motif))
  
  # Save ggseqlogo plot
  ggsave(filename = paste0("motif/outputs/", entropy_level, "_ig_plot.svg"),
         plot = p_ig, device = "svg", width = 6, height = 4)
  # Return results for saving
  list(ig = ig, df25 = df25, plot = p)
}

# Vector of entropy levels
entropy_levels <- names(motif.vec)
enrtopy_level <- "lowest_entropy"
# Initialize lists to store results
ig_list <- list()
df_list <- list()
plot_list <- list()

# Loop over entropy levels and perform analysis
for (level in entropy_levels) {
  results <- perform_analysis(level)
}

averaged_match <- function(entropy_level) {
  set.seed(123)  # Ensure reproducibility for sampling
  
  # Load the model
  model_path <- paste0("motif/checkpoints/synthetic_", entropy_level)
  model <- load_cp(model_path, cp_filter = "last_ep")
  
  # Load the motif
  motif <- motif.vec[[entropy_level]]
  
  # List files
  folder <- paste0("motif/synthetic_", entropy_level, "/motif/test")
  files_in_folder <- list.files(folder, full.names = TRUE)
  
  # Initialize vector to store match rates
  matches <- numeric(100)
  
  for (i in 1:100) {
    # Randomly sample a file without replacement
    sample_file <- sample(files_in_folder, 1)
    files_in_folder <- setdiff(files_in_folder, sample_file)  # Remove the sampled file
    
    # Load and process instance
    instance <- microseq::readFasta(sample_file)$Sequence[1]
    onehot_instance <- seq_encoding_label(char_sequence = instance,
                                          maxlen = 600,
                                          start_ind = 1,
                                          vocabulary = c("A", "C", "G", "T"))
    
    # Find motif position and calculate IG
    motif_pos <- find_motif(instance, motif)[1]
    motif_end <- motif_pos + nchar(motif) - 1
    onehot_baseline_25 <- onehot_instance * 0 + 0.25
    ig <- ig_modified(m_steps = 400,
                      input_seq = onehot_instance,
                      baseline_type = "modify",
                      baseline_onehot = onehot_baseline_25,
                      target_class_idx = 2,
                      model = model,
                      num_baseline_repeats = 1)
    
    # Extract motif region and analyze match
    igmat <- as.data.frame(t(as.matrix(ig)))
    rownames(igmat) <- c("A", "C", "G", "T")
    motifregion <- igmat[, motif_pos:(motif_end)]  # Adjust to ensure proper column indexing
    row_max_in_cols <- apply(motifregion, 2, function(col) rownames(igmat)[which.max(col)])
    interpret <- unname(row_max_in_cols)
    motif_letters <- strsplit(motif, split = "")[[1]]
    matches[i] <- sum(interpret == motif_letters) / length(motif_letters)
  }
  
  return(matches)
}

# Vector of entropy levels
entropy_levels <- names(motif.vec)

# DataFrame to store all match vectors
matches_df <- setNames(as.data.frame(matrix(ncol = length(entropy_levels), nrow = 100)),
                       entropy_levels)

# Perform analysis for each level and populate the DataFrame
for (level in entropy_levels) {
  matches_df[, level] <- averaged_match(level)
}

matches_df

library(tidyr)   # For data transformation
library(ggplot2) # For plotting

# Assuming matches_df is already created and populated

matches_long <- matches_df %>%
  pivot_longer(cols = everything(), names_to = "Entropy_Level", values_to = "Match")
entropy_order <- c("lowest_entropy", "very_low_entropy", "low_entropy", "mod_low_entropy", "medium_entropy", "mod_high_entropy", "high_entropy")

matches_long$Entropy_Level <- factor(matches_long$Entropy_Level, levels = entropy_order)

# Create a boxplot
ggplot(matches_long, aes(x = Entropy_Level, y = Match)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  theme_bw() +
  labs(title = "Average Nucleotides Match Values by Entropy Level (100 samples)",
       x = "Entropy Level",
       y = "Match Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("motif/outputs/match_boxplot.svg", width = 10, height = 6, units = "in")

