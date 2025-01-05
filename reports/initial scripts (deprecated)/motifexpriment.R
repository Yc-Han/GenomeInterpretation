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

### Experiment One: Different lengths
motif.vec <- c(len33 = "TATAGCGCAGCTGTCAGCATGAGTCCATGATCA",
               len30 = "TATAGCGCAGCTGTCAGCATGAGTCCATGA",
               len27 = "TATAGCGCAGCTGTCAGCATGAGTCCA",
               len24 = "TATAGCGCAGCTGTCAGCATGAGT",
               len21 = "TATAGCGCAGCTGTCAGCATG",
               len18 = "TATAGCGCAGCTGTCAGC",
               len15 = "TATAGCGCAGCTGTC",
               len12 = "TATAGCGCAGCT",
               len09  = "TATAGCGCA")
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
  file_name <- sprintf("motif/len/synthetic_%s.csv", length_name)
  write.csv(df, file_name, row.names = FALSE)
}

## data processing
library(caret)

# Specify the directory containing your CSV files
data_directory <- "motif/len"

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
    maxlen = lenseq, # lenseq = 600
    layer_lstm = NULL,
    layer_dense = c(2L),
    vocabulary_size = 4,
    kernel_size = c(5, 7, 9),
    filters = c(32, 64, 128),
    pool_size = c(3, 3, 3),
    learning_rate = 0.0005
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

saveRDS(evaluation_results, file = "motif/outputs/evaluation_results2.RDS")

accuracies <- numeric(length(evaluation_results))

# Extract accuracy for each model and store in the vector
for (i in seq_along(evaluation_results)) {
  accuracies[i] <- evaluation_results[[i]][[1]]$accuracy
}

# The last element needs to be moved to the first position
# Rotate the vector so that the last element becomes the first
# Rotate the vector so that the last element becomes the first
# accuracies <- c(accuracies[length(accuracies)], accuracies[-length(accuracies)])

acc_df <- data.frame(Motif = rev(names(motif.vec)), Accuracy = accuracies)
ggplot(acc_df, aes(x = Motif, y = Accuracy)) +
  geom_line(group = 1, color = "blue") +  # Connect points with a line
  geom_point(size = 3, color = "red") +  # Red points
  theme_minimal() +  # Minimal theme
  labs(title = "Accuracy by Motif Length",
       x = "Motif",
       y = "Model Accuracy")
