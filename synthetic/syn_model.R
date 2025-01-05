library(tidyverse)
library(keras)
library(reticulate)
# use_condaenv("tf", required = TRUE)
library(deepG)
library(keras)
library(magrittr)
library(data.table)
set.seed(42)

filters <- c(32,32,64,64,128,128,256)

model <- create_model_lstm_cnn(
  maxlen = 1000000,
  layer_lstm = NULL,
  layer_dense = c(256,2),
  learning_rate = 0.000001,
  filters = filters,
  kernel_size = rep(24, length(filters)),
  pool_size = rep(4, length(filters)))

target_from_csv <- "synthetic/file_labels.csv"
target_df <- read.csv(target_from_csv)
label_names <- names(target_df)[names(target_df) != "file"]
head(target_df)
print(label_names)
train_path <- "synthetic/data/train"
val_path <- "synthetic/data/validation"
check_path <- "synthetic/checkpoints"

hist <- train_model(
  train_type = "label_csv",
  target_from_csv = target_from_csv,
  model = model,
  path = train_path,
  path_val = val_path,
  # path_checkpoint = check_path,
  train_val_ratio = 0.5,
  run_name = "synthetic_test",
  batch_size = 5,
  epochs = 2,
  steps_per_epoch = 10,
  max_samples = 5,
  format = "fasta",
  concat_seq = "",
  vocabulary_label = label_names)

plot(hist)