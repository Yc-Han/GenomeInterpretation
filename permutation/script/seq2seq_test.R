generate_mask <- function(one_hot_encoded_matrix, model, maxlen = 5000) {
  if (dim(one_hot_encoded_matrix)[3] != 4) {
    stop("Input matrix must have shape [1, x, 4].")
  }
  sequence_length <- dim(one_hot_encoded_matrix)[2]
  # Initialize an empty list to store predictions
  predictions <- list()
  # Loop through the input matrix in chunks of size maxlen
  start_idx <- 1
  while (start_idx <= sequence_length) {
    # End index for the current chunk
    end_idx <- min(start_idx + maxlen - 1, sequence_length)
    # Extract the chunk [1, maxlen, 4]
    chunk <- one_hot_encoded_matrix[, start_idx:end_idx, ]
    # If the chunk is shorter than maxlen, pad with zeros
    if ((end_idx - start_idx + 1) < maxlen) {
      padded_chunk <- array(0, dim = c(1, maxlen, 4))
      padded_chunk[, 1:(end_idx - start_idx + 1), ] <- chunk
      chunk <- padded_chunk
    }
    reshaped_chunk <- array_reshape(chunk, dim = c(1, 5000, 4))
    # Predict using the model
    chunk_prediction <- model %>% predict(reshaped_chunk, verbose=0)
    # Trim padding in the output if applicable
    chunk_prediction <- chunk_prediction[, 1:(end_idx - start_idx + 1), ]
    # Append the prediction to the list
    predictions[[length(predictions) + 1]] <- chunk_prediction
    # Update the start index for the next reshaped_chunk
    start_idx <- end_idx + 1
  }
  # Concatenate all predictions into a single array
  concatenated_predictions <- do.call(abind::abind, c(predictions, list(along = 2)))
  # Ensure the output shape is [x, 1]
  concatenated_predictions <- array_reshape(concatenated_predictions,
   dim = c(dim(concatenated_predictions)[1] * dim(concatenated_predictions)[2], 1))
  # Ensure the output shape is [1, x, 1]
  return(concatenated_predictions)
}

# Example usage:
loaded_model <- create_model_seq2seq()
summary(loaded_model)
loaded_model %>% load_model_weights_hdf5("permutation/seq2seq_cps/trial2_best_model_epoch_04_val_loss_1.34.h5")

#### Set up the environment ####
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "keras", "openxlsx", "zoo", "rtracklayer")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}
install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

#### Load the data ####
maxlen <- 1000000
spore.model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
#### Single FASTA ####
# example: Bacillus subtilis, first 1 Mbp.
fasta_path <- "sporulation/genome/GCF_029537415.1_ASM2953741v1_genomic.fasta"
instance <- microseq::readFasta(fasta_path)$Sequence[1]
instance_sub <- substr(instance, 1, maxlen)
pop.dense <- nt_density(instance_sub)
sequence_length <- nchar(instance)
if (sequence_length < maxlen) {
  instance_sub <- paste0(instance_sub, paste(rep("N", maxlen - sequence_length), collapse = ""))
}
onehot_instance <-  seq_encoding_label(char_sequence = instance_sub,
                                       maxlen = maxlen,
                                       start_ind = 1,
                                       vocabulary = c("A", "C", "G", "T"))
pred <- predict(spore.model, onehot_instance)
onehot_baseline <- array(rep(nt_density(instance_sub), each = dim(onehot_instance)[2]), 
                         dim = dim(onehot_instance))

concatenated_mask <- generate_mask(onehot_instance, loaded_model, maxlen = 5000)
print(dim(concatenated_mask))
print(concatenated_mask[1:100, ])

binary_mask <- ifelse(concatenated_mask > 0.677, TRUE, FALSE)
print(sum(binary_mask))

masked_onehot <- onehot_instance
masked_onehot[, binary_mask, ] <- matrix(pop.dense, 
                                          nrow = sum(binary_mask), 
                                          ncol = length(pop.dense), 
                                          byrow = TRUE)
pred_masked <- predict(spore.model, masked_onehot)
