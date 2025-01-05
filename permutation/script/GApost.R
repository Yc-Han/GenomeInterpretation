library(keras)
library(reticulate)
library(tensorflow)
library(purrr)
source("permutation/script/gentoma.R")
spore.model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
fullmatches <- c()
count <- 0
predictions <- data.frame(high = numeric(), low = numeric(), all_low = numeric(), random = numeric(), allzero = numeric())
save_grouped_rds_batches <- function(processed_list, item.name, fragment_length = 5000,
                                     group_size = 20, output_dir = "permutation/tensors/") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  batch_counter <- 1
  grouped_inputs <- list()
  grouped_masks <- list()
  
  for (i in seq_along(processed_list)) {
    instance <- processed_list[[i]]$instance
    mask <- processed_list[[i]]$mask
    
    # One-hot encode
    encoded_instance <- deepG::seq_encoding_label(
      char_sequence = instance,
      maxlen = nchar(instance),
      start_ind = 1,
      vocabulary = c("A", "C", "G", "T")
    )
    
    # Split into fragments
    num_fragments <- ceiling(nchar(instance) / fragment_length)
    for (frag_idx in seq_len(num_fragments)) {
      start_ind <- (frag_idx - 1) * fragment_length + 1
      end_ind <- min(frag_idx * fragment_length, nchar(instance))
      input_fragment <- encoded_instance[,start_ind:end_ind, , drop = FALSE]
      mask_fragment <- mask[start_ind:end_ind]
      
      # Append to grouped lists
      grouped_inputs <- c(grouped_inputs, list(input_fragment))
      grouped_masks <- c(grouped_masks, list(mask_fragment))
      
      # Save a group when it reaches the group size
      if (length(grouped_inputs) == group_size) {
        saveRDS(
          list(inputs = grouped_inputs, masks = grouped_masks),
          file = paste0(output_dir, item.name, "_Batch_", batch_counter, ".rds")
        )
        grouped_inputs <- list()
        grouped_masks <- list()
        batch_counter <- batch_counter + 1
      }
    }
  }
  
  # Save any remaining fragments
  if (length(grouped_inputs) > 0) {
    saveRDS(
      list(inputs = grouped_inputs, masks = grouped_masks),
      file = paste0(output_dir, item.name, "_Batch_", batch_counter, ".rds")
    )
  }
}

# Process all groups in the folder
folder <- "permutation/singles/"
file.remove(list.files(folder, pattern = ".out", full.names = TRUE))  # Remove unnecessary files

files <- list.files(folder, full.names = TRUE)         
files <- split(files, gsub("(.*)\\.fasta.*", "\\1", files))  # Group files by base name

# Define the processing function for each group
process_group <- function(file_group, output_dir) {
  file_group <- unlist(file_group)
  item.list <- lapply(file_group, readRDS)
  item.name <- sub(".*singles/(.*)\\.fasta.*", "\\1", file_group[1])
  
  # Filter item.list
  item.list <- map(item.list, ~ keep(.x, ~ !is.null(.x) && .x$best_individual$eff.pred >= 0.5))
  
  # Rearrange the list
  rearranged_list <- lapply(seq_along(item.list[[1]]), function(i) {
    valid_elements <- lapply(item.list, function(x) {
      if (length(x) >= i && !is.null(x[[i]])) x[[i]] else NULL
    })
    Filter(Negate(is.null), valid_elements)
  })
  # if list is empty, to the next file
  if (length(rearranged_list) == 0) {
    cat("No valid elements in file: ", item.name, "\n")
    return()
  }
  rearranged_list <- Filter(length, rearranged_list)
  names(rearranged_list) <- paste0("Instance", seq_along(rearranged_list))
  
  # Process each instance
  processed_list <- lapply(rearranged_list, function(elements) {
    instance <- gsub("N", "", elements[[1]]$instance)  # Remove all "N"s
    if (nchar(instance) > 1e6) stop("Instance length is not 1e6 after cleaning.")
    
    instance_length <- nchar(instance)
    mask_matrix <- matrix(0, nrow = length(elements), ncol = instance_length)
    
    for (i in seq_along(elements)) {
      element <- elements[[i]]
      starts <- element$best_individual$starts
      ends <- element$best_individual$ends
      intervals <- data.frame(starts, ends)
      for (j in seq_len(nrow(intervals))) {
        start <- max(1, intervals$starts[j])
        end <- min(instance_length, intervals$ends[j])
        
        if (start <= end) {  # Only update if the interval is valid
          mask_matrix[i, start:end] <- 1
        } else {
          cat("At file: ", item.name, "\n")
          cat("Invalid interval: ", start, end, "\n")
        }
      }
    }
    mask_mean <- colMeans(mask_matrix)
    list(instance = instance, mask = mask_mean)
  })
  names(processed_list) <- paste0("Instance_", seq_along(processed_list))
  
  # Save grouped batches
  # save_grouped_rds_batches(processed_list, item.name = item.name, output_dir = output_dir)
  # save the proportion of mask_mean==1 to a vector
  for (proc_item in processed_list) {
    # if only not all fullmatches > 0 are == 1
    if (any(proc_item$mask > 0 & proc_item$mask != 1)) {
      # Save the mask to fullmatches
      count <<- count + 1
      fullmatches <<- c(fullmatches, proc_item$mask)
      # false prediction of:
      # 1. masking the >0.5 elements
      # 2. masking the <0.5 elements#
      instance <- proc_item$instance
      sequence_length <- nchar(instance)
      if (sequence_length < 1e6) {
        instance <- paste0(instance, paste(rep("N", 1e6 - sequence_length), collapse = ""))
        proc_item$mask <- c(proc_item$mask, rep(0, 1e6 - sequence_length))
      }
      onehot_instance <- deepG::seq_encoding_label(
        char_sequence = instance,
        maxlen = 1e6,
        start_ind = 1,
        vocabulary = c("A", "C", "G", "T")
      )
      # Calculate the relative frequencies for each channel
      pop.dense <- nt_density(proc_item$instance)
      # Update the one-hot encoded instance
      num_high_spots <- sum(proc_item$mask >= 0.35)
      
      # Sample the same number of spots from proc_item$mask < 0.35
      low_indices <- which(proc_item$mask < 0.35 & proc_item$mask > 0)
      if (length(low_indices) > num_high_spots) {
        sampled_low_indices <- sample(low_indices, num_high_spots)
      } else {
        sampled_low_indices <- low_indices  # If fewer, use all
      }
      random_indices <- sample(which(proc_item$mask > 0), num_high_spots)
      zero_indices <- sample(which(proc_item$mask == 0), num_high_spots)
      
      # Update high_instance
      high_instance <- onehot_instance
      high_instance[, proc_item$mask >= 0.35, ] <- matrix(pop.dense, 
                                                         nrow = sum(proc_item$mask >= 0.35), 
                                                         ncol = length(pop.dense), 
                                                         byrow = TRUE)
      
      # Update low_instance with sampled spots
      low_instance <- onehot_instance
      low_instance[, sampled_low_indices, ] <- matrix(pop.dense, 
                                                      nrow = length(sampled_low_indices), 
                                                      ncol = length(pop.dense), 
                                                      byrow = TRUE)
      random <- onehot_instance
      random[, random_indices, ] <- matrix(pop.dense, 
                                         nrow = length(random_indices), 
                                         ncol = length(pop.dense), 
                                         byrow = TRUE)
      all_low <- onehot_instance
      all_low[, low_indices, ] <- matrix(pop.dense, 
                                        nrow = length(low_indices), 
                                        ncol = length(pop.dense), 
                                        byrow = TRUE)
      allzero <- onehot_instance
      allzero[, zero_indices, ] <- matrix(pop.dense, 
                                         nrow = length(zero_indices), 
                                         ncol = length(pop.dense), 
                                         byrow = TRUE)
      # Predict the spore formation
      high_pred <- predict(spore.model, high_instance, verbose=0)[1,1]
      low_pred <- predict(spore.model, low_instance, verbose=0)[1,1]
      all_low_pred <- predict(spore.model, all_low, verbose=0)[1,1]
      random_pred <- predict(spore.model, random, verbose=0)[1,1]
      allzero_pred <- predict(spore.model, allzero, verbose=0)[1,1]
      predictions <<- rbind(predictions, data.frame(high = high_pred, low = low_pred,
                                                    all_low = all_low_pred, random = random_pred,
                                                    allzero = allzero_pred))
    }
  }
}

# Loop through all file groups
output_dir <- "permutation/tensors/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

walk(files, ~ process_group(.x, output_dir))
count
### proportion of >0 elements
sum(fullmatches > 0) / length(fullmatches)
### among the >0 elements, the proportion of 1
sum(fullmatches[fullmatches > 0] == 1) / sum(fullmatches > 0)
### among the >0 elements, the distribution of the values
### histogramm with relative frequencies
library(ggplot2)

# Filter fullmatches
filtered_matches <- fullmatches[fullmatches > 0]

# Create a data frame for plotting
data <- data.frame(value = filtered_matches)

# Plot histogram with percentages on top of bars
ggplot(data, aes(x = value)) +
  geom_histogram(aes(y = ..count../sum(..count..) * 100), binwidth = 0.05, fill = "white", color = "black") +
  scale_y_continuous(name = "Relative Frequency", labels = scales::percent_format(scale = 1)) +
  labs(x = "Mean Mask Value", title = "Histogram of Mean Mask Values", subtitle = ">0: 13.26%, in 56 Species") +
  geom_text(
    stat = "bin",
    aes(
      label = ifelse(..count.. > 0, scales::percent(..count../sum(..count..), accuracy = 0.1), ""),
      y = ..count../sum(..count..) * 100
    ),
    binwidth = 0.05, vjust = -0.5, size = 2
  ) +
  theme_bw()

#### boxplot for all three cols in predictions
library(tidyverse)
predictions_long <- pivot_longer(predictions, cols = everything(), names_to = "Category", values_to = "Value")

library(car)
shapiro.test(predictions$high)
shapiro.test(predictions$low)
leveneTest(c(predictions$high, predictions$low), 
           group = rep(c("high", "low"), each = nrow(predictions)))
pairwise.wilcox.test(
  predictions_long$Value,
  predictions_long$Category,
  p.adjust.method = "bonferroni"
)

# Create boxplot
library(ggsignif)

# Boxplot with significance annotations using ggsignif
ggplot(predictions_long %>%
         mutate(Category = factor(Category, 
                                  levels = c("high", "low", "all_low", "random", "allzero"))),
       aes(x = Category, y = Value)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Masking criteria",
    y = "False Prediction"
  ) +
  scale_x_discrete(labels = c(
    "high" = "≥ 0.35", 
    "low" = "sampled < 0.35", 
    "all_low" = "all < 0.35", 
    "random" = "sampled > 0",
    "allzero" = "sampled = 0"
  )) +
  geom_signif(
    comparisons = list(c("high", "low"), c("high", "all_low"), c("low", "random"),
                       c("random", "allzero"), c("low", "allzero")),
    annotations = c("***", "ns","***", "***", "***"),  # Use significance levels or exact p-values
    y_position = c(0.7, 0.8, 1.0, 0.65, 0.55) * max(predictions_long$Value),
    tip_length = 0.03,
    color = "red"
  )

