#### SETUP ####

# install.packages("tensorflow")
# tensorflow::install_tensorflow()
# devtools::install_github("GenomeNet/deepG")


# install.packages("BiocManager")
# BiocManager::install("rtracklayer")

required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "caret", "keras", "magrittr", "patchwork", 
                       "ggseqlogo", "openxlsx", "zoo",
                       "future", "future.apply", "rtracklayer")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)

# Load all the packages
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)

# source("genepermutation.R")
# source("ig_modified.R")
# source("seqsyn.R")

#### BASE VARS ####

maxlen <- 1000000
class_weight <- list("0" = 1, "1" = 3.66)
target_from_csv <- "sporulation/sporeinfo.csv"
target_df <- read.csv(target_from_csv)
label_names <- names(target_df)[names(target_df) != "file"]
print(label_names)
data_path <- "sporulation/genome"
gff_path <- "sporulation/gff/gff_new"
fasta_path <- "sporulation/genome/GCF_029537415.1_ASM2953741v1_genomic.fasta"
# scale-up: for all file with target_df$ability_TRUE == 1
file_label <- target_df %>% dplyr::filter(file == basename(fasta_path))
print(file_label)
gff_path_1 <- paste0(gff_path, "/GCF_029537415.1_ASM2953741v1_genomic.gff")
# scale-up: for all file with target_df$ability_TRUE ==1, change .fasta to .gff
gff <- readGFF(gff_path_1, version=0,
               columns=NULL, tags=NULL, filter=NULL, nrows=-1,
               raw_data=FALSE)
list_data <- gff@listData
gff.df <- data.frame(
  start = list_data$start,
  end = list_data$end,
  Name = list_data$Name,
  gene = list_data$gene,
  product = list_data$product
) %>%
  filter(end < maxlen) %>%
  filter(!is.na(product))
gff.df<- gff.df %>%
  mutate(length = end - start)
# name hypothetical genes
gff.df <- gff.df %>%
  arrange(start) %>% 
  mutate(gene = ifelse(is.na(gene), paste0("hypo_", row_number()), gene))
spore.model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
instance <- microseq::readFasta(fasta_path)$Sequence[1]
instance_sub <- substr(instance, 1, maxlen) # fehlende Länge: N ausfüllen
sequence_length <- nchar(instance)
if (sequence_length < maxlen) {
  instance_sub <- paste0(instance_sub, paste(rep("N", maxlen - sequence_length), collapse = ""))
}
onehot_instance <-  seq_encoding_label(char_sequence = instance_sub,
                                       maxlen = maxlen,
                                       start_ind = 1,
                                       vocabulary = c("A", "C", "G", "T"))
pred.ori <- predict(spore.model, onehot_instance, verbose = 0)[1,1]

# how many steps are needed?
# Function to calculate the probability or minimum iterations (steps)
calculate_sampling <- function(q_single_sample, step = NULL, prob = NULL) {
  # If step is provided, calculate the probability
  if (!is.null(step)) {
    p_not_included_one_sample <- 1 - q_single_sample
    prob_included <- 1 - (p_not_included_one_sample ^ step)
    return(prob_included)
  }
  
  # If probability is provided, calculate the minimum steps required
  if (!is.null(prob)) {
    p_not_included_one_sample <- 1 - q_single_sample
    min_steps <- log(1 - prob) / log(p_not_included_one_sample)
    return(ceiling(min_steps))
  }
  
  # If neither step nor prob is provided, return an error message
  stop("Either 'step' or 'prob' must be provided.")
}

# calculate_sampling(q_single_sample = 0.05, step = 1000)

#### GIVING WINDOWS ####

# parallel setup
set.seed(42)
plan(multicore, workers = 15)
window.size <- 0.40
total.sub <- round(window.size * maxlen)
seg.len <- 1000 #round(0.001 * maxlen)
seg.num <- total.sub / seg.len
iterations <- 100
run_iteration <- function(i) {
  # Randomly sample seg.num start positions
  samps <- sample(1:(maxlen - seg.len), seg.num, replace = FALSE)
  
  # Create a default mask of zeros (mask will be same size as onehot_instance)
  mask <- array(0, dim = c(1, maxlen, 4))  # Default all zeros
  
  # Set the selected segments to 1 (keeping original information)
  for (samp in samps) {
    mask[1, samp:(samp + seg.len), ] <- 1  # Set the region to 1 (keep)
  }
  
  # Apply the mask to the one-hot instance
  new_instance <- onehot_instance * mask
  
  # Perform the prediction using the model
  pred <- spore.model %>% predict(new_instance, verbose = 0)
  
  # Return a list: samps (start positions) and pred[1, 2]
  if (length(pred) > 1) {
    return(list(samps = samps, pred = pred[1, 2]))
  } else {
    return(NULL)  # Handle cases where the prediction is invalid
  }
}

results <- future_lapply(1:iterations, run_iteration, future.seed = TRUE)

# Extract predictions from results
predictions <- unlist(lapply(results, function(res) if (!is.null(res)) res$pred))

# Sanity check: Compute the average of the highest 10% of predictions
if (length(predictions) > 0) {
  top_10_percent <- sort(predictions, decreasing = TRUE)[1:ceiling(0.1 * length(predictions))]
  avg_top_10_percent <- mean(top_10_percent)
  
  # Check if the average of the top 10% is less than 0.75 and give a warning if true
  if (avg_top_10_percent < 0.75) {
    warning("The average of the highest 10% of predictions is less than 0.75!")
  }
  
  # Print the average for debugging
  cat("Average of top 10% predictions:", avg_top_10_percent, "\n")
} else {
  warning("No valid predictions were made.")
}
max(predictions)
best_samps <- results[[which.max(predictions)]]$samps

#### TRIMMING WINDOWS ####

# for the best performing windowed sequence, we iteratively delete more q% of its samples
# for each q, it is repeated 100 times
# the maximum of the predictions is calculated
# stop if the max pred drops under 0.8
# report everything
set.seed(42)
interatively_remove <- function(q) {
  # Remove q% of the best_samps randomly
  num_to_keep <- round((1 - q) * length(best_samps))  # Number of samples to keep
  samps <- sample(best_samps, num_to_keep, replace = FALSE)
  mask <- array(0, dim = c(1, maxlen, 4))  # Default all zeros
  for (samp in samps) {
    mask[1, samp:(samp + seg.len), ] <- 1  # Set the region to 1 (keep)
  }
  new_instance <- onehot_instance * mask
  pred <- spore.model %>% predict(new_instance, verbose = 0)
  if (length(pred) > 1) {
    return(list(samps = samps, pred = pred[1, 2]))  # Return the correct prediction
  } else {
    return(NULL)
  }
}

# Initialize the iterative removal process
q <- 0  # Initial percentage removal
max_pred <- Inf
list_results <- list()

while (max_pred > 0.75) {
  iterations <- 100  # Number of iterations for each q value
  q <- q + 0.05  # Increment q by 1%
  print(paste0("Removing q: ", q*100, "%"))
  
  # Perform the iterative removal in parallel
  results <- future_lapply(1:iterations, function(x) interatively_remove(q), future.seed = TRUE)
  
  # Extract predictions from the results
  predictions <- unlist(lapply(results, function(res) if (!is.null(res)) res$pred))
  samps <- unlist(lapply(results, function(res) if (!is.null(res)) res$samps))
  
  # Get the maximum prediction for this q
  max_pred <- max(predictions, na.rm = TRUE)
  # Store the results for this q
  list_results[[as.character(q)]] <- list(predictions = predictions, max_pred = max_pred, samps = samps)
  
  print(paste("Max pred:", max_pred))
}

# extract predictions and plot boxplot, grouped by q
predictions <- unlist(lapply(list_results, function(res) res$pred))
q_values <- as.numeric(names(list_results))
q_values <- q_values[order(q_values)]
q_values <- round(q_values, 2)
q_values <- q_values[order(q_values)]
q_values <- as.character(q_values)
# repeate each 50
q_values <- rep(q_values, each = 100)
boxplot(predictions ~ q_values, xlab = "q", ylab = "Prediction", main = "Prediction vs. q")
# add hline for pred = 0.75
abline(h = 0.75, col = "red", lty = 2)

# extract best_samp as q=10, pred=max
predictions <- list_results[["0.15"]]$predictions
num_to_keep <- round((1 - 0.15) * length(best_samps))
index <- which(predictions == max(predictions))
best_samps <- list_results[["0.15"]]$samps[num_to_keep * (index - 1) + 1:num_to_keep*index]

chunk_size <- 10000  # Number of positions per line
num_chunks <- ceiling(maxlen / chunk_size)  # Total number of chunks to plot

# Setup an empty plot with enough space for each chunk
par(mar = c(2, 1, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(1, chunk_size), ylim = c(0, num_chunks), axes = FALSE, main = "")

# Add x-axis only on the first chunk
axis(1, at = seq(1, maxlen, by = 1000), labels = c(1, seq(1000, maxlen-1, by = 1000)))

# Loop over each chunk and plot the sequence segments
for (i in 1:num_chunks) {
  # Calculate the start and end positions of the current chunk
  chunk_start <- (i - 1) * chunk_size + 1
  chunk_end <- min(i * chunk_size, maxlen)
  # Plot the black baseline for the current chunk
  segments(1, num_chunks - i + 1, chunk_end - chunk_start + 1, num_chunks - i + 1, col = "black", lwd = 1)
  
  # Plot the red segments corresponding to the "given" parts
  for (j in best_samps) {
    if (j >= chunk_start && j <= chunk_end) {
      # Calculate the start and end of the red segment relative to the chunk
      seg_start <- j - chunk_start + 1
      seg_end <- min(j + seg.len - chunk_start + 1, chunk_size)
      segments(seg_start, num_chunks - i + 1, seg_end, num_chunks - i + 1, col = "yellow", lwd = 1)
    }
  }
}
#### INFERENCE ####
windows <- data.table(samp = best_samps, samp_end = best_samps + seg.len)
saveRDS(windows, "permutation/windows.rds")
setDT(gff.df)
setkey(gff.df, start, end)
setkey(windows, samp, samp_end)
overlaps <- foverlaps(windows, gff.df, by.x = c("samp", "samp_end"), by.y = c("start", "end"), nomatch = NA)
saveRDS(overlaps, "permutation/windows_overlaps.rds")
# get difference in unique(overlaps$gene) and unique(gff.df$gene)
unique_genes <- unique(overlaps$gene)
unique_genes_gff <- unique(gff.df$gene)
unique_genes_not_in_gff <- setdiff(unique_genes_gff, unique_genes)
unique_genes_not_in_gff
# only 271/1135 (24%) genes are included, coverage=19.6%

# extract samp:samp_end in instance_sub (a string), then ggseqlogo
# extract the sequence for each samp:samp_end as text
extract_sequence <- function(instance_sub, samp, samp_end) {
  return(substr(instance_sub, samp, samp_end))  # Extracts sequence between samp and samp_end
}
windows[, sequence := mapply(extract_sequence, instance_sub, samp, samp_end)]
# Plot the sequence logo
ggseqlogo(windows$sequence) +
  # show x ticks at every 10th position
  scale_x_continuous(breaks = c(1, seq(100, seg.len-1, by = 100)))
# concatenate the sequences
synthetic <- paste(windows$sequence, collapse = "")
nucleotide_percentages <- prop.table(table(strsplit(synthetic, "")[[1]])) * 100
print(nucleotide_percentages)
# repeat the synthetic until maxlen
repeated_synthetic <- strrep(synthetic, ceiling(maxlen / nchar(synthetic)))
repeated_synthetic <- substr(repeated_synthetic, 1, maxlen)
synthetic_onehot <- seq_encoding_label(char_sequence = repeated_synthetic,
                                       maxlen = maxlen,
                                       start_ind = 1,
                                       vocabulary = c("A", "C", "G", "T"))
predict(spore.model, synthetic_onehot, verbose = 0)

# create a TTTTT and a AAAA sequence for maxlen each prediction
# TTTTT
t_seq <- strrep("T", maxlen)
t_onehot <- seq_encoding_label(char_sequence = t_seq,
                              maxlen = maxlen,
                              start_ind = 1,
                              vocabulary = c("A", "C", "G", "T"))
predict(spore.model, t_onehot, verbose = 0)

# AAAA
a_seq <- strrep("A", maxlen)
a_onehot <- seq_encoding_label(char_sequence = a_seq,
                              maxlen = maxlen,
                              start_ind = 1,
                              vocabulary = c("A", "C", "G", "T"))
predict(spore.model, a_onehot, verbose = 0) # pred = 1!

# AATT...
an_seq <- strrep("AATT", maxlen / 4)
an_onehot <- seq_encoding_label(char_sequence = an_seq,
                                  maxlen = maxlen,
                                  start_ind = 1,
                                  vocabulary = c("A", "C", "G", "T"))
predict(spore.model, an_onehot, verbose = 0) # pred = 1!

plan(sequential)