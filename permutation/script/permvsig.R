required_packages <- c("deepG", "ggplot2", "microseq", "data.table", "seqinr", 
                       "dplyr", "caret", "keras", "magrittr", "patchwork", 
                       "ggseqlogo", "openxlsx", "zoo", 
                       "reticulate", "future", "future.apply", "rtracklayer")
invisible(lapply(required_packages, library, character.only = TRUE))
source("ig_modified.R")
set.seed(42)
maxlen <- 10000
# load model
checkpoint_path <- file.path("permutation/checkpoints")
run_name <- "permsyn"
model <- load_cp(paste0(checkpoint_path, "/", run_name), 
                 cp_filter = "last_ep")

# load data: target_val.rds
target.val <- readRDS("permutation/syndata/target_val.rds")
gene1 <- readRDS("permutation/syndata/gene1.rds")
gene2 <- readRDS("permutation/syndata/gene2.rds")
instance <- substr(target.val, 1, maxlen)
onehot_instance <-  seq_encoding_label(char_sequence = instance,
                                       maxlen = maxlen,
                                       start_ind = 1,
                                       vocabulary = c("A", "C", "G", "T"))
find_genes <- function(sequence, gene1, gene2) {
  # Use fixed string matching for efficiency
  pos1 <- gregexpr(gene1, sequence, fixed = TRUE)[[1]][1]
  pos2 <- gregexpr(gene2, sequence, fixed = TRUE)[[1]][1]
  
  # Handle cases where the gene is not found
  if (pos1 == -1) pos1 <- NA
  if (pos2 == -1) pos2 <- NA
  
  # Return a named vector with starting positions
  return(c(gene1 = pos1, gene2 = pos2))
}
motif_pos <- find_genes(instance, gene1, gene2)
motif_end <- c(motif_pos[1] + nchar(gene1) - 1, motif_pos[2] + nchar(gene2) - 1)

onehot_baseline_25 <- onehot_instance * 0 + 0.25
pred <- predict(model, onehot_instance, verbose = 0)
pred

null_instance <- onehot_instance * 0
null_pred <- predict(model, null_instance, verbose = 0)
null_pred

#### Modified IG
ig25 <- ig_modified(m_steps = 400,
                    input_seq = onehot_instance,
                    baseline_type = "modify",
                    baseline_onehot = onehot_baseline_25,
                    target_class_idx = 2,
                    model = model,
                    num_baseline_repeats = 1,
                    pred_stepwise = TRUE)
abs_sum <- rowSums(abs(as.array(ig25)))
df25 <- data.frame(abs_sum = abs_sum, position = 1 : maxlen)

df_zoo <- zoo(df25$abs_sum, order.by = df25$position)
roll_mean <- rollmean(df_zoo, k = 50, fill = NA)
png("permutation/plots/syn_ig.png", width = 8, height = 6, units = "in", res = 300)
plot(roll_mean, col = "red", xlab = "Position (bp)", ylab = "Mean abs sum of IG / 50 bp")
current_ymin <- par("usr")[3]
current_ymax <- par("usr")[4]
for(i in 1:length(motif_pos)) {
  rect(motif_pos[i], current_ymin, motif_end[i], current_ymax, col = rgb(0, 0, 1, 0.2), border = NA)
}
# Save the plot
dev.off()

##### Windowing Algorithm #####
source("permutation/script/window_helpers.R")
set.seed(56)
window.size <- 0.25 # start with 40% of the sequence
total.sub <- round(window.size * maxlen)
seg.len <- 1000 # 1000 bp ~ avg. length of a gene
seg.num <- total.sub / seg.len
iterations <- 100
pop.dense <- nt_density(instance)
results <- future_lapply(1:iterations, give_window, future.seed = TRUE)
predictions <- unlist(lapply(results, function(res) if (!is.null(res)) res$pred))
# Sanity check with skipping condition
max_prediction <- max(predictions)
cat("Max false prediction:", max_prediction, "\n")

best_samps <- results[[which.max(predictions)]]$samps
#### SKIPPED TRIMMING WINDOWS
chunk_size <- 2000  # Number of positions per line
num_chunks <- ceiling(maxlen / chunk_size)  # Total number of chunks to plot

# Setup an empty plot with enough space for each chunk
png("permutation/plots/syn_remove_window.png", width = 8, height = 6, units = "in", res = 300)
par(mar = c(2, 1, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(1, chunk_size), ylim = c(0, num_chunks), axes = FALSE, main = "")
# Add x-axis only on the first chunk
axis(1, at = seq(1, maxlen, by = 100), labels = c(1, seq(100, maxlen-1, by = 100)))
# Loop over each chunk and plot the sequence segments
for (i in 1:num_chunks) {
  # Calculate the start and end positions of the current chunk
  chunk_start <- (i - 1) * chunk_size + 1
  chunk_end <- min(i * chunk_size, maxlen)
  
  # Plot the black baseline for the current chunk
  segments(1, num_chunks - i + 1, chunk_end - chunk_start + 1, num_chunks - i + 1, col = "darkgrey", lwd = 5)
  
  # Plot the blue regions corresponding to the "motif" positions
  for (j in seq_along(motif_pos)) {
    if (motif_pos[j] <= chunk_end && motif_end[j] >= chunk_start) {
      # Calculate the start and end of the blue region relative to the chunk
      blue_start <- max(motif_pos[j], chunk_start) - chunk_start + 1
      blue_end <- min(motif_end[j], chunk_end) - chunk_start + 1
      rect(blue_start, num_chunks - i + 0.8, blue_end, num_chunks - i + 1.2, col = rgb(0, 0, 1, 0.3), border = NA)
      
      # Handle any overflow into the next chunk(s)
      if (motif_end[j] > chunk_end) {
        remaining_len <- motif_end[j] - chunk_end
        next_chunk <- i + 1
        cat(j, remaining_len, next_chunk, "\n")
        while (remaining_len > 0 && next_chunk <= num_chunks) {
          next_chunk_start <- (next_chunk - 1) * chunk_size + 1
          next_chunk_end <- min(next_chunk * chunk_size, maxlen)
          next_blue_start <- 1
          next_blue_end <- min(remaining_len, next_chunk_end - next_chunk_start + 1)
          cat(next_blue_start, next_blue_end, "\n")
          rect(next_blue_start, num_chunks - next_chunk + 0.8, next_blue_end, num_chunks - next_chunk + 1.2, col = rgb(0, 0, 1, 0.1), border = NA)
          
          remaining_len <- remaining_len - (next_blue_end - next_blue_start + 1)
          next_chunk <- next_chunk + 1
        }
      }
    }
  }
  
  # Plot the red/yellow segments corresponding to the "best_samps" parts
  for (j in best_samps) {
    if (j + seg.len >= chunk_start && j <= chunk_end) {
      # Part of the segment is within this chunk
      seg_start <- max(j, chunk_start) - chunk_start + 1
      seg_end <- min(j + seg.len - 1, chunk_end) - chunk_start + 1
      segments(seg_start, num_chunks - i + 1, seg_end, num_chunks - i + 1, col = "red", lwd = 5)
      
      # If the segment extends beyond the current chunk, plot the overflow in the next chunk(s)
      if (j + seg.len > chunk_end) {
        # Calculate how much of the segment spills over into the next chunk(s)
        remaining_len <- j + seg.len - chunk_end
        
        # Loop through the remaining chunks and plot the remaining parts
        next_chunk <- i + 1
        while (remaining_len > 0 && next_chunk <= num_chunks) {
          next_chunk_start <- (next_chunk - 1) * chunk_size + 1
          next_chunk_end <- min(next_chunk * chunk_size, maxlen)
          next_seg_start <- 1
          next_seg_end <- min(remaining_len, next_chunk_end - next_chunk_start + 1)
          segments(next_seg_start, num_chunks - next_chunk + 1, next_seg_end, num_chunks - next_chunk + 1, col = "red", lwd = 5)
          
          # Update remaining length and move to the next chunk
          remaining_len <- remaining_len - (next_seg_end - next_seg_start + 1)
          next_chunk <- next_chunk + 1
        }
      }
    }
  }
}
dev.off()