library(deepG)
set.seed(123)
# create data
seq_length <- 100000
maxlen <- 600
voc <- c('A', 'C', 'G', 'T')
s1 <- deepG:::random_seq(vocabulary = voc, seq_length = seq_length)
s2 <- deepG:::random_seq(vocabulary = voc, seq_length = seq_length)
s1_val <- deepG:::random_seq(vocabulary = voc, seq_length = seq_length)
s2_val <- deepG:::random_seq(vocabulary = voc, seq_length = seq_length)

# add pattern to s1
pattern <- 'GGTCGACACGAGTATAGCTAGTATACTCCG'
freq <- 100
starts <- sample(1:(seq_length - nchar(pattern)), floor(seq_length)/freq)
starts_val <- sample(1:(seq_length - nchar(pattern)), floor(seq_length)/freq)

# Define the function to insert pattern into sequence
insert_pattern <- function(sequence, pattern, starts) {
  # Loop through each start position and insert the pattern
  for (start in starts) {
    substr(sequence, start, start + nchar(pattern) - 1) <- pattern
  }
  return(sequence)
}

# Insert pattern into s1 at the start positions
s1 <- insert_pattern(s1, pattern, starts)
s1_val <- insert_pattern(s1_val, pattern, starts_val)

find_pattern_positions <- function(sequence, pattern) {
  
  positions <- list()
  count <- 1
  pattern_length <- nchar(pattern)
  
  # Loop through the sequence to find matches
  for (i in 1:(nchar(sequence) - pattern_length + 1)) {
    substring <- substr(sequence, i, i + pattern_length - 1)
    
    # Check if the substring matches the pattern
    if (substring == pattern) {
      positions[[count]] <- i
      count <- count + 1
    }
  }
  
  return(unlist(positions))
}

find_pattern_positions(s1, pattern)
find_pattern_positions(s2, pattern)
find_pattern_positions(s1_val, pattern)
find_pattern_positions(s2_val, pattern)

microseq::writeFasta(data.frame(Sequence = s1, Header = 's1'), 'ig_data/s1.fasta')
microseq::writeFasta(data.frame(Sequence = s2, Header = 's2'), 'ig_data/s2.fasta')
microseq::writeFasta(data.frame(Sequence = s1_val, Header = 's1'), 'ig_data/s1_val.fasta')
microseq::writeFasta(data.frame(Sequence = s2_val, Header = 's2'), 'ig_data/s2_val.fasta')

# 

model <- create_model_lstm_cnn(
  maxlen = maxlen,
  layer_lstm = NULL,
  layer_dense = c(2L),
  vocabulary_size = 4,
  kernel_size = c(5, 7, 9),  # Smaller and varied filter sizes
  filters = c(32, 64, 128),
  pool_size = c(3, 3, 3),    # Smaller pooling sizes
  learning_rate = 0.0005     # Slightly increased learning rate
)

train_model(
  model = model,
  path = c('ig_data/s1.fasta', 'ig_data/s2.fasta'),
  path_val = c('ig_data/s1_val.fasta', 'ig_data/s2_val.fasta'),
  vocabulary_label = c("s1", "s2"),
  batch_size = 64,
  steps_per_epoch = 300,
  epochs = 10,
  step = 250,
  seed = sample(1:1000, 2)
)

#random_seq <- sample(voc, maxlen, replace = TRUE) %>% paste(collapse = '')
#random_seq <- strrep('T', maxlen)

#start <- 10
#m_steps <- 50
#baseline_type <- 'zero' # 'shuffle' # 'zero'
#num_baseline_repeats <- ifelse(baseline_type == 'shuffle', 30, 1)
#substr(random_seq, start, start + nchar(pattern) - 1) <- pattern

#random_seq <- SeqSyn(maxlen, 1, codon.dict, by.codon=FALSE)
#random_seq <- insert_motif(random_seq, pattern)
set.seed(42)
random_seq <- sample(voc, maxlen, replace = TRUE) %>% paste(collapse = '')
starts <- c(10, 110, 220, 300, 400, 500)
for (start in starts) {
  substr(random_seq, start, start + nchar(pattern) - 1) <- pattern
  input_seq <- seq_encoding_label(char_sequence = random_seq, start_ind = 1, maxlen = maxlen, vocabulary = voc)
}
#onehot_instance <- seq_encoding_label(char_sequence = input_seq[1,,], start_ind = 1, maxlen = maxlen, vocabulary = voc)
onehot_instance <- input_seq
onehot_baseline_25 <- onehot_instance*0 + 0.25
predict(model, onehot_instance, verbose = 0)
ig25 <- ig_modified(m_steps = 100,
                    input_seq = onehot_instance,
                    baseline_type = "modify",
                    baseline_onehot = onehot_baseline_25,
                    target_class_idx = 1,
                    model = model,
                    num_baseline_repeats = 1)
motif_pos <- find_motif(random_seq, pattern)
motif_end <- motif_pos + nchar(pattern)
sum <- rowSums(as.array(ig25))
abs_sum <- rowSums(abs(as.array(ig25)))
df25 <- data.frame(abs_sum = abs_sum, sum=sum, position = 1 : maxlen)
pattern_length<- nchar(pattern)
ggplot(df25, aes(x = position, y = sum))+
  geom_rect(aes(xmin = 10, xmax = 10 + pattern_length, ymin = -Inf, ymax = Inf), 
            fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = 110, xmax = 110 + pattern_length, ymin = -Inf, ymax = Inf), 
            fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = 220, xmax = 220 + pattern_length, ymin = -Inf, ymax = Inf), 
            fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = 300, xmax = 300 + pattern_length, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = 400, xmax = 400 + pattern_length, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2) +
  geom_rect(aes(xmin = 500, xmax = 500 + pattern_length, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2) +
  geom_point() + theme_bw() + labs(subtitle = "Baseline 0.25, Pred ~1")

igmat <- as.data.frame(t(as.matrix(ig25)))
rownames(igmat) <- c("A", "C", "G", "T")
igmat <- as.matrix(igmat)
ggseqlogo(igmat, method='custom', seq_type='dna') + xlim(5,155) + labs(x="bp", y="IG") +  #geom_rect(aes(xmin = motif_pos, xmax = motif_end, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.2) +
  # add motif as text
  geom_text(aes(x = motif_pos + 6, y = 0.02, label = pattern), size = 3.5, color = "blue")

motifregion <- igmat[, motif_pos:(motif_end -1)]
# get row name of max in each column
row_max_in_cols <- apply(motifregion, 2, function(col) rownames(igmat)[which.max(col)])

# Print the result
interpret <- unname(row_max_in_cols)

motif_letters <- strsplit(pattern, split = "")[[1]]
sum(interpret == motif_letters) / length(motif_letters)

