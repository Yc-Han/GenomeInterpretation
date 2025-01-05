### GENTOMA Usage ###

## Sporulation Example
library(reticulate)

#### Set up the environment ####
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "caret", "keras", "magrittr", "patchwork", 
                       "ggseqlogo", "openxlsx", "zoo",
                       "future", "future.apply", "rtracklayer")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))
source("permutation/script/gentoma.R")

#### Load the data ####
maxlen <- 1000000
target_from_csv <- "sporulation/sporeinfo.csv"
target_df <- read.csv(target_from_csv)
label_names <- names(target_df)[names(target_df) != "file"]
print(label_names)
folder_path <- "sporulation/genome"
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)

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
pred <- predict(model, onehot_instance)
onehot_baseline <- array(rep(nt_density(instance_sub), each = dim(onehot_instance)[2]), 
                         dim = dim(onehot_instance))
# Set parameters
set.seed(123) # 123, 42, 777, 1200, 9999
window.size <- 0.10
total.sub <- round(window.size * maxlen)
seg.len <- round(0.01 * maxlen)
seg.num <- total.sub / seg.len
iterations <- 200
fidx <- 1
N <- 20    # Max generations
M <- 3     # Stagnant generations
lambda <- 1

# GOM step
result <- genetically_optimized_masking(instance, onehot_instance, pop.dense, maxlen, seg.len,
                                        seg.num, iterations, model, fidx, N=N, M=M,
                                        lambda=lambda)
individual <- result$best_individual

source("ig_modified.R")
gff_path <- "sporulation/gff/gff_new"
gff_path_1 <- paste0(gff_path, "/GCF_029537415.1_ASM2953741v1_genomic.gff")
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
gff.df_fil <- gff.df[grep("spore|sporulation", gff.df$product, ignore.case = TRUE), ]
gff.df_fil <- gff.df_fil %>%
  mutate(length = end - start)
igs <- ig_modified(
  input_seq = onehot_instance,
  baseline_type = "modify",
  baseline_onehot = onehot_baseline,
  target_class_idx = 2,
  model = model,
  num_baseline_repeats = 1,
  pred_stepwise = TRUE
)
abs_sum <- rowSums(abs(as.array(igs)))
sum <- rowSums(as.array(igs))
df <- data.frame(sum = sum, abs_sum = abs_sum, position = 1:maxlen)
df$smoothed_sum <- ksmooth(df$position, df$sum, kernel = "box", bandwidth = 100)$y
df$smoothed_abs_sum <- ksmooth(df$position, df$abs_sum, kernel = "box", bandwidth = 100)$y
df_zoo <- zoo(df$abs_sum, order.by = df$position)
roll_mean <- rollmean(df_zoo, k = 100, align = "left", fill = NA)

png("permutation/plots/Bacillus_Subtilis_zoo.png", width = 3800, height = 2400, res = 600)
plot(roll_mean, col = "red", lwd = 0.1, xlab = "Position (bp)", ylab = "Mean abs sum of IG / 100 bp",
     main = "Bacillus subtilis, Baseline: nt_density")
current_ymin <- par("usr")[3]
current_ymax <- par("usr")[4]
for(i in 1:nrow(gff.df_fil)) {
  rect(gff.df_fil$start[i], current_ymin, gff.df_fil$end[i], current_ymax, col = rgb(0, 0, 1, 0.5), border = NA)
}
legend("topright", legend = c("Spore-related genes"),
       fill = rgb(0, 0, 1, 0.5), border = NA, cex = 0.75)
dev.off()

png("permutation/plots/Bacillus_Subtilis_Coverage_Plot_seed9999.png", width = 3800, height = 2400, res = 600)
plot(roll_mean, col = "grey", lwd = 0.1, xlab = "Position (bp)", ylab = "Mean sum of IG / 100 bp")
current_ymin <- par("usr")[3]
current_ymax <- par("usr")[4]
for (i in 1:length(individual$starts)) {
  # adjust color to orange
  rect(individual$starts[i], current_ymin, individual$ends[i], current_ymax,
       col = rgb(1, 0, 0, 0.5), border = NA)
}
for(i in 1:nrow(gff.df_fil)) {
  rect(gff.df_fil$start[i], current_ymin, gff.df_fil$end[i], current_ymax,
       col = rgb(0, 0, 1, 0.5), border = NA)
}
# Highlight overlaps with green borders and add 'x' marks
for (i in 1:nrow(gff.df_fil)) {
  for (j in 1:length(individual$starts)) {
    overlap_start <- max(gff.df_fil$start[i], individual$starts[j])
    overlap_end <- min(gff.df_fil$end[i], individual$ends[j])
    if (overlap_start < overlap_end) {  # Overlap exists
      rect(overlap_start, current_ymin, overlap_end, current_ymax, 
           col = NA, border = "green", lwd = 1)  # Green border with bold width
      
      # Add 'x' mark at the top of the y-axis for overlap
      points(x = mean(c(overlap_start, overlap_end)), y = 0, hjust = 1,
             pch = 4, col = "darkgreen", cex = 0.7, lwd = 1)  # 'x' mark
    }
  }
}
# add legend: blue: spore-related genes, orange: GOM selected regions
#legend("topright", legend = c("Masked regions", "Spore-related genes", "Overlap"),
#       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), "green"), border = NA, cex = 0.75)
dev.off()

# excel tabel for matched genes
setDT(gff.df)
setDT(df)
setDT(gff.df_fil)

# Set key for faster subsetting
setkey(df, position)

# Function to get mean of 'sum' based on start and end positions
get_mean_sum <- function(start, end) {
  mean(df[position >= start & position <= end, abs_sum])
}

# Apply the function to each row in gff.df and create a new column
gff.df[, mean_sum := get_mean_sum(start, end), by = 1:nrow(gff.df)]
gff.df_fil[, mean_sum := get_mean_sum(start, end), by = 1:nrow(gff.df_fil)]
ecdf_sum <- ecdf(df$abs_sum)

# Apply ECDF to each value in mean_sum to create a new column quantile_position
gff.df[, quantile_position := ecdf_sum(mean_sum)]
gff.df_fil[, quantile_position := ecdf_sum(mean_sum)]
write.xlsx(gff.df, "permutation/plots/bs_importance.xlsx")

#### Multiple FASTAs ####
# Set parameters for multiple FASTAs
window.size <- 0.25
total.sub <- round(window.size * maxlen)
seg.len <- 1000  # Segments of 1000 bp
seg.num <- total.sub / seg.len
iterations <- 150
fidx <- 1
N <- 10    # Max generations
M <- 3     # Stagnant generations
sample.count <- 10
multiple_fasta_results <- gom_all_fasta(
  target_df = target_df,
  folder_path = folder_path,
  maxlen = maxlen,
  seg.len = seg.len,
  seg.num = seg.num,
  iterations = iterations,
  model = model,
  fidx = fidx,
  N = N,
  M = M,
  lambda = 1,             # Penalty factor for coverage
  update_rate = 0.75,    # 75% of population is replaced each iteration
  mutation_prob = 0.2,   # Mutation probability per segment
  sample.count = sample.count
)
