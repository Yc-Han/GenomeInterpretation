### Phyla-stratified GENTOMA

library(reticulate)
#### Set up the environment ####
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "caret", "keras", "magrittr", "openxlsx", "parallel", "rtracklayer")

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))
# Source necessary scripts
source("permutation/script/gentoma.R")
source("permutation/script/fixed_encoder.R")
# Get the file name from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("No arguments supplied.")
}

fasta_file <- args[1]

# Load the model
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)

# Set parameters
set.seed(123)  # For reproducibility
num_seeds <- 10
maxlen <- 1000000
window.size <- 0.15
total.sub <- round(window.size * maxlen)
seg.len <- 5000
seg.num <- total.sub / seg.len
iterations <- 300
fidx <- 1
N <- 20    # Max generations
M <- 3     # Stagnant generations
lambda <- 1
update_rate <- 0.75
mutation_prob <- 0.2

# Print setup information
print_setup(maxlen, seg.len, seg.num, iterations, fidx, N, M,
            update_rate, mutation_prob, pop.dense = NA, lambda, pred.orignal = NA)

# Process the file
seeds <- sample(1e6, num_seeds)
cat("Processing file:", fasta_file, "\n")

# Initialize results list for this file
results <- list()

for (seed in seeds) {
  cat("Processing file:", fasta_file, "with seed:", seed, "\n")
  # Set the seed
  set.seed(seed)
  result <- tryCatch({
    gom_fasta(fasta_file, maxlen, seg.len, seg.num, iterations, model,
              fidx, N, M, lambda, update_rate, mutation_prob)
  }, error = function(e) {
    cat("Error processing file:", fasta_file, "with seed:", seed, "\n", e$message, "\n")
    NULL
  })
  output_filename <- paste0("permutation/singles_longer/results_", basename(fasta_file), "_", seed, ".rds")
  saveRDS(result, output_filename)
}

