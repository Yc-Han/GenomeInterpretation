### Phyla-stratified GENTOMA
library(reticulate)
#### Set up the environment ####
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "caret", "keras", "magrittr","openxlsx", "parallel", "rtracklayer")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

microbes <- read.xlsx("sporulation/microbe.cards table S1.xlsx")
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
microbes <- microbes %>%
  mutate(Gff.file = gsub(".fasta", ".gff", Fasta.file)) %>%
  filter(Member.of.WA.subset == TRUE & !is.na(Fasta.file) & !is.na(Gff.file)) %>%
  filter(Spore.formation == TRUE)

Order_counts <- table(microbes$Phylum)
phyla_to_keep <- names(Order_counts)
# Sample up to 5 genomes per Order
set.seed(123)  # For reproducibility
sampled_files <- list()
for (phy in phyla_to_keep) {
  phy_data <- microbes %>%
    filter(Phylum == phy)
  num_genomes <- nrow(phy_data)
  sampled_indices <- if (num_genomes >= 36) sample(1:num_genomes, 36) else 1:num_genomes
  sampled_files[[phy]] <- paste0("sporulation/genome/", phy_data$Fasta.file[sampled_indices])
}
file.names <- unlist(sampled_files)
writeLines(file.names, "permutation/selected.csv")

source("permutation/script/gentoma.R")
source("permutation/script/fixed_encoder.R")
set.seed(123)  # For reproducibility
num_seeds <- 5
maxlen <- 1000000
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
window.size <- 0.10
total.sub <- round(window.size * maxlen)
seg.len <- 2000 #round(0.001 * maxlen)
seg.num <- total.sub / seg.len
iterations <- 300
fidx <- 1
N <- 20    # Max generations
M <- 3     # Stagnant generations
lambda <- 1
update_rate <- 0.75
mutation_prob <- 0.2
print_setup(maxlen, seg.len, seg.num, iterations, fidx, N, M,
            update_rate, mutation_prob, pop.dense=NA, lambda, pred.orignal=NA)

library(parallel)

# Detect the number of cores available
numCores <- 20

# Create a cluster with the available cores
cl <- makeCluster(numCores)

# Export necessary variables and functions to the cluster nodes
clusterExport(cl, varlist = c("maxlen", "seg.len", "seg.num", "iterations", "fidx",
                              "N", "M", "lambda", "update_rate", "mutation_prob", "num_seeds"))

# Define the function to process each file
process_file <- function(fasta_file) {
  # Load required libraries within the function
  library(reticulate)
  required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                         "caret", "keras", "magrittr","openxlsx", "parallel", "rtracklayer")
  invisible(lapply(required_packages, library, character.only = TRUE))
  
  # Source necessary scripts
  source("permutation/script/gentoma.R")
  source("permutation/script/fixed_encoder.R")
  
  # Load the model within the function to ensure it's available on each node
  model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)
  
  fasta_path <- file.path(fasta_file)
  seeds <- sample(1e6, num_seeds)
  cat("Processing file:", fasta_file, "\n")
  
  # Initialize results list for this file
  results <- list()
  
  for (seed in seeds) {
    cat("Processing file:", fasta_file, "with seed:", seed, "\n")
    # Set the seed
    set.seed(seed)
    result <- tryCatch({
      gom_fasta(fasta_path, maxlen, seg.len, seg.num, iterations, model,
                fidx, N, M, lambda, update_rate, mutation_prob)
    }, error = function(e) {
      cat("Error processing file:", fasta_file, "with seed:", seed, "\n", e$message, "\n")
      NULL
    })
    # Store the result with the seed as the key
    results[[paste0("seed_", seed)]] <- result
  }
  # Return results for this file
  return(results)
}

# Export the 'process_file' function to cluster nodes
clusterExport(cl, varlist = c("process_file"))

# Process files in parallel using parLapply
output_list <- parLapply(cl, file.names, process_file)

# Assign names to output list
names(output_list) <- file.names

# Save the output
saveRDS(output_list, "permutation/strat_gentoma_incremental_cl.rds")

# Stop the cluster after processing
stopCluster(cl)