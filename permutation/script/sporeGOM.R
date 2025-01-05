### GENTOMA Usage ###

## Sporulation Example

#### Set up the environment ####
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "keras", "keras3", "magrittr", "openxlsx", "future", "future.apply")

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))
source("permutation/script/gentoma.R")
source("permutation/script/fixed_encoder.R")

#### Load the data ####
maxlen <- 1000000
target_from_csv <- "sporulation/sporeinfo.csv"
target_df <- read.csv(target_from_csv)
label_names <- names(target_df)[names(target_df) != "file"]
print(label_names)
folder_path <- "sporulation/genome"
model <- keras::load_model_hdf5("sporulation/spore_model.h5", compile = FALSE)

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
  gamma = 1,             # Penalty factor for coverage
  update_rate = 0.75,    # 75% of population is replaced each iteration
  mutation_prob = 0.2,   # Mutation probability per segment
  sample.count = 200
)

saveRDS(multiple_fasta_results, "permutation/spore_fasta_GOM.rds")