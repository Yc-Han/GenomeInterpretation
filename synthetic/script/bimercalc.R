### Setup ###

required_packages <- c("microseq", "HMM", "parallel", "openxlsx", "tidyverse", "doParallel", "foreach", "rtracklayer", "data.table")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)

# Load all the packages
invisible(lapply(required_packages, library, character.only = TRUE))

### Parallel setup ###
# Detect the number of available cores
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)

### Order-Specific Transmission Probabilities ###

# Read and filter the microbes data
microbes <- read.xlsx("sporulation/microbe.cards table S1.xlsx")
microbes <- microbes %>%
  mutate(Gff.file = gsub(".fasta", ".gff", Fasta.file)) %>%
  filter(Member.of.WA.subset == TRUE & !is.na(Fasta.file) & !is.na(Gff.file))

Order_counts <- table(microbes$Order)
phyla_to_keep <- names(Order_counts[Order_counts > 5])  # Orders with more than 5 genomes
microbes_filtered <- microbes %>%
  filter(Order %in% phyla_to_keep)

# Sample up to 5 genomes per Order
set.seed(123)  # For reproducibility
sampled_files <- list()
for (phy in phyla_to_keep) {
  phy_data <- microbes_filtered %>%
    filter(Order == phy)
  num_genomes <- nrow(phy_data)
  sampled_indices <- if (num_genomes >= 5) sample(1:num_genomes, 5) else 1:num_genomes
  sampled_files[[phy]] <- paste0("sporulation/genome/", phy_data$Fasta.file[sampled_indices])
}

# Prepare the list of file paths
phyla_file_list <- lapply(names(sampled_files), function(phy) {
  list(Order = phy, fasta = sampled_files[[phy]])
})

# Update GFF file paths
phyla_file_list <- lapply(phyla_file_list, function(x) {
  x$gff <- gsub(".fasta", ".gff", x$fasta)
  x$gff <- gsub("genome/", "gff/gff_new/", x$gff)
  x
})

# Parallelized version of the loop using foreach
result_list <- foreach(i = seq_along(phyla_file_list), .packages = c("microseq", "data.table", "rtracklayer")) %dopar% {
  # Initialize accumulators
  total_transition_counts <- matrix(0, nrow = 2, ncol = 2,
                                    dimnames = list(c("coding", "non-coding"),
                                                    c("coding", "non-coding")))
  # Initialize nucleotide counts
  total_nucleotide_counts_coding <- setNames(rep(0, 4), c("A", "C", "G", "T"))
  total_nucleotide_counts_noncoding <- setNames(rep(0, 4), c("A", "C", "G", "T"))
  
  # Initialize bimer counts
  possible_bimers <- as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0))
  bimer_counts_coding <- setNames(rep(0, length(possible_bimers)), possible_bimers)
  bimer_counts_noncoding <- setNames(rep(0, length(possible_bimers)), possible_bimers)
  
  # Initialize total positions
  total_positions_coding <- 0
  total_positions_noncoding <- 0
  
  sequence_lengths <- numeric()  # Collect all sequence lengths
  # Get the number of sequences in this Order
  num_sequences <- length(phyla_file_list[[i]]$fasta)
  
  # Loop over each sequence in the Order
  for (j in seq_len(num_sequences)) {
    fasta_path <- phyla_file_list[[i]]$fasta[j]
    gff_path <- phyla_file_list[[i]]$gff[j]
    
    # Read the sequence
    seq <- microseq::readFasta(fasta_path)$Sequence[1]
    seq <- gsub("[^ACGT]", "", seq)
    sequence_lengths <- c(sequence_lengths, nchar(seq))
    # Read the GFF data
    gff <- readGFF(gff_path)
    list_data <- gff@listData
    gff.dt <- data.table(
      start = list_data$start,
      end = list_data$end,
      Name = list_data$Name,
      gene = list_data$gene,
      product = list_data$product
    ) %>%
      filter(!is.na(product) & !is.na(gene)) %>%
      mutate(length = end - start + 1) %>%
      filter(length > 100 & length < 2000)
    
    # Create the nucleotide data table
    seq.dt <- data.table(
      nucleotide = unlist(strsplit(seq, ""))
    )
    seq.dt[, position := .I]
    seq.dt[, `:=`(coding = 0, non.coding = 1)]
    seq.dt[, `:=`(start = position, end = position)]
    
    # Set keys for foverlaps
    setkey(seq.dt, start, end)
    setkey(gff.dt, start, end)
    
    # Find overlaps to identify coding regions
    overlap <- foverlaps(seq.dt, gff.dt, nomatch = 0)
    seq.dt[overlap$position, `:=`(coding = 1, non.coding = 0)]
    seq.dt[, `:=`(start = NULL, end = NULL)]
    
    # Create bimers
    seq.dt[, bimer := paste0(nucleotide, data.table::shift(nucleotide, type = "lead"))]
    seq.dt <- seq.dt[!is.na(bimer)]  # Remove last row where bimer is NA
    
    # Count bimers in coding and noncoding regions
    # For coding regions
    bimer_counts_seq_coding <- seq.dt[coding == 1, .N, by = bimer]
    bimer_counts_coding[bimer_counts_seq_coding$bimer] <- bimer_counts_coding[bimer_counts_seq_coding$bimer] + bimer_counts_seq_coding$N
    
    # For noncoding regions
    bimer_counts_seq_noncoding <- seq.dt[coding == 0, .N, by = bimer]
    bimer_counts_noncoding[bimer_counts_seq_noncoding$bimer] <- bimer_counts_noncoding[bimer_counts_seq_noncoding$bimer] + bimer_counts_seq_noncoding$N
    
    # Update total positions (number of bimers)
    positions_coding <- sum(seq.dt$coding == 1)
    positions_noncoding <- sum(seq.dt$coding == 0)
    total_positions_coding <- total_positions_coding + positions_coding
    total_positions_noncoding <- total_positions_noncoding + positions_noncoding
    
    # Count nucleotides in coding and noncoding regions
    # For coding regions
    nucleotide_counts_coding <- seq.dt[coding == 1, .N, by = nucleotide]
    total_nucleotide_counts_coding[nucleotide_counts_coding$nucleotide] <- total_nucleotide_counts_coding[nucleotide_counts_coding$nucleotide] + nucleotide_counts_coding$N
    
    # For noncoding regions
    nucleotide_counts_noncoding <- seq.dt[coding == 0, .N, by = nucleotide]
    total_nucleotide_counts_noncoding[nucleotide_counts_noncoding$nucleotide] <- total_nucleotide_counts_noncoding[nucleotide_counts_noncoding$nucleotide] + nucleotide_counts_noncoding$N
    
    # Calculate state transitions
    seq.dt[, state := ifelse(coding == 1, "coding", "non-coding")]
    seq.dt[, next_state := data.table::shift(state, type = "lead")]
    seq.dt <- seq.dt[!is.na(next_state)]
    
    # Count transitions for this sequence
    transition_counts <- seq.dt[, .N, by = .(state, next_state)]
    
    # Initialize a counts matrix for this sequence
    counts_matrix <- matrix(0, nrow = 2, ncol = 2,
                            dimnames = list(c("coding", "non-coding"),
                                            c("coding", "non-coding")))
    # Populate the counts matrix
    for (k in seq_len(nrow(transition_counts))) {
      from_state <- transition_counts$state[k]
      to_state <- transition_counts$next_state[k]
      counts <- transition_counts$N[k]
      counts_matrix[from_state, to_state] <- counts_matrix[from_state, to_state] + counts
    }
    
    # Accumulate transition counts
    total_transition_counts <- total_transition_counts + counts_matrix
  }
  
  # Calculate average transition probabilities for the Order
  total_transitions_from_state <- rowSums(total_transition_counts)
  total_transitions_from_state[total_transitions_from_state == 0] <- 1  # Avoid division by zero
  transprob_mat <- sweep(total_transition_counts, 1, total_transitions_from_state, "/")
  
  # Calculate nucleotide densities for coding and noncoding regions
  nt_density_coding <- total_nucleotide_counts_coding / sum(total_nucleotide_counts_coding, na.rm = TRUE)
  nt_density_noncoding <- total_nucleotide_counts_noncoding / sum(total_nucleotide_counts_noncoding, na.rm = TRUE)
  
  # Calculate bimer densities for coding and noncoding regions
  bimer_density_coding <- bimer_counts_coding / sum(bimer_counts_coding, na.rm = TRUE)
  bimer_density_noncoding <- bimer_counts_noncoding / sum(bimer_counts_noncoding, na.rm = TRUE)
  avg_sequence_length <- mean(sequence_lengths)
  # Store the results in the result list
  list( avg_sequence_length = avg_sequence_length,
        transprob_mat = transprob_mat,
        nt_density_coding = nt_density_coding,
        nt_density_noncoding = nt_density_noncoding,
        bimer_density_coding = bimer_density_coding,
        bimer_density_noncoding =
          bimer_density_noncoding )
}

# optional: name the list
# paste the corresponding phylum to names(Order_counts[Order_counts > 5])
phylum_of_order <- microbes_filtered %>%
  filter(Order %in% phyla_to_keep) %>%
  select(Order, Phylum) %>%
  distinct() %>%
  arrange(Order)
phy_order_names <- paste0(phylum_of_order$Phylum, "_", phylum_of_order$Order)


names(result_list) <- phy_order_names

stopCluster(cl)

clean_result_list <- function(result) {
  # Clean each sublist by removing NAs in all the components
  result$transprob_mat <- result$transprob_mat[complete.cases(result$transprob_mat), ]
  result$nt_density_coding <- result$nt_density_coding[complete.cases(result$nt_density_coding)]
  result$nt_density_noncoding <- result$nt_density_noncoding[complete.cases(result$nt_density_noncoding)]
  result$bimer_density_coding <- result$bimer_density_coding[complete.cases(result$bimer_density_coding)]
  result$bimer_density_noncoding <- result$bimer_density_noncoding[complete.cases(result$bimer_density_noncoding)]
  
  return(result)
}

# Apply the cleaning function to each element of result_list
result_list <- lapply(result_list, clean_result_list)
saveRDS(result_list, "permutation/estimates_parallel.rds")
