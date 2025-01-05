# creating entirely synthetic sequences

# source genepermutation.R
source("genepermutation.R")

# SeqSyn: a function that takes length, n, dict as arguments and synthesizes a
# gene sequence of length using codons from dict.
SeqSyn <- function(length, n = 1, dict = codon.dict, by.codon = TRUE) {
  sequences <- vector("list", n)
  for (j in 1:n) {
    sequence <- c()
    if (by.codon) {
      for (i in 1:(length / 3)) {
        amino_acid <- sample(names(dict), 1)
        codon <- sample(dict[[amino_acid]], 1)
        sequence <- c(sequence, codon)
      }
    } else {
      bases <- c("A", "C", "G", "T")
      sequence <- sample(bases, length, replace = TRUE)
    }
    sequences[[j]] <- paste(sequence, collapse = "")
  }
  return(sequences)
}

