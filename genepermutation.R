codon.dict <- list(
  "A" = c("GCT", "GCC", "GCA", "GCG"),  # Alanine
  "C" = c("TGT", "TGC"),                # Cysteine
  "D" = c("GAT", "GAC"),                # Aspartic acid
  "E" = c("GAA", "GAG"),                # Glutamic acid
  "F" = c("TTT", "TTC"),                # Phenylalanine
  "G" = c("GGT", "GGC", "GGA", "GGG"),  # Glycine
  "H" = c("CAT", "CAC"),                # Histidine
  "I" = c("ATT", "ATC", "ATA"),         # Isoleucine
  "K" = c("AAA", "AAG"),                # Lysine
  "L" = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),  # Leucine
  "M" = c("ATG"),                       # Methionine
  "N" = c("AAT", "AAC"),                # Asparagine
  "P" = c("CCT", "CCC", "CCA", "CCG"),  # Proline
  "Q" = c("CAA", "CAG"),                # Glutamine
  "R" = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),  # Arginine
  "S" = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),  # Serine
  "T" = c("ACT", "ACC", "ACA", "ACG"),  # Threonine
  "V" = c("GTT", "GTC", "GTA", "GTG"),  # Valine
  "W" = c("TGG"),                       # Tryptophan
  "Y" = c("TAT", "TAC")                 # Tyrosine
)

stop.codons <- list("*" = c("TAA", "TAG", "TGA"))

all.dict <- c(codon.dict, stop.codons)

tokenize_triplets <- function(sequence) {
  # Preallocate space and avoid sapply
  n_triplets <- floor(length(sequence) / 3)
  triplets <- character(n_triplets)
  
  # Use vectorized operations for efficiency
  for (i in seq_len(n_triplets)) {
    start <- (i - 1) * 3 + 1
    triplets[i] <- paste(sequence[start:(start + 2)], collapse = "")
  }
  return(triplets)
}

triplets_keying <- function(triplets, dict=all.dict) {
  # Preallocate and vectorize key matching
  keys <- character(length(triplets))
  
  # Create a reverse lookup table for faster matching
  lookup <- unlist(lapply(names(dict), function(key) {
    setNames(rep(key, length(dict[[key]])), dict[[key]])
  }))
  
  keys <- lookup[triplets]
  keys[is.na(keys)] <- "X"
  return(keys)
}

permute_sequence <- function(sequence, type="ok",
                             min.subs, max.subs,
                             dict=codon.dict, spec.cond=FALSE,
                             spec.region=NULL) {
  mutable_sequence <- sequence
  keyed_sequence <- triplets_keying(mutable_sequence, all.dict)
  num.sub <- sample(min.subs:max.subs, 1)
  if (spec.cond) {
    sub.indices <- sample(spec.region, num.sub)
  } else {
    all_indices <- 1:length(mutable_sequence)
    # Check if '*' exists in keyed_sequence and find its index
    if ('*' %in% keyed_sequence) {
      k <- which(keyed_sequence == '*')
      # Remove the index k from all_indices
      all_indices <- all_indices[-k]
    }
    eligible_indices <- setdiff(all_indices, spec.region)
    sub.indices <- sample(eligible_indices, num.sub)
      #sample from not spec.region, and not the last one of the vector
      #vorher: sample(1:(length(mutable_sequence)-1), num.sub)
  }
  if (type == "ok") {
    replacements <- sapply(sub.indices, function(i) {
      sample(dict[[keyed_sequence[i]]], 1)
    })
  } else if (type == "func") {
    replacements <- sapply(sub.indices, function(i) {
      sample(unlist(dict), 1)
      # Exclude the current key's values from the list
      #available_elements <- unlist(dict[-which(names(dict) == keyed_sequence[i])])
      #sample(available_elements, 1)
    })
  } else {
    stop("Invalid type")
  }
  for (i in 1:length(sub.indices)) {
    mutable_sequence[sub.indices[i]] <- replacements[i]
  }
  return(mutable_sequence)
}

GenePermutation <- function(sequence, num.perm,
                            min.subs, max.subs,
                            dict=codon.dict,
                            spec.region=NULL) {
  permuted <- data.frame(
    seq = character(),
    label = character()
  )
  for (i in 1:num.perm) {
    label <- sample(c("normal", "abnormal", "special"), size = 1) #, prob = c(0.6, 0.2, 0.2))
    permuted_seq <- permute_sequence(sequence,
                                     type="ok", min.subs=min.subs,
                                     max.subs=max.subs,
                                     dict=dict, spec.cond=FALSE,
                                     spec.region=NULL)
    if (label == "abnormal") {
      permuted_seq <- permute_sequence(permuted_seq,
                                       type="func", min.subs=min.subs,
                                       max.subs=max.subs,
                                       dict=dict, spec.cond=FALSE,
                                       spec.region=spec.region)
    } else if (label == "special") {
      permuted_seq <- permute_sequence(permuted_seq,
                                       type="func", min.subs=min.subs,
                                       max.subs=max.subs,
                                       dict=dict, spec.cond=TRUE,
                                       spec.region=spec.region)
    }
    permuted <- rbind(permuted, data.frame(seq = I(list(permuted_seq)), label = label))
  }
  return(permuted)
}

### Function for inserting a motif
insert_motif <- function(sequence, motif) {
  seq_length <- nchar(sequence)
  motif_length <- nchar(motif)
  if (motif_length %% 3 != 0) {
    stop("Motif length must be a multiple of 3 to fit into codons.")
  }
  possible_start_positions <- seq(1, seq_length - motif_length + 1, by = 3)
  start_pos <- sample(possible_start_positions, 1)
  pre_motif <- substr(sequence, 1, start_pos - 1)
  post_motif <- substr(sequence, start_pos + motif_length, seq_length)
  new_sequence <- paste0(pre_motif, motif, post_motif)
  return(new_sequence)
}

# function that takes "repeats" as argument, replicates input for repeats times
# and performs ok-mutate on each replicate, and then inserts the motif
GeneMotif <- function(sequence, motif, repeats, min.subs, max.subs) {
  result <- data.frame(
    seq = character(),
    label = character()
  )
  for (i in 1:repeats) {
    motif_seq <- sequence
    motif_split <- strsplit(sequence, "")[[1]]
    label <- sample(c("motif", "normal"), 1)
    motif_trip <- tokenize_triplets(motif_split)
    motif_seq <- permute_sequence(motif_trip, type="ok", min.subs=min.subs,
                                  max.subs=max.subs, dict=codon.dict, spec.cond=FALSE,
                                  spec.region=NULL)
    motif_seq <- paste(motif_seq, collapse = "")
    if (label == "motif") {
      motif_seq <- insert_motif(motif_seq, motif)
    }
    result <- rbind(result, data.frame(seq = I(list(motif_seq)), label))
  }
  return(result)
}

# finding the motif:
find_motif <- function(sequence, motif) {
  seq_length <- nchar(sequence)
  motif_length <- nchar(motif)
  
  positions <- c()
  
  for (i in 1:(seq_length - motif_length + 1)) {
    sub_seq <- substr(sequence, i, i + motif_length - 1)
    if (sub_seq == motif) {
      positions <- c(positions, i)
    }
  }
  return(positions)
}
