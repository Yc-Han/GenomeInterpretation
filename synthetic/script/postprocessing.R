library(rentrez)
library(seqinr)

## gram-staining genes
# Relevant genes related to Gram-staining characteristics
# gene_names <- c("murA", "murB", "murC", "murD", "murE", "murF", "pbp1", "pbp2", "pbp3", "ddl", "ftsI")
gene_names <- c(
  # Ribosomal Protein Genes (Small Subunit)
  "rpsA", "rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG",
  # Ribosomal Protein Genes (Large Subunit)
  "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplG",
  # Molecular Chaperones
  "groEL", "groES", "dnaK", "dnaJ", "grpE", "hsp60", "hsp70"
)

# Vector of bacterial orders
bacterial_orders <- c(
  "Actinomycetota_Actinomycetales",     "Pseudomonadota_Aeromonadales",
  "Pseudomonadota_Alteromonadales",     "Bacillota_Bacillales",
  "Bacteroidota_Bacteroidales",         "Actinomycetota_Bifidobacteriales",
  "Pseudomonadota_Burkholderiales",     "Campylobacterota_Campylobacterales",
  "Pseudomonadota_Caulobacterales",     "Pseudomonadota_Cellvibrionales",
  "Pseudomonadota_Chromatiales",        "Pseudomonadota_Enterobacterales",
  "Bacillota_Eubacteriales",            "Bacteroidota_Flavobacteriales",
  "Euryarchaeota_Halobacteriales",      "Pseudomonadota_Hyphomicrobiales",
  "Actinomycetota_Kitasatosporales",    "Bacillota_Lachnospirales",
  "Bacillota_Lactobacillales",          "Pseudomonadota_Legionellales",
  "Spirochaetota_Leptospirales",        "Pseudomonadota_Lysobacterales",
  "Actinomycetota_Micrococcales",       "Actinomycetota_Micromonosporales",
  "Pseudomonadota_Moraxellales",        "Actinomycetota_Mycobacteriales",
  "Mycoplasmatota_Mycoplasmoidales",    "Myxococcota_Myxococcales",
  "Pseudomonadota_Neisseriales",        "Pseudomonadota_Oceanospirillales",
  "Pseudomonadota_Pasteurellales",      "Actinomycetota_Propionibacteriales",
  "Pseudomonadota_Pseudomonadales",     "Actinomycetota_Pseudonocardiales",
  "Pseudomonadota_Rhodospirillales",    "Pseudomonadota_Rickettsiales",
  "Pseudomonadota_Sphingomonadales",    "Spirochaetota_Spirochaetales",
  "Actinomycetota_Streptosporangiales", "Pseudomonadota_Vibrionales"
)

# Create a named vector of iconic species for each order
# Note: The species are selected as common representatives of their respective orders.
iconic_species <- c(
  "Streptomyces coelicolor",    # Actinomycetota_Actinomycetales
  "Aeromonas hydrophila",       # Pseudomonadota_Aeromonadales
  "Alteromonas macleodii",      # Pseudomonadota_Alteromonadales
  "Bacillus subtilis",          # Bacillota_Bacillales
  "Bacteroides fragilis",       # Bacteroidota_Bacteroidales
  "Bifidobacterium longum",     # Actinomycetota_Bifidobacteriales
  "Burkholderia cepacia",       # Pseudomonadota_Burkholderiales
  "Campylobacter jejuni",       # Campylobacterota_Campylobacterales
  "Caulobacter crescentus",     # Pseudomonadota_Caulobacterales
  "Vibrio cholerae",            # Pseudomonadota_Cellvibrionales
  "Allochromatium vinosum",     # Pseudomonadota_Chromatiales
  "Escherichia coli",           # Pseudomonadota_Enterobacterales
  "Clostridium perfringens",    # Bacillota_Eubacteriales
  "Flavobacterium johnsoniae",  # Bacteroidota_Flavobacteriales
  "Halobacterium salinarum",    # Euryarchaeota_Halobacteriales
  "Rhizobium leguminosarum",    # Pseudomonadota_Hyphomicrobiales
  "Kitasatospora setae",        # Actinomycetota_Kitasatosporales
  "Lachnospira multipara",      # Bacillota_Lachnospirales
  "Lactobacillus acidophilus",  # Bacillota_Lactobacillales
  "Legionella pneumophila",     # Pseudomonadota_Legionellales
  "Leptospira interrogans",     # Spirochaetota_Leptospirales
  "Lysobacter enzymogenes",     # Pseudomonadota_Lysobacterales
  "Micrococcus luteus",         # Actinomycetota_Micrococcales
  "Micromonospora aurantiaca",  # Actinomycetota_Micromonosporales
  "Acinetobacter baumannii",    # Pseudomonadota_Moraxellales
  "Mycobacterium tuberculosis", # Actinomycetota_Mycobacteriales
  "Mycoplasma pneumoniae",      # Mycoplasmatota_Mycoplasmoidales
  "Myxococcus xanthus",         # Myxococcota_Myxococcales
  "Neisseria meningitidis",     # Pseudomonadota_Neisseriales
  "Marinomonas mediterranea",   # Pseudomonadota_Oceanospirillales
  "Haemophilus influenzae",     # Pseudomonadota_Pasteurellales
  "Propionibacterium acnes",    # Actinomycetota_Propionibacteriales
  "Pseudomonas aeruginosa",     # Pseudomonadota_Pseudomonadales
  "Pseudonocardia dioxanivorans",# Actinomycetota_Pseudonocardiales
  "Rhodospirillum rubrum",      # Pseudomonadota_Rhodospirillales
  "Rickettsia rickettsii",      # Pseudomonadota_Rickettsiales
  "Sphingomonas paucimobilis",  # Pseudomonadota_Sphingomonadales
  "Treponema pallidum",         # Spirochaetota_Spirochaetales
  "Streptosporangium roseum",   # Actinomycetota_Streptosporangiales
  "Vibrio parahaemolyticus"     # Pseudomonadota_Vibrionales
)

# Ensure that the species vector is named with the corresponding orders
names(iconic_species) <- bacterial_orders

# Loop over each species and download the genes
for (order in names(iconic_species)) {
  species <- iconic_species[order]
  cat("Processing species:", species, "from order:", order, "\n")
  
  # Replace spaces with "+" for URL encoding in search terms
  species_query <- gsub(" ", "+", species)
  
  # Create a directory for the species
  species_dir <- file.path("permutation/ncbi", gsub(" ", "_", order))
  dir.create(species_dir, showWarnings = FALSE)
  
  for (gene in gene_names) {
    cat("  Downloading gene:", gene, "\n")
    
    # Construct the search term
    search_term <- paste0(species_query, "[Organism] AND ", gene, "[Gene] AND RefSeq[Filter] AND CDS[Feature Key]")
    
    # Search NCBI Nucleotide database for the specific gene CDS, not the whole genome
    search_result <- entrez_search(db = "nucleotide", term = search_term, retmax = 1)
    
    if (length(search_result$ids) > 0) {
      # Fetch the gene sequence in FASTA format
      sequence <- entrez_fetch(db = "nucleotide", id = search_result$ids[1], rettype = "fasta_cds_na")
      
      # Save the sequence to a file
      gene_filename <- paste0(gene, ".fasta")
      gene_filepath <- file.path(species_dir, gene_filename)
      write(sequence, file = gene_filepath)
    } else {
      cat("    Gene", gene, "not found for species", species, "\n")
    }
  }
}

# Function to process a saved FASTA file and retain only the correct sequence
process_fasta_file <- function(file_path, gene_name) {
  # Read the entire FASTA file
  fasta_lines <- readLines(file_path)
  
  # Initialize an empty vector to hold the processed lines
  processed_fasta <- c()
  
  # Variable to keep track of whether the correct sequence is found
  is_correct_sequence <- FALSE
  
  # Loop through the FASTA file line by line
  for (line in fasta_lines) {
    # If the line starts with '>', it's a header line
    if (startsWith(line, ">")) {
      # Check if the header contains the correct gene name
      if (grepl(gene_name, line, ignore.case = F)) {
        # If it contains the correct gene, set flag to TRUE and start collecting the sequence
        is_correct_sequence <- TRUE
        processed_fasta <- c(processed_fasta, line)  # Include the header line
      } else {
        # If it's a different gene, set flag to FALSE and skip the sequence
        is_correct_sequence <- FALSE
      }
    } else {
      # If it's a sequence line and the correct sequence flag is set to TRUE, include it
      if (is_correct_sequence) {
        processed_fasta <- c(processed_fasta, line)
      }
    }
  }
  
  # If the correct sequence was found, overwrite the file with the processed version
  if (length(processed_fasta) > 0) {
    writeLines(processed_fasta, file_path)
  } else {
    # Optionally, you can print a warning if no correct sequence was found
    cat("Warning: No correct sequence found in", file_path, "\n")
  }
}

# Example usage: Process all FASTA files for their corresponding gene names
# List all files in the directory where gene sequences are saved
fasta_files <- list.files("permutation/ncbi/", recursive = TRUE, pattern = "\\.fasta$", full.names = TRUE)

# Apply the function to each file
for (file_path in fasta_files) {
  # Extract the gene name from the filename (assuming it's in the filename)
  gene_name <- tools::file_path_sans_ext(basename(file_path))  # Removes extension to get the gene name
  
  # Process the file to retain only the correct sequence
  process_fasta_file(file_path, gene_name)
}

# delete all fasta files that are > 5 kb
fasta_files <- list.files("permutation/ncbi/", recursive = TRUE, pattern = "\\.fasta$", full.names = TRUE)
# Loop over the files and check their sizes
for (file in fasta_files) {
  # Get the file size in bytes
  file_info <- file.info(file)
  file_size_kb <- file_info$size / 1024  # Convert to KB
  
  # If the file is larger than 5 KB, delete it
  if (file_size_kb > 5) {
    cat("Deleting file:", file, "Size:", file_size_kb, "KB\n")
    file.remove(file)  # Delete the file
  }
}

set.seed(42) 

labelled_yes_orders <- sample(bacterial_orders, 20)

startcodon <- "ATG"

# list all files in permutation/hmmdata/
sourceseqs <- list.files("permutation/hmmdata/", full.names = TRUE)

source("genepermutation.R")
# Mutator:
# 1. each trimer has 0.0001 probability of being synonymously mutated
# 2. random insertion (of a random length randomly sampled ACGT seq)
# 3. random deletion (of a random length)
# 4. random replacement (of a random length randomly sampled ACGT seq)
# 2,3,4 combined 0.5 prob of happening for every 1 Mbp of sequence
# (seqs is usually much longer than 1 Mbp)

mutate_sequence <- function(
    sequence,
    codon.dict,  # The dictionary you provided
    codon_mutation_rate = 1e-3,       # Probability of synonymous mutation per codon
    event_rate_per_Mbp  = 0.5,      # Probability for insertion/deletion/replacement per Mb
    max_indel_length    = 5000         # Max length for insertion/deletion/replacement
) {
  #----------------------------#
  # 1) SYNONYMOUS CODON MUTATION
  #----------------------------#
  # Tokenize into codons
  triplets <- tokenize_triplets(sequence)
  # Identify amino-acid keys for each triplet
  keys <- triplets_keying(triplets, dict = all.dict)
  
  # Mutate each codon with probability codon_mutation_rate
  n_codons <- length(triplets)
  mutate_flags <- runif(n_codons) < codon_mutation_rate
  
  # For codons flagged to mutate:
  #  - Sample from the same amino-acid set in codon.dict[[key]]
  #  - Exclude the original codon itself if possible
  for (i in which(mutate_flags)) {
    aa_key <- keys[i]
    # If we have a valid amino acid key in codon.dict
    if (aa_key %in% names(codon.dict)) {
      original_codon <- triplets[i]
      possible_syns <- codon.dict[[aa_key]]
      if (length(possible_syns) > 1) {
        # Pick a new synonymous codon different from the original
        new_codon <- original_codon
        while (new_codon == original_codon) {
          new_codon <- sample(possible_syns, 1)
        }
        triplets[i] <- new_codon
      }
    }
  }
  
  # Re-assemble mutated sequence (codon string -> single-character vector)
  mutated_seq_vector <- unlist(strsplit(paste(triplets, collapse = ""), split = ""))
  
  #----------------------------#
  # 2) INSERTIONS, DELETIONS, REPLACEMENTS
  #----------------------------#
  # For each 1 Mb (1e6 bp), there's event_rate_per_Mbp chance
  # => Expected number of each event is:
  seq_length <- length(mutated_seq_vector)
  expected_events <- event_rate_per_Mbp * (seq_length / 1e6)
  
  # Draw how many insertions, deletions, replacements to perform (Poisson)
  n_insert   <- rpois(1, lambda = expected_events)
  n_delete   <- rpois(1, lambda = expected_events)
  n_replace  <- rpois(1, lambda = expected_events)
  
  # Helper to generate random nucleotides
  random_nucleotides <- function(k) {
    sample(c("A","C","G","T"), size = k, replace = TRUE)
  }
  
  #----- Perform Insertions -----#
  for (i in seq_len(n_insert)) {
    if (seq_length == 0) {
      # If sequence is empty, insertion is just appended
      mutated_seq_vector <- c(mutated_seq_vector, random_nucleotides(sample.int(max_indel_length, 1)))
      seq_length <- length(mutated_seq_vector)
      next
    }
    pos <- sample.int(seq_length + 1, 1)  # insertion can be between positions
    len_ins <- sample.int(max_indel_length, 1)
    insert_bases <- random_nucleotides(len_ins)
    
    mutated_seq_vector <- c(
      mutated_seq_vector[1:(pos - 1)],
      insert_bases,
      mutated_seq_vector[pos:seq_length]
    )
    seq_length <- length(mutated_seq_vector)
  }
  
  #----- Perform Deletions -----#
  for (i in seq_len(n_delete)) {
    if (seq_length == 0) break
    pos <- sample.int(seq_length, 1)
    len_del <- sample.int(max_indel_length, 1)
    end_del <- min(pos + len_del - 1, seq_length)
    
    mutated_seq_vector <- c(
      mutated_seq_vector[1:(pos - 1)],
      mutated_seq_vector[(end_del + 1):seq_length]
    )
    seq_length <- length(mutated_seq_vector)
    if (seq_length <= 0) break
  }
  
  #----- Perform Replacements -----#
  for (i in seq_len(n_replace)) {
    if (seq_length == 0) break
    pos <- sample.int(seq_length, 1)
    len_rep <- sample.int(max_indel_length, 1)
    end_rep <- min(pos + len_rep - 1, seq_length)
    # Generate random replacement
    replacement <- random_nucleotides(end_rep - pos + 1)
    
    mutated_seq_vector <- c(
      mutated_seq_vector[1:(pos - 1)],
      replacement,
      mutated_seq_vector[(end_rep + 1):seq_length]
    )
    seq_length <- length(mutated_seq_vector)
    if (seq_length <= 0) break
  }
  
  # Return final mutated sequence (as a character vector)
  return(mutated_seq_vector)
}

# Inserter
insert_downloaded_genes <- function(
    seqs,
    seqpath,
    labelled_yes_orders,
    startcodon = "ATG",
    ribosomal_genes = c(
      # Ribosomal Protein Genes (Small Subunit)
      "rpsA", "rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG",
      # Ribosomal Protein Genes (Large Subunit)
      "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplG"
    ),
    chaperonin_genes = c("groEL", "groES", "dnaK", "dnaJ", "grpE", "hsp60", "hsp70")
) {
  #--- 2. Determine the order name from the filename (remove ".rds")
  check <- gsub("\\.rds$", "", basename(seqpath))
  is_yes_label <- check %in% labelled_yes_orders
  
  #--- 3. List all .fasta files for this order
  #     (assuming you've already downloaded them into "permutation/ncbi/<order>/*.fasta")
  fasta_dir <- file.path("permutation", "ncbi", check)
  if (!dir.exists(fasta_dir)) {
    # If no directory, simply return the original seq + empty log
    return(list(
      mutated_seq = seqs,
      insertion_log = data.frame(
        seqpath   = character(),
        gene_file = character(),
        start     = integer(),
        end       = integer(),
        stringsAsFactors = FALSE
      )
    ))
  }
  all_fasta_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  # Separate into two families by matching filenames against known gene names
  ribo_files <- all_fasta_files[grepl(paste(ribosomal_genes, collapse = "|"), all_fasta_files)]
  chap_files <- all_fasta_files[grepl(paste(chaperonin_genes, collapse = "|"), all_fasta_files)]
  
  #--- 4. We will insert genes *per 1 Mbp* of sequence
  seq_length <- nchar(seqs)
  n_megablocks <- floor(seq_length / 1e6)  # how many full 1 Mbp segments
  
  # For tracking where we've already inserted (to avoid overlap)
  inserted_regions <- integer(0)
  
  # A data frame to record insertions (start, end, gene, etc.)
  insertion_log <- data.frame(
    seqpath   = character(),
    gene_file = character(),
    start     = integer(),
    end       = integer(),
    stringsAsFactors = FALSE
  )
  
  # Helper function to find a random non-overlapping start-codon position
  pick_start_position <- function(seqs, inserted_regions, startcodon) {
    start_positions <- gregexpr(startcodon, seqs, fixed = TRUE)[[1]]
    if (all(start_positions < 0)) return(NA)  # No matches at all
    
    # Exclude any start positions that overlap with prior insertions
    valid_positions <- setdiff(start_positions, inserted_regions)
    if (length(valid_positions) == 0) return(NA)
    
    sample(valid_positions, 1)
  }
  
  #--- 5. Process each 1 Mbp chunk
  for (mb in seq_len(n_megablocks)) {
    # (Re)compute boundaries in the CURRENT mutated sequence
    current_len <- nchar(seqs)
    mb_start <- (mb - 1)*1e6 + 1
    mb_end   <- min(mb * 1e6, current_len)
    
    # If our start goes beyond the current sequence, nothing to do
    if (mb_start > mb_end) {
      break
    }
    
    # Decide how many total genes to insert in this Mb
    if (is_yes_label) {
      # 2–6 total, half from ribosomal, half from chaperonin
      total_to_insert <- sample(2:6, 1)
      n_ribo <- floor(total_to_insert / 2)
      n_chap <- total_to_insert - n_ribo
    } else {
      # 0–3 from one of the two families
      total_to_insert <- sample(0:3, 1)
      if (total_to_insert == 0) {
        next  # no insertions in this chunk
      }
      # pick randomly which family to use
      if (runif(1) < 0.5) {
        n_ribo <- total_to_insert
        n_chap <- 0
      } else {
        n_ribo <- 0
        n_chap <- total_to_insert
      }
    }
    
    # A helper to pick one ATG in the current chunk that does NOT overlap with inserted_regions
    pick_start_in_chunk <- function(seqs, chunk_start, chunk_end, inserted_regions, startcodon) {
      # Extract substring for the chunk
      chunk_substr <- substring(seqs, chunk_start, chunk_end)
      
      # Find all local startcodon positions
      local_positions <- gregexpr(startcodon, chunk_substr, fixed = TRUE)[[1]]
      local_positions <- local_positions[local_positions != -1]
      if (!length(local_positions)) {
        return(NA_integer_)  # No start codons in this chunk
      }
      
      # Convert local positions to absolute positions in the full sequence
      absolute_positions <- local_positions + chunk_start - 1
      
      # Exclude any positions already used by previous insertions
      valid_positions <- setdiff(absolute_positions, inserted_regions)
      if (!length(valid_positions)) {
        return(NA_integer_)
      }
      # Pick one
      sample(valid_positions, 1)
    }
    
    #--------------------#
    # Insert from RIBOSOMAL family
    #--------------------#
    if (n_ribo > 0 && length(ribo_files) > 0) {
      chosen_ribo <- sample(ribo_files, n_ribo, replace = (length(ribo_files) < n_ribo))
      for (fpath in chosen_ribo) {
        pos <- pick_start_in_chunk(
          seqs, mb_start, mb_end, inserted_regions, startcodon
        )
        if (is.na(pos)) {
          # No valid position in chunk; skip
          next
        }
        
        # Read the gene to insert
        fasta_data <- read.fasta(fpath, as.string = TRUE)
        fasta_seq <- paste(getSequence(fasta_data)[[1]], collapse = "")
        gene_len <- nchar(fasta_seq)
        
        # Perform the insertion
        seqs <- paste0(
          substring(seqs, 1, pos - 1),
          fasta_seq,
          substring(seqs, pos + gene_len, nchar(seqs))
        )
        
        # Mark new region: [pos, pos + gene_len - 1]
        new_region <- seq(pos, length.out = gene_len)
        inserted_regions <- c(inserted_regions, new_region)
        
        # Log the insertion
        insertion_log <- rbind(insertion_log, data.frame(
          seqpath   = seqpath,
          gene_file = basename(fpath),
          start     = pos,
          end       = pos + gene_len - 1,
          stringsAsFactors = FALSE
        ))
        
        # Because sequence length changed, update mb_end if needed
        # (to keep the same chunk size in the *mutated* sequence)
        delta_len <- gene_len  # how much we inserted
        # If pos <= mb_end, everything from pos onward is shifted by delta_len
        if (pos <= mb_end) {
          mb_end <- mb_end + delta_len
        }
      }
    }
    
    #--------------------#
    # Insert from CHAPERONIN family
    #--------------------#
    if (n_chap > 0 && length(chap_files) > 0) {
      chosen_chap <- sample(chap_files, n_chap, replace = (length(chap_files) < n_chap))
      for (fpath in chosen_chap) {
        pos <- pick_start_in_chunk(
          seqs, mb_start, mb_end, inserted_regions, startcodon
        )
        if (is.na(pos)) {
          next
        }
        fasta_data <- read.fasta(fpath, as.string = TRUE)
        fasta_seq <- paste(getSequence(fasta_data)[[1]], collapse = "")
        gene_len <- nchar(fasta_seq)
        
        seqs <- paste0(
          substring(seqs, 1, pos - 1),
          fasta_seq,
          substring(seqs, pos + gene_len, nchar(seqs))
        )
        
        new_region <- seq(pos, length.out = gene_len)
        inserted_regions <- c(inserted_regions, new_region)
        
        insertion_log <- rbind(insertion_log, data.frame(
          seqpath   = seqpath,
          gene_file = basename(fpath),
          start     = pos,
          end       = pos + gene_len - 1,
          stringsAsFactors = FALSE
        ))
        
        # Adjust mb_end if insertion happened before the current chunk boundary
        if (pos <= mb_end) {
          mb_end <- mb_end + gene_len
        }
      }
    }
  }
  
  #--- 6. Return final mutated sequence + insertion log
  list(
    mutated_seq   = seqs,
    insertion_log = insertion_log
  )
}

# read in all files
# Read in all files and process each one
for (seqpath in sourceseqs) {
  # Read the .rds file and extract sequences
  source <- readRDS(seqpath)
  seqs <- source$observation
  # split seqs into a single-nucleotide vector
  seqs <- unlist(strsplit(seqs, ""))
  # mutate the sequence 4 times, append to the original sequence
  # so seqs will in the end be 5 times the original sequence long
  for (i in 1:4) {
    mutated_seq <- mutate_sequence(seqs, codon.dict)
    seqs <- c(seqs, mutated_seq)
    # print i-th mutation done
    cat("Mutation", i, "done for:", basename(seqpath), "\n")
  }
  seqs <- paste(seqs, collapse = "")
  # insert genes into the sequence
  result <- insert_downloaded_genes(seqs, seqpath, labelled_yes_orders)
  # print insertion done
  cat("Insertion done for:", basename(seqpath), "\n")
  # save the insertion log to a new .rds file in synthetic/insertion_log/
  saveRDS(result$insertion_log, file = paste0("synthetic/insertion_log/", basename(seqpath)))
  # save the mutated sequence to a new .fasta file in synthetic/data/
  write.fasta(sequences = result$mutated_seq, names = basename(seqpath), file.out = paste0("synthetic/data/", basename(seqpath), ".fasta"))
}
