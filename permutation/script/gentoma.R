### Setup ###
required_packages <- c("deepG", "tidyverse", "microseq", "data.table", "seqinr",
                       "caret", "keras", "magrittr","openxlsx", "parallel", "rtracklayer")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

# required parameters:
# instance: DNA sequence
# onehot_instance: one-hot encoded instance sequence
# maxlen: maximum length of the model
# seg.len: length of the segments in the initial population
# seg.num: number of segments in the initial population
# iterations: number of iterations / individuals in the population
# model: deepG model
# fidx: index of the false prediction

#### Genetic Algorithm ####

## Function for detecting nucleotide density
nt_density <- function(gene_seq) {
  nt_bases <- c("A", "C", "G", "T")
  # if "N" is present in the sequence, remove it
  gene_seq <- gsub("N", "", gene_seq)
  gene_chars <- unlist(strsplit(gene_seq, NULL))
  gene_factor <- factor(gene_chars, levels = nt_bases)
  nt_counts <- table(gene_factor)
  pop.dense <- nt_counts / length(gene_chars)
  pop.dense <- as.numeric(pop.dense)
  return(pop.dense) # the background nucleotide density / noise
}

## Function for calculating Eff.pred
## Eff.pred: the effective prediction value,
## I(flipped prediction) / (1 + lambda * sqrt(coverage))
## coverage: sum of window lengths / maxlen
coverage <- function(ends, starts, maxlen) {
  return(sum(ends - starts) / maxlen)
}
eff_pred <- function(false.pred, starts, ends, maxlen, lambda = 1) {
  cov <- coverage(ends, starts, maxlen)
  adjusted_pred_false <- ifelse(false.pred > 0.5, 1, false.pred)
  return(adjusted_pred_false / (1 + lambda * sqrt(cov)))
}
##### Initial Population #####
## Function to generate the initial population
initial_population <- function(onehot_instance, pop.dense) {
  # Randomly select windows of seg.len from the sequence
  samps <- sample(seq_len(maxlen - seg.len), seg.num, replace = FALSE)
  idx_matrix <- do.call(rbind, lapply(samps, function(s) s:(s + seg.len)))
  new_instance <- onehot_instance
  # Replace the selected windows with background noise
  new_instance[1, idx_matrix, ] <- matrix(pop.dense, 
                                          nrow = length(idx_matrix), 
                                          ncol = length(pop.dense), 
                                          byrow = TRUE)
    #matrix(rep(pop.dense, length(idx_matrix)),
    #                                      nrow = length(idx_matrix), byrow = TRUE)
  # also tried: outer(idx_matrix, pop.dense, FUN = function(i, j) j)
  pred <- predict(model, new_instance, verbose = 0)
  eff.preds <- eff_pred(pred[1, fidx], samps, samps + seg.len, maxlen)
  if (length(pred) > 1) {
    return(list(starts = samps,
                ends = samps + seg.len,
                eff.pred = eff.preds))
  } else {
    return(NULL)
  }
}

##### Population Evaluation #####
#### Function to evaluate the population
#### Step 1: extract eff.pred values
#### Step 2: rank the individuals based on eff.pred
#### Step 3: assign selection probs based on roulette wheel selection
evaluate_population <- function(population) {
  # Extract eff.pred from each individual in the population
  eff.preds <- unlist(lapply(population, function(res) if (!is.null(res)) res$eff.pred))
  ranked_idx <- order(eff.preds, decreasing = TRUE)
  ranked_population <- population[ranked_idx]
  # Calculate total sum of eff.preds
  total_eff.pred <- sum(eff.preds, na.rm = TRUE)
  eff.preds <- eff.preds[ranked_idx]
  # Assign selection probabilities proportional to eff.preds (Roulette Wheel)
  selection_probs <- eff.preds / total_eff.pred
  selection_probs[selection_probs < 0] <- 0
  
  result <- list(
    ranked_population = ranked_population,
    selection_probs = selection_probs
  )
  return(result)
}

##### Selection and Crossover #####
## Number of offspring: 75% of the population will be replaced with new individuals

perform_crossover <- function(parent1, parent2, maxlen) {
  # Sample a cutoff index
  cutoff <- sample(1:maxlen, 1)
  
  # Function to split into front and behind based on cutoff
  split_parent <- function(parent) {
    front_idx <- which(parent$ends <= cutoff)
    behind_idx <- which(parent$starts > cutoff)
    list(
      front = list(starts = parent$starts[front_idx], ends = parent$ends[front_idx]),
      behind = list(starts = parent$starts[behind_idx], ends = parent$ends[behind_idx])
    )
  }
  
  # Split both parents into front and behind parts
  parent1_parts <- split_parent(parent1)
  parent2_parts <- split_parent(parent2)
  
  # Combine front from one parent with behind from the other
  offspring1 <- list(
    starts = c(parent1_parts$front$starts, parent2_parts$behind$starts),
    ends = c(parent1_parts$front$ends, parent2_parts$behind$ends)
  )
  
  offspring2 <- list(
    starts = c(parent2_parts$front$starts, parent1_parts$behind$starts),
    ends = c(parent2_parts$front$ends, parent1_parts$behind$ends)
  )
  
  return(list(offspring1, offspring2))
}

generate_offspring <- function(population, selection_probs, update.rate = 0.75, maxlen) {
  offsprings <- list()
  num_offspring <- round(update.rate * iterations)
  num_crossovers <- floor(num_offspring / 2)  # Number of crossovers needed for pairs of offspring
  # Select parents and generate offspring pairs
  for (i in 1:num_crossovers) {
    parent_idx <- sample(seq_along(population), 2, prob = selection_probs, replace = TRUE)
    parent1 <- population[[parent_idx[1]]]
    parent2 <- population[[parent_idx[2]]]
    
    # Perform crossover to generate two offspring
    new_offspring <- perform_crossover(parent1, parent2, maxlen)
    
    # Append both offspring
    offsprings <- append(offsprings, new_offspring)
  }
  
  # If num_offspring is odd, generate one more offspring
  if (num_offspring %% 2 == 1) {
    parent_idx <- sample(seq_along(population), 2, prob = selection_probs, replace = TRUE)
    parent1 <- population[[parent_idx[1]]]
    parent2 <- population[[parent_idx[2]]]
    
    # Perform crossover but only take the first offspring
    new_offspring <- perform_crossover(parent1, parent2, maxlen)
    offsprings <- append(offsprings, new_offspring[1])  # Add only the first offspring
  }
  return(offsprings)
}

# Replace the least fit 75% of the population with the new offspring
replace_population <- function(population, offspring, update.rate = 0.75) {
  num_to_replace <- round(update.rate * iterations)
  ranked_population <- population[1:(length(population) - num_to_replace)]  # Retain the top (1 - 0.75) portion
  new_population <- append(ranked_population, offspring)  # Add new offspring
  return(new_population)
}

##### Mutation #####
## two types of mutations: elongation and deletion
## for each individual, each start:end pair has a probability of 0.2 of being mutated
## => sample from "no", "elongation", "deletion" with probabilities 0.8, 0.1, 0.1
## deletion: this pair is removed from the vectors
## elongation: the pair is extended by 100 bp in both directions, use max,min to control indexing errors
perform_mutation <- function(individual, maxlen, mutation_prob = 0.2) {
  # Iterate over each start:end pair and decide whether to mutate
  starts <- individual$starts
  ends <- individual$ends
  
  # Mutate each start:end pair with probabilities 0.8 (no mutation), 0.1 (elongation), 0.1 (deletion)
  mutation_choices <- sample(c("no", "elongation", "deletion"), 
                             size = length(starts), 
                             replace = TRUE, 
                             prob = c(0.8, 0.1, 0.1))
  
  new_starts <- c()
  new_ends <- c()
  
  for (i in seq_along(starts)) {
    mutation_type <- mutation_choices[i]
    
    if (mutation_type == "no") {
      # No mutation, keep the start:end pair as it is
      new_starts <- c(new_starts, starts[i])
      new_ends <- c(new_ends, ends[i])
      
    } else if (mutation_type == "elongation") {
      # Elongation: Extend by 100 bp in both directions
      # this can also happen reversely (shortening by 100 bp)
      if (runif(1) < 0.5) {
        new_start <- max(1, starts[i] - 1000)
        new_end <- min(maxlen, ends[i] + 1000)
      } else {
        new_start <- max(1, starts[i] + 1000)
        new_end <- min(maxlen, ends[i] - 1000)
      }
      new_starts <- c(new_starts, new_start)
      new_ends <- c(new_ends, new_end)
      
    } else if (mutation_type == "deletion") {
      # 50% chance of replacement by another start:end pair
      if (runif(1) < 0.5) {
        new_start <- sample(1:(maxlen - seg.len), 1)
        new_starts <- c(new_starts, new_start)
        new_ends <- c(new_ends, new_start + seg.len)
      } else {
        # Deletion: Do nothing, i.e., remove this start:end pair from the new vectors
        next
      }
    }
  }
  # each individual has a 0.1 probability of receiving an additional start:end pair (window)
  #if (runif(1) < 0.1) {
  #  new_start <- sample(1:(maxlen - seg.len), 1)
  #  new_starts <- c(new_starts, new_start)
  #  new_ends <- c(new_ends, new_start + seg.len)
  #}
  ### NB: Addition makes performance worse!
  
  # Update the individual with the new mutated start:end pairs
  individual$starts <- new_starts
  individual$ends <- new_ends
  
  return(individual)
}

# Apply mutation to the population
mutate_population <- function(population, maxlen) {
  mutated_population <- lapply(population, function(individual) {
    perform_mutation(individual, maxlen)
  })
  return(mutated_population)
}

##### Next Generation #####
next_population <- function(population, i, onehot_instance, pop.dense) {
  # Get starts and ends for the i-th individual
  samps <- population[[i]]$starts
  seg.len <- population[[i]]$ends - population[[i]]$starts
  # Build the index matrix using starts and ends
  idx_matrix <- unlist(mapply(seq, samps, population[[i]]$ends, SIMPLIFY = FALSE))
  # Create new instance and replace the selected windows with background noise
  new_instance <- onehot_instance
  new_instance[1, idx_matrix, ] <- matrix(pop.dense, 
                                          nrow = length(idx_matrix), 
                                          ncol = length(pop.dense), 
                                          byrow = TRUE) # hier tensorflow tensor, aber kopie explizit machen
    #matrix(rep(pop.dense, length(idx_matrix)),
    #                                      nrow = length(idx_matrix), byrow = TRUE)
  pred <- predict(model, new_instance, verbose = 0)
  # Calculate Eff.pred
  eff.preds <- eff_pred(pred[1, fidx], samps, population[[i]]$ends, maxlen)

  if (length(pred) > 1) {
    return(list(starts = samps,
                ends = population[[i]]$ends,
                eff.pred = eff.preds))
  } else {
    return(NULL)
  }
}

##### Stop Condition #####
# 1. N generations produced. or
# 2. best Eff.pred has not significantly improved in the last M generations
# Let M = 3, N = 50

##### Summary Function #####
print_setup <- function(maxlen, seg.len, seg.num, iterations, fidx, N, M,
                        update_rate, mutation_prob, pop.dense, lambda, pred.orignal) {
  cat("\n#### Genetic Algorithm Setup ####\n")
  cat("Max Length of Sequence (maxlen): ", maxlen, "\n")
  cat("Initial Segment Length:          ", seg.len, "\n")
  cat("Initial Number of Segments:      ", seg.num, "\n")
  cat("Iterations (Population Size):    ", iterations, "\n")
  cat("False Prediction Index:          ", fidx, "\n")
  cat("Original False Prediction:       ", "\n")
  cat("Coverage penalization (lambda):  ", lambda, "\n")
  cat("Stop Criteria:\n")
  cat("   - Max Generations:            ", N, "\n")
  cat("   - Stagnant Generations:       ", M, "\n")
  cat("Update Rate:                     ", update_rate * 100, "% of population replaced\n")
  cat("Elitism:                          1 individual\n")
  cat("Mutation Probability:            ", mutation_prob * 100, "% per start:end pair\n")
  cat("\nBackground Nucleotide Density:\n")
  cat("   A: ", round(pop.dense[1], 4), "\n")
  cat("   C: ", round(pop.dense[2], 4), "\n")
  cat("   G: ", round(pop.dense[3], 4), "\n")
  cat("   T: ", round(pop.dense[4], 4), "\n")
  cat("###################################\n")
}
genetically_optimized_masking <- function(instance, onehot_instance, pop.dense, maxlen,
                                         seg.len, seg.num, iterations, model,
                                         fidx, N = 10, M = 3, lambda = 1,
                                         update_rate = 0.75, mutation_prob = 0.2,
                                         print_summary = TRUE) {
  pop.dense <- nt_density(instance)
  pred.orignal <- predict(model, onehot_instance, verbose = 0)
  pred.orignal <- pred.orignal[1, fidx]
  if (print_summary) {
    print_setup(maxlen, seg.len, seg.num, iterations, fidx, N, M, update_rate,
                mutation_prob, pop.dense, lambda, pred.orignal)
  }
  ### return NULL and a warning if pred.orignal > 0.5
  if (pred.orignal > 0.5) {
    warning("Original false prediction is already above 0.5")
    return(NULL)
  }
  ### initial population
  population <- lapply(1:iterations, function(i) initial_population(onehot_instance, pop.dense))
  best_eff_pred <- pred.orignal
  best_population <- population
  stagnant_generations <- 0
  generation <- 1
  population.eval <- evaluate_population(population)
  while (generation <= N && stagnant_generations < M) {
    ### Selection and Crossover
    offspring_population <- generate_offspring(population.eval$ranked_population, 
                                               population.eval$selection_probs, 
                                               update_rate, maxlen)

    new_population <- replace_population(population.eval$ranked_population, 
                                         offspring_population, 
                                         update_rate)
    ### Mutation, Elitism = 1
    elite <- new_population[1]
    mutated_population <- mutate_population(new_population[-1], maxlen)
    mutated_population <- c(elite, mutated_population)
    next.population <- lapply(1:iterations, function(i) next_population(mutated_population, i, onehot_instance, pop.dense))
    population.eval <- evaluate_population(next.population)
    # Check the best Eff.pred of the current generation
    current_best_eff_pred <- max(unlist(lapply(population.eval$ranked_population, function(ind) if (!is.null(ind)) ind$eff.pred)))
    # Stop Condition Checks
    if (current_best_eff_pred > best_eff_pred) {
      best_eff_pred <- current_best_eff_pred
      best_population <- population.eval$ranked_population
      stagnant_generations <- 0  # Reset stagnant generation counter if improvement is seen
    } else {
      stagnant_generations <- stagnant_generations + 1  # Increment if no significant improvement
    }
    if (print_summary) {
      cat("Generation:", generation, "Best Eff.pred:", current_best_eff_pred, "\n")
      cat("Stagnant Generations:", stagnant_generations, "\n")
    }
    # Move to the next generation
    population <- next.population
    generation <- generation + 1
  }
  cat("Stopped after", generation - 1, "generations\n")
  cat("Best Eff.pred:", best_eff_pred, "\n")
  final_population <- best_population
  # calculate coverage for the best individual
  best_individual <- final_population[[1]]
  best_coverage <- coverage(best_individual$ends, best_individual$starts, maxlen)
  cat("Best Coverage:", best_coverage, "\n")
  cat("###################################\n")
  # Return the final population and the best Eff.pred after N generations or when convergence is reached
  list(instance = instance, best_individual = best_individual)
}

# Run the genetic algorithm
# result <- genetically_optimized_masking(instance, onehot_instance,
#                         maxlen, seg.len, seg.num, iterations, model, fidx)

##### Processing multiple sequences #####
gom_all_fasta <- function(target_df, folder_path, maxlen, seg.len, seg.num, iterations,
                          model, fidx, N = 10, M = 3, lambda = 1, update_rate = 0.75,
                          mutation_prob = 0.2, sample.count = NULL) {
  # Filter the files with ability_TRUE == 1
  files_to_process <- target_df %>% filter(ability_TRUE == 1)
  
  # if sample.count is not NULL, sample the files
  if (!is.null(sample.count)) {
    files_to_process <- files_to_process %>% sample_n(sample.count)
  }
  print_setup(maxlen, seg.len, seg.num, iterations, fidx, N, M, update_rate,
              mutation_prob, pop.dense=NA, lambda, pred.orignal=NA)
  # Apply the genetic optimization to each file
  results <- lapply(files_to_process$file, function(fasta_file) {
    fasta_path <- file.path(folder_path, fasta_file)
    cat("Processing file:", fasta_file, "\n")
    result <- gom_fasta(fasta_path, maxlen,seg.len, seg.num, iterations, model,
                        fidx, N, M, lambda, update_rate, mutation_prob)
    return(list(file_name = fasta_file, result = result))
  })
  return(results)
}

#### Function to process a single FASTA file
gom_fasta <- function(fasta_path, maxlen, seg.len, seg.num, iterations, model,
                      fidx, N = 10, M = 3, lambda = 1, update_rate = 0.75,
                      mutation_prob = 0.2) {
  instance <- microseq::readFasta(fasta_path)$Sequence[1]
  sequence_length <- nchar(instance)
  num_segments <- floor(sequence_length / maxlen)
  if (sequence_length < maxlen) {
    num_segments <- 1
    # if shorter than 0.5 * maxlen, return NULL and a warning
    if (sequence_length < 0.1 * maxlen) {
      cat("Sequence length is less than 10% of the maximum length \n")
      return(NULL)
    }
  }
  cat("Total sequence length:", sequence_length, "\n")
  cat("Number of full-length segments:", num_segments, "\n")
  segment_results <- lapply(1:num_segments, function(i) {
    # Extract each maxlen segment from the full sequence
    start_ind <- (i - 1) * maxlen + 1
    instance_sub <- substr(instance, start_ind, start_ind + maxlen - 1)
    if (length(instance_sub) < maxlen) {
      instance_sub <- paste0(instance_sub, paste(rep("N", maxlen - length(instance_sub)), collapse = ""))
    }
    # Create the one-hot encoding for this segment
    onehot_instance <- seq_encoding_label(char_sequence = instance_sub,
                                          maxlen = maxlen,
                                          start_ind = 1,
                                          vocabulary = c("A", "C", "G", "T"))
    pop.dense <- nt_density(instance_sub)
    # Apply genetically optimized masking to the segment
    result <- genetically_optimized_masking(instance_sub, onehot_instance, pop.dense, maxlen = maxlen, 
                                           seg.len = seg.len, seg.num = seg.num, iterations = iterations, 
                                           model = model, fidx = fidx, N = N, M = M, 
                                           lambda = lambda, update_rate = update_rate, 
                                           mutation_prob = mutation_prob, print_summary = FALSE)
    
    return(result)
  })

  return(segment_results)
}

