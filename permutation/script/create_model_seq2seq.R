library(keras)
library(reticulate)
library(tensorflow)
library(bbotk)
library(mlrintermbo)
library(paradox)

trial <- 4
batch_size <- 16

# Ensure the checkpoint directory exists
checkpoint_dir <- "permutation/seq2seq_cps"
if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

create_model_seq2seq <- function(
    ### THE CURRENT BEST SETTING IS SAVED AS DEFAULT HERE
  maxlen = 5000,
  vocab_size = 4,
  
  # Convolutional layer parameters
  conv_filters = c(32, 64),
  conv_kernels = c(16, 8),
  leaky_relu_alpha = 0.85,
  
  # Residual block parameters (2 convolutional layers inside the block)
  residual_filters = 32, 
  residual_kernel_size = 4, 
  
  # LSTM parameters
  lstm_units = 64,
  lstm_layers = 1,
  dropout_lstm = 0.30,
  recurrent_dropout_lstm = 0.2,
  
  # Final output layer
  final_units = 1,
  final_activation = "sigmoid",
  
  # Regularization
  batch_norm_momentum = 0.1288,
  kernel_regularizer_rate = 0,
  dropout_rate = 0.3,
  
  # Training-related parameters
  learning_rate = 0.001,
  loss_fn = balanced_cross_entropy(),
  optimizer_type = "adam"
) {
  
  input <- layer_input(shape = c(maxlen, vocab_size))
  
  # Add multiple convolutional layers dynamically based on conv_filters and conv_kernels
  x <- input
  for (i in seq_along(conv_filters)) {
    x <- x %>%
      layer_conv_1d(filters = conv_filters[i], kernel_size = conv_kernels[i], padding = "same",
                    kernel_regularizer = regularizer_l2(kernel_regularizer_rate)) %>%
      layer_batch_normalization(momentum = batch_norm_momentum) %>%
      layer_activation_leaky_relu(alpha = leaky_relu_alpha)
  }
  
  # Residual block with regularization and dropout
  res <- x
  
  x <- x %>%
    layer_conv_1d(filters = residual_filters, kernel_size = residual_kernel_size, padding = "same",
                  kernel_regularizer = regularizer_l2(kernel_regularizer_rate)) %>%
    layer_batch_normalization(momentum = batch_norm_momentum) %>%
    layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
    layer_conv_1d(filters = residual_filters, kernel_size = residual_kernel_size, padding = "same",
                  kernel_regularizer = regularizer_l2(kernel_regularizer_rate)) %>%
    layer_batch_normalization(momentum = batch_norm_momentum)
  
  if (residual_filters != conv_filters[length(conv_filters)]) {
    res <- res %>% 
      layer_conv_1d(filters = residual_filters, kernel_size = 1, padding = "same",
                    kernel_regularizer = regularizer_l2(kernel_regularizer_rate))
  }
  
  x <- layer_add(list(x, res)) %>%
    layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
    layer_dropout(rate = dropout_rate)
  
  # LSTM layers with dropout and L2 regularization
  for (i in seq_len(lstm_layers)) {
    x <- x %>%
      bidirectional(layer_lstm(units = lstm_units, 
                               return_sequences = TRUE, 
                               dropout = dropout_lstm, 
                               recurrent_dropout = recurrent_dropout_lstm,
                               kernel_regularizer = regularizer_l2(kernel_regularizer_rate)))
  }
  
  # TimeDistributed Dense for per-timestep output
  x <- x %>%
    time_distributed(layer_dense(units = final_units, activation = final_activation,
                                 kernel_regularizer = regularizer_l2(kernel_regularizer_rate)))
  
  model <- keras_model(inputs = input, outputs = x)
  
  # Choose optimizer
  opt <- switch(
    optimizer_type,
    "adam" = optimizer_adam(learning_rate = learning_rate),
    "rmsprop" = optimizer_rmsprop(learning_rate = learning_rate),
    "sgd" = optimizer_sgd(learning_rate = learning_rate),
    stop("Unsupported optimizer type. Choose 'adam', 'rmsprop', or 'sgd'.")
  )
  
  model %>% compile(
    optimizer = opt,
    loss = loss_fn
  )
  
  return(model)
}

balanced_cross_entropy <- function() {
  function(y_true, y_pred) {
    epsilon <- k_epsilon()  # Small constant to avoid log(0)
    
    # Clip predictions to avoid log(0) errors
    y_pred <- k_clip(y_pred, epsilon, 1 - epsilon)
    
    # Calculate the class weights dynamically
    pos_weight <- k_sum(1 - y_true) / k_sum(y_true)
    neg_weight <- 1
    
    # Calculate the weighted losses
    pos_loss <- -pos_weight * y_true * k_log(y_pred)
    neg_loss <- -neg_weight * (1 - y_true) * k_log(1 - y_pred)
    
    # Combine losses
    loss <- pos_loss + neg_loss
    return(k_mean(loss))
  }
}

model <- create_model_seq2seq(
conv_filters = c(60, 120, 240),
  conv_kernels = c(6, 5, 4),
  residual_filters = 64,
  residual_kernel_size = 4,
  lstm_units = 120,
  lstm_layers = 1,
  dropout_lstm = 0.30,
  recurrent_dropout_lstm = 0.2,
  final_units = 1,
  final_activation = "sigmoid",
  batch_norm_momentum = 0.15,
  kernel_regularizer_rate = 0.0001,
  learning_rate = 0.0005,
  loss_fn = balanced_cross_entropy(),
  optimizer_type = "adam"
)
summary(model)

#model %>% load_model_weights_hdf5("permutation/seq2seq_cps/trial2_best_model_epoch_04_val_loss_1.34.h5")

set.seed(42)

# List all grouped files
file_paths <- list.files("permutation/tensors/", pattern = "\\.rds$", full.names = TRUE)

# Randomly shuffle and split into train and validation sets
file_paths <- sample(file_paths)
val_fraction <- 0.1
val_size <- ceiling(length(file_paths) * val_fraction)

val_file_paths <- file_paths[1:val_size]
train_file_paths <- file_paths[(val_size + 1):length(file_paths)]

#### For cluster slow loading
data_generator <- function(file_paths, batch_size = 32, maxlen = 5000, vocab_size = 4) {
  iterator <- list(
    file_index = 1,
    batch_index = 1,
    current_data = NULL,
    num_batches = 0,
    inputs = NULL,
    masks = NULL
  )
  
  reset_iterator <- function() {
    iterator$file_index <- 1
    iterator$batch_index <- 1
    iterator$current_data <- NULL
    iterator$num_batches <- 0
    iterator$inputs <- NULL
    iterator$masks <- NULL
  }
  
  next_batch <- function() {
    # Generator is called
    if (is.null(iterator$current_data) || iterator$batch_index > iterator$num_batches) {
      if (iterator$file_index > length(file_paths)) {
        cat("Returning NULL from generator.\n")
        reset_iterator()
        return(NULL)
      }
      # Load new data
      iterator$current_data <- readRDS(file_paths[iterator$file_index])
      inputs <- iterator$current_data$inputs
      masks <- iterator$current_data$masks
      
      iterator$num_batches <- ceiling(length(inputs) / batch_size)
      iterator$inputs <- inputs
      iterator$masks <- masks
      iterator$batch_index <- 1
      iterator$file_index <- iterator$file_index %% length(file_paths) + 1  # Cycle through files
    }
    
    start_ind <- (iterator$batch_index - 1) * batch_size + 1
    end_ind <- min(iterator$batch_index * batch_size, length(iterator$inputs))
    
    if (start_ind > length(iterator$inputs)) {
      return(NULL)  # No more data for this epoch
    }
    
    input_batch <- iterator$inputs[start_ind:end_ind]
    mask_batch <- iterator$masks[start_ind:end_ind]
    iterator$batch_index <- iterator$batch_index + 1
    
    input_tensor <- array_reshape(do.call(rbind, input_batch), dim = c(length(input_batch), maxlen, vocab_size))
    mask_tensor <- array_reshape(do.call(rbind, mask_batch), dim = c(length(mask_batch), maxlen, 1))
    
    list(input_tensor, mask_tensor)
  }
  
  py_iterator(next_batch)
}

# Create training and validation generators
# train_generator <- data_generator(train_file_paths, batch_size = batch_size)
train_generator <- data_generator(train_file_paths, batch_size = batch_size)

ResettableDataGenerator <- function(generator_function) {
  env <- new.env()
  env$generator <- generator_function()
  # Reset function
  reset_generator <- function() {
    env$generator <- generator_function()
  }
  
  # Iterator wrapper
  next_batch <- function() {
    tryCatch(
      iter_next(env$generator),
      error = function(e) {
        if (inherits(e, "stop_iteration")) {
          reset_generator()
          return(env$generator())
        } else {
          stop(e)
        }
      }
    )
  }
  
  return(list(
    reset = reset_generator,
    `next` = next_batch
  ))
}

# Wrap your validation generator
val_generator <- ResettableDataGenerator(function() data_generator(val_file_paths, batch_size = batch_size))

model_path <- file.path(checkpoint_dir, "model_test.keras")
total_val_samples <- sum(sapply(val_file_paths, function(file) length(readRDS(file)$inputs)))
val_steps <- floor(total_val_samples / batch_size)

# Define the checkpoint callback
checkpoint_callback <- callback_model_checkpoint(
  filepath = file.path(checkpoint_dir, "trial_", trial, "_model_epoch_{epoch:02d}_val_loss_{val_loss:.2f}.h5"),
  save_weights_only = TRUE,
  save_best_only = TRUE,
  monitor = "val_loss",
  mode = "min",
  verbose = 1
)

callback_reduce_lr <- callback_reduce_lr_on_plateau(
  monitor = "val_loss",
  factor = 0.5,
  patience = 2,
  verbose = 1
)

#### Bayesian HPO

# define feasible values
domain <- ps(
  residual_filters = p_int(lower = 1),
  batch_norm_momentum = p_dbl(0, 1),
  dropout_rate = p_dbl(0, 1)
)

# define search space
# values are to_tune('lower bound', 'upper bound', [other options])
# if lower / upper not given: search entire space.
domain$values = list(
  residual_filters = to_tune(16, 256, logscale = TRUE),
  batch_norm_momentum = to_tune(0.1, 0.9),
  dropout_rate = to_tune(0, 0.5)
)

obj <- ObjectiveRFun$new(
  fun = function(xs) {
    model <- create_model_seq2seq(
      residual_filters = xs$residual_filters,
      dropout_rate = xs$dropout_rate,
      batch_norm_momentum = xs$batch_norm_momentum
    )
    history <- model %>% fit(
      x = train_generator,
      steps_per_epoch = 200,
      validation_data = val_generator$`next`,
      validation_steps = val_steps,
      epochs = 1,
      callbacks = list(
        callback_lambda(on_epoch_begin = function(epoch, logs) val_generator$reset())
      ),
      verbose = 0
    )
    val_loss <- tail(history$metrics$val_loss, 1)
    list(y = val_loss)
  },
  domain = domain,
  codomain = ps(y = p_dbl(tags = "minimize"))
)

number.of.evals <- 100

optimizer <- opt("intermbo", on.surrogate.error = "warn", infill.opt = "focussearch", infill.opt.focussearch.maxit = 20)

oi <- suppressMessages(OptimInstanceSingleCrit$new(obj, terminator = trm("evals", n_evals = number.of.evals)))

optimizer$optimize(oi)

oi$archive$data

