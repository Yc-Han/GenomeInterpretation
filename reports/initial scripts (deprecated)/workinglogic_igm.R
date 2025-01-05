baseline_seq <- array(onehot_baseline, dim=dim(onehot_instance))
baseline_seq <- tensorflow::tf$cast(baseline_seq, 
                                    dtype = "float32")

path_gradients <- compute_gradients(model = model, input_tensor = baseline_seq,
                  target_class_idx = 3, input_idx = NULL, pred_stepwise = FALSE)

avg_grads <- integral_approximation(gradients = path_gradients)

igw <- (onehot_instance - baseline_seq)[1,,] * avg_grads[,]

abs_sum <- rowSums(abs(as.array(igw)))
df <- data.frame(abs_sum = abs_sum, position = 1 : 477)
ggplot(df, aes(x = position, y = abs_sum)) + geom_rect(aes(xmin = 90, xmax = 180, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.2) + geom_smooth(method = "auto") + geom_point() + theme_bw()

heatmaps_integrated_grad(integrated_grads = igw,
                         input_seq = onehot_instance)

####
py_run_string("import tensorflow as tf")
py$gradients <- path_gradients
# riemann_trapezoidal
py_run_string("
grads_mid = (gradients[:, :-1, :] + gradients[:, 1:, :]) / tf.constant(2.0)
grads_with_endpoints = tf.concat([grads_mid, gradients[:, -1:, :]], axis=1)
")
py_run_string("integrated_gradients = tf.math.reduce_mean(grads_with_endpoints, axis=0)")


####
py_run_string("import tensorflow as tf")
py$input_tensor <- baseline_seq
py$input_idx <- as.integer(NULL - 1)
py$target_class_idx <- as.integer(3 - 1)
py$model <- model

if (!is.null(input_idx)) {
  py_run_string(
    "with tf.GradientTape() as tape:
             tape.watch(input_tensor[input_idx])
             probs = model(input_tensor)[:, target_class_idx]
    ")
} else {
  py_run_string(
    "with tf.GradientTape() as tape:
             tape.watch(input_tensor)
             probs = model(input_tensor)[:, target_class_idx]
    ")
}

grad <- py$tape$gradient(py$probs, py$input_tensor)
if (!is.null(input_idx)) {
  return(grad[[input_idx]])
} else {
  return(grad)
}
