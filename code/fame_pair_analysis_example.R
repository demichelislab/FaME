# Finds the path where to load the source file
orig_path <- getwd()
source_file <- head(grep('^-?-f(?:ile)=', commandArgs(trailingOnly = FALSE), value = TRUE, perl = TRUE), 1)
if (length(source_file)) {
  new_wd <- sub('^-[^=]+=(.*?)/?[^/]+$', '\\1', source_file)
  print(paste('Changing directory to', new_wd))
  setwd(new_wd)
} else {
  stop('Cannot determine source directory')
}
while (!dir.exists('code')) {
  setwd('..')
  if (getwd() == '/') {
    stop('Cannot find the root folder of the project')
  }
}
project_path <- getwd()
setwd(orig_path)
source(paste0(project_path, '/code/fame_core.R'))

nn <- 100
in1 <- matrix(sample(c(TRUE, FALSE, NA), nn^2, replace = TRUE, prob = c(10,10, 1)), nrow = nn)
in2 <- matrix(sample(c(TRUE, FALSE, NA), nn^2, replace = TRUE, prob = c(10,10, 1)), nrow = nn)

(results_with_two_different_matrices <- run_fame_on_matrix_pair(in1, in2, use_names = FALSE))
(results_with_a_single_matrix <- run_fame_on_matrix_pair(in1, use_names = FALSE))


