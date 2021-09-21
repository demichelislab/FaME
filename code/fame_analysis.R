librarian::shelf(
  parallel
)

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

num_cores <- detectCores()

# Reads all the files with genomic data from the specified folder
data_per_gene <- pipeline({
  list.files('data/genomic', full.names = TRUE)
  map_dfr(fread)
  group_by(sample_id)
  filter(any(as_cn_disc != 'nd'))
  ungroup
  as.data.table
})

# All the possible combinations of these genes will be tested
interest_genes <- pipeline({
  fread('data/resources/gene_sets/cancer_and_druggable_genes.tsv.gz')
  pull(gene)
})

# The genomic data of the interest genes
interest_genes_data <- pipeline({
  data_per_gene[hugo %in% interest_genes]
  nest(data = -dataset)
  arrange(stri_order(dataset))
  (setNames(.$data, .$dataset))
  map( ~ pipeline({
      .x
      arrange(hugo, sample_id)
      as.data.table
    })
  )
})

# This matrix contains the allele specific data of the loaded datasets
as_cn_matrix <- pipeline({
  interest_genes_data
  map(~ dcast(.x, hugo ~ sample_id, value.var = 'as_cn_disc', fill = NA_character_))
})

# This builds the binary aberration matrix for a specific allele specific copy
#   number aberration. This matrix will have a value of 1 if a sample have
#   an Homozygous deletion at a specific gene 0 otherwise.
homo_del_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'homo_del'))
})

# This matrix will have a value of 1 if a sample have an Hemizygous deletion at
#   a specific gene 0 otherwise.
hemi_del_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'hemi_del'))
})

# This matrix will have a value of 1 if a sample have an Copy Neutral LOH at a
# specific gene 0 otherwise.
cnnl_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'cnnl'))
})

# This matrix will have a value of 1 if a sample have an amplification at a
# specific gene 0 otherwise.
amp_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'amp'))
})

# This matrix will have a value of 1 if a sample have an SNV at a
# specific gene 0 otherwise.
snv_matrix <- pipeline({
  interest_genes_data
  map(~ {
    pipeline({
      dcast(.x, hugo ~ sample_id, value.var = 'count_snvs_deleterious', fill = 0L, fun.aggregate = function (x) as.integer(sign(x)))
    })
  })
  map(~ as.matrix(.x[, -1], rownames.value = .x[[1]]))
})

all_datasets <- pipeline({
  data_per_gene
  (dataset)
  unique
  stri_sort
  setNames(.)
})

# Named list with the simple aberrations to test
aberrations_bases <- list(
  'hemi_del' = hemi_del_matrix,
  'cnnl'     = cnnl_matrix,
  'snv'      = snv_matrix
)

# Named list with the combined aberrations to test. In this case the aberrations
#   are combined with an OR operation.
aberration_combinations <- pipeline({
  list(
    hemi_cnnl     = c('hemi_del', 'cnnl'),
    hemi_cnnl_snv = c('hemi_del', 'cnnl', 'snv')
  )
  map(~ {
    selected_abs <- aberrations_bases[.x]
    map(all_datasets, function (nn) {
      cc <- map(selected_abs, ~ .x[[nn]])
      reduce(cc, combine_mats)
    })
  })
})

aberrations <- c(
  aberrations_bases,
  aberration_combinations
)

# This combines all the matrices for all the dataset in an unique matrix in
#   order to test all the dataset in a pan-cancer fashion.
pancancer_aberrations <- pipeline({
  aberrations
  map(~ list(pancancer = reduce(.x, cbind)))
})

# Generates all the possible combinations of aberrations. FAME efficiency enable
#   the possibility to test many aberrations combinations.
aberration_comparisons <- pipeline({
  c(
    'hemi_del',
    'cnnl',
    'hemi_cnnl',
    'hemi_cnnl_snv',
    'snv'
  )
  (expand.grid(
    a2               = .,
    a1               = .,
    stringsAsFactors = FALSE
  ))
  filter(
    a1 <= a2
  )
  select(a1, a2)
  as.data.table
})

# Tests all the aberrations combinations on all the pairs of genes separately
#   for each dataset.
results_per_project <- pipeline({
  aberration_comparisons
  pmap(c)
  map_dfr(~ {
    cat(.x, '\n')
    pipeline({
      aberrations[.x]
      map(~ { .x[all_datasets] })
      pmap(list)
      map(~ compute_counts(.x[1], .x[2]))
      mcmapply(names(.), FUN = function (dd, nn) {
        print(nn)
        res <- run_tests(dd)
        gc()
        res
      }, mc.cores = num_cores, SIMPLIFY = FALSE)
      bind_rows(.id = 'dataset')
      mutate(
        a1 = .x[[1]],
        a2 = .x[[2]]
      )
      select(a1, a2, everything())
      as.data.table
    })
  })
  (? gc())
})

fwrite(results_per_project, 'data/result_pairs.tsv', sep = '\t')

# Tests all the aberrations combinations on all the pairs of genes on the
#   pancancer dataset (all the datasets combined).
results_pancancer <- pipeline({
  aberration_comparisons
  pmap(c)
  map_dfr(~ {
    cat(.x, '\n')
    pipeline({
      pancancer_aberrations[.x]
      pmap(list)
      map(~ compute_counts(.x[1], .x[2]))
      mcmapply(names(.), FUN = function (dd, nn) {
        print(nn)
        res <- run_tests(dd)
        gc()
        res
      }, mc.cores = num_cores, SIMPLIFY = FALSE)
      bind_rows(.id = 'dataset')
      mutate(
        a1 = .x[[1]],
        a2 = .x[[2]]
      )
      select(a1, a2, everything())
      as.data.table
      (? gc())
    })
  })
  (? gc())
})

fwrite(results_pancancer, 'data/result_pairs_pancancer.tsv', sep = '\t')
