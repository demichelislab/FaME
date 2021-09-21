# Installs or imports the required libraries
librarian::shelf(
  stringi,
  pipeR,
  data.table,
  dplyr,
  tidyr,
  purrr,
  juliangehring/HighSpeedStats
)

# Combines two binary matrices with a specific operator (defaults to '|')
combine_mats <- function (mm1, mm2, FUN = `|`) {
  res <- FUN(mm1, mm2)
  storage.mode(res) <- 'integer'
  res
}

matrix_to_vector <- function (mm) {
  dim(mm) <- NULL
  mm
}

# Construct a binary aberration matrix based on the aberration passed as input.
#   Undefined samples are taken in consideration
build_aberration_matrix <- function (dat, aberration_to_extract, na_aberration = 'nd') {
  mm                       <- dat[, -1]
  res                      <- mm == aberration_to_extract
  res[mm == na_aberration] <- NA
  rownames(res)            <- dat[[1]]
  storage.mode(res)        <- 'integer'
  res
}

# Computes the contingency tables for all the genes in the matrices passed in
#   in input. Expected format of the matrices is genes on the rows and samples
#   on the columns
compute_counts <- function (mat1, mat2, use_names = TRUE) {
  are_matrices_identical <- identical(mat1, mat2)
  mms <- list(
    m1 = mat1[[1]],
    m2 = mat2[[1]]
  )
  nns <- pipeline({
    setNames(mms, c('row', 'col'))
    map(rownames)
  })
  if (use_names) {
    if (any(map_lgl(nns, is.null))) {
      stop('ERROR: All the matrices must have rownames to proceed')
    }
  } else {
    nns <- pipeline({
      setNames(mms, c('row', 'col'))
      map(nrow)
      map(~ seq.int(1, .x))
    })
  }
  ll <- map(nns, length)
  def_elems <- pipeline({
    mms
    map(~ {
      xx               <- !is.na(.x)
      storage.mode(xx) <- 'integer'
      xx
    })
  })
  mms_na <- pipeline({
    mms
    map2(
      names(mms),
      function(mm, nn) pipeline({
        mm
        replace_na(0L)
        list(1L - .)
        setNames(paste0(nn, c('', '_n')))
        map(~ .x * def_elems[[nn]])
      })
    )
    flatten
  })
  res <- with(
    mms_na,
    data.table(
      g1   = rep(nns$row, ll$col),
      g2   = rep(nns$col, each = ll$row),
      n_11 = matrix_to_vector(tcrossprod(  m1,   m2)),
      n_10 = matrix_to_vector(tcrossprod(  m1, m2_n)),
      n_01 = matrix_to_vector(tcrossprod(m1_n,   m2)),
      n_00 = matrix_to_vector(tcrossprod(m1_n, m2_n))
    )
  )
  if (are_matrices_identical) {
    res <- if (use_names) {
      pipeline({
        ll$row
        diag
        upper.tri
        matrix_to_vector
        (res[.])
      })
    } else {
      res[g2 > g1]
    }
  }
  res
}

# Creates an empty results table with the correct types.
create_empty_dt <- function (use_names = TRUE) {
  gene_type = if (use_names) character else numeric
  data.table(
    g1         = gene_type(0),
    g2         = gene_type(0),
    p_value    = numeric(0),
    fdr        = numeric(0),
    odds_ratio = numeric(0),
    f_min      = numeric(0),
    n_11       = integer(0),
    n_10       = integer(0),
    n_01       = integer(0),
    n_00       = integer(0),
    n_tot      = integer(0),
    n_1        = integer(0),
    n_2        = integer(0),
    f_1        = numeric(0),
    f_2        = numeric(0),
    a1         = character(0),
    a2         = character(0)
  )
}

# Computes additional information to be added to the results table
add_other_info <- function (dat) {
  columns_order <- c(
    'g1',
    'g2',
    'p_value',
    'fdr',
    'odds_ratio',
    'f_min'
  )
  res <- dat
  if (nrow(res)) {
    res <- pipeline({
      dat[,
        `:=`(
          n_tot     = n_11 + n_10 + n_01 + n_00,
          n_1       = (n_11 + n_10),
          n_2       = (n_11 + n_01)
        )
      ][,
        `:=`(
          f_1       = n_1 / n_tot,
          f_2       = n_2 / n_tot
        )
      ][,
        `:=`(
          f_min      = pmin(f_1, f_2),
          odds_ratio = (n_11 / n_01) / (n_10 / n_00)
        )
      ][
        order(fdr, p_value, f_min)
      ]
      (~ gc())
    })
  } else {
    res <- create_empty_dt(use_names = is.character(dat$g1))
  }
  setcolorder(res, c(columns_order, setdiff(colnames(res), columns_order)))
  res
}

# Runs the fisher tests on the count data
run_tests <- function (dat, filter_freq = 0, fisher_fun = ultrafastfet) {
  pipeline({
    dat[
      (n_11 + pmin(n_10, n_01)) / (n_11 + n_10 + n_01 + n_00) > filter_freq
    ][,
      `:=`(
        p_value    = 1,
        fdr        = 1
      )
    ]
    (~ gc())
    (dd ~ {
      if (nrow(dd) > 0) {
        pipeline({
          dd[,
              p_value := fisher_fun(n_11, n_10, n_01, n_00)
          ]
          (~ gc())
          (.[,
              fdr := p.adjust(p_value, method = 'BH')
          ])
          (~ gc())
        })
      } else dd
    })
  })
}

# Runs the fisher tests on the count data
run_fame_on_matrix_pair <- function (
  mat1,
  mat2 = NULL,
  use_names = TRUE,
  filter_freq = 0,
  fisher_fun  = HighSpeedStats::ultrafastfet
) {
  mm1 <- list(m1 = mat1)
  mm2 <- mm1
  if (!is.null(mat2)) {
    mm2 <- list(m2 = mat2)
  }
  pipeline({
    compute_counts(mm1, mm2, use_names = use_names)
    run_tests(filter_freq = filter_freq, fisher_fun = fisher_fun)
    add_other_info()
  })
}

