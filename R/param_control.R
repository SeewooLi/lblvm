build_loading_matrix <- function(formula_string, variable_order = NULL) {
  # Split the formula string into lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)  # remove leading/trailing whitespace
  lines <- lines[lines != ""]  # remove empty lines
  lines <- lines[grepl("^f\\w*\\s*~", lines)]

  # Initialize list to store mappings
  mapping <- list()

  for (line in lines) {
    # Split formula into LHS and RHS
    parts <- strsplit(line, "~")[[1]]
    factor <- trimws(parts[1])
    variables <- unlist(strsplit(parts[2], "\\+"))
    variables <- trimws(variables)

    for (var in variables) {
      mapping[[var]] <- c(mapping[[var]], factor)
    }
  }

  # Determine variable and factor order
  all_vars <- unique(names(mapping))
  all_factors <- sort(unique(unlist(mapping)))

  # Use user-specified order if given
  if (!is.null(variable_order)) {
    if (!all(all_vars %in% variable_order)) {
      stop("Some variables in the formula are not in the data.")
    }
    all_vars <- variable_order
  } else {
    all_vars <- sort(all_vars)
  }

  # Create matrix
  mat <- matrix(0, nrow = length(all_vars), ncol = length(all_factors),
                dimnames = list(all_vars, all_factors))

  # Fill in 1s where variables load onto factors
  for (var in names(mapping)) {
    for (fac in mapping[[var]]) {
      mat[var, fac] <- 1
    }
  }

  return(mat)
}
build_constraint_matrix <- function(formula_string, variable_order, factor_order) {
  # Split and trim lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)
  constraint_lines <- lines[grepl("==", lines)]

  # Initialize zero matrix
  mat <- matrix(0, nrow = length(variable_order), ncol = length(factor_order),
                dimnames = list(variable_order, factor_order))

  # Assign unique constraint group IDs
  for (i in seq_along(constraint_lines)) {
    terms <- trimws(unlist(strsplit(constraint_lines[i], "==")))
    for (term in terms) {
      split_term <- strsplit(term, "\\.")[[1]]
      var <- split_term[1]
      fac <- split_term[2]
      mat[var, fac] <- i
    }
  }

  return(mat)
}
apply_fixed_values <- function(formula_string, init_matrix, loading_matrix) {
  # Extract and trim assignment lines
  lines <- unlist(strsplit(formula_string, "\n"))
  lines <- trimws(lines)
  assign_lines <- lines[grepl("<-", lines, fixed = TRUE)]

  for (line in assign_lines) {
    parts <- unlist(strsplit(line, "<-", fixed = TRUE))
    term <- trimws(parts[1])  # e.g., "r4.f1"
    value <- as.numeric(trimws(parts[2]))  # e.g., 3

    split_term <- strsplit(term, "\\.")[[1]]
    var <- split_term[1]
    fac <- split_term[2]

    # Set value and mark as fixed
    if (var %in% rownames(init_matrix) && fac %in% colnames(init_matrix)) {
      init_matrix[var, fac] <- value
      loading_matrix[var, fac] <- 0
    } else {
      warning(sprintf("Invalid variable/factor name: %s.%s", var, fac))
    }
  }

  return(list(initial = init_matrix, loading = loading_matrix))
}

apply_equal_constraints_seq <- function(param_id_mat, eq_mat) {
  stopifnot(all(dim(param_id_mat) == dim(eq_mat)))

  rn <- rownames(param_id_mat)
  cn <- colnames(param_id_mat)

  # Convert to long form
  df <- as.data.frame(as.table(param_id_mat))
  colnames(df) <- c("row", "col", "param_id")
  df$eq <- as.vector(eq_mat)

  # Only keep cells with non-NA param_id
  df <- df[!is.na(df$param_id), ]

  # Track groups based on equal constraints
  edge_list <- do.call(rbind, lapply(split(df, df$eq), function(group) {
    ids <- unique(group$param_id)
    if (group$eq[1] != 0 && length(ids) > 1) t(utils::combn(ids, 2)) else NULL
  }))

  # Build graph of equal constraints
  if (!is.null(edge_list)) {
    g <- igraph::graph_from_edgelist(apply(edge_list, 2, as.character), directed = FALSE)
  } else {
    g <- igraph::make_empty_graph()
  }

  all_ids <- unique(as.character(df$param_id))
  g <- igraph::add_vertices(g, nv = length(setdiff(all_ids, igraph::V(g)$name)),
                            name = setdiff(all_ids, igraph::V(g)$name))

  comps <- igraph::components(g)
  grouped_ids <- split(as.integer(igraph::V(g)$name), comps$membership)

  # Assign sequential new IDs to grouped values
  id_map <- list()
  next_id <- 1
  for (group in grouped_ids) {
    for (id in group) {
      id_map[[as.character(id)]] <- next_id
    }
    next_id <- next_id + 1
  }

  # Assign new IDs to ungrouped param_ids (not in constraint graph)
  remaining_ids <- setdiff(as.character(df$param_id), names(id_map))
  for (id in remaining_ids) {
    id_map[[id]] <- next_id
    next_id <- next_id + 1
  }

  # Rebuild the output matrix
  out_mat <- param_id_mat
  for (i in seq_len(nrow(out_mat))) {
    for (j in seq_len(ncol(out_mat))) {
      id <- param_id_mat[i, j]
      if (!is.na(id)) {
        out_mat[i, j] <- id_map[[as.character(id)]]
      }
    }
  }

  rownames(out_mat) <- rn
  colnames(out_mat) <- cn
  return(out_mat)
}

group_parameters <- function(df) {
  # Flatten all values, get unique parameter IDs
  vals <- stats::na.omit(as.numeric(unlist(df)))
  parent <- stats::setNames(as.list(vals), vals)  # union-find parent map

  # Union-find helpers
  find <- function(x) {
    while (parent[[as.character(x)]] != x) {
      parent[[as.character(x)]] <- parent[[as.character(parent[[as.character(x)]])]]
      x <- parent[[as.character(x)]]
    }
    x
  }

  union <- function(x, y) {
    root_x <- find(x)
    root_y <- find(y)
    if (root_x != root_y) {
      parent[[as.character(root_y)]] <<- root_x
    }
  }

  # For each row, union all values
  apply(df, 1, function(row) {
    row_vals <- stats::na.omit(as.numeric(row))
    if (length(row_vals) > 1) {
      for (i in 2:length(row_vals)) {
        union(row_vals[1], row_vals[i])
      }
    }
  })

  # Group by connected component (root)
  groups <- list()
  for (v in vals) {
    root <- find(v)
    root_str <- as.character(root)
    if (!root_str %in% names(groups)) {
      groups[[root_str]] <- c()
    }
    groups[[root_str]] <- c(groups[[root_str]], v)
  }

  # Clean and sort
  groups <- lapply(groups, function(x) sort(unique(x)))
  names(groups) <- paste0("[[", seq_along(groups), "]]")
  return(groups)
}


efa_helper_matrices <- function(d, data, eq_interval = FALSE){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  formula_string <- paste(paste(paste0("\n f",1:d), paste(cns, collapse = " + "), sep = " ~ "), collapse = "\n")

  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  init_mat <- matrix(nrow = ncol(data), ncol = max(category)+1)
  init_mat[,1] <- 2
  for(i in 1:nrow(init_mat)){
    init_mat[i, 2:(category[i]+1)] <- 0
  }

  init_mat <- cbind(matrix(rep(c((d:1)/d), each=ncol(data)), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2),
                    init_mat)
  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  load_mat[,1:d][upper.tri(load_mat[,1:d])] <- 0
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))
  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  init_mat <- init_mat * load_mat
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index
  ))
}

efa_helper_matrices_cont <- function(d, data){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  formula_string <- paste(paste(paste0("\n f",1:d), paste(cns, collapse = " + "), sep = " ~ "), collapse = "\n")

  # initial matrix for thresholds
  init_mat <- cbind(matrix(rep(c((d:1)/d), each=ncol(data)), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE), center = TRUE, scale = TRUE)/2),
                    2)
  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  load_mat[,1:d][upper.tri(load_mat[,1:d])] <- 0
  fns <- c(colnames(load_mat), "c", "nu")
  load_mat <- cbind(load_mat,1,1)
  init_mat <- init_mat * load_mat
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index
  ))
}

cfa_helper_matrices <- function(formula_string, data, eq_interval = FALSE){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  d <- ncol(load_mat)
  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))

  init_mat2 <- matrix(nrow = ncol(data), ncol = max(category)+1)
  init_mat2[,1] <- 2
  for(i in 1:nrow(init_mat2)){
    init_mat2[i, 2:(category[i]+1)] <- 0
  }


  init_mat <- cbind(matrix(rep(1/d, each=ncol(data)*d), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2)
  )* cbind(load_mat, 1)
  init_mat <- cbind(init_mat, init_mat2)

  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index,
    d = d
  ))
}

cfa_helper_matrices_cont <- function(formula_string, data){
  if(is.null(colnames(data))){
    cns <- paste0("v", 1:ncol(data))
  } else {
    cns <- colnames(data)
  }

  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  d <- ncol(load_mat)
  # initial matrix for thresholds
  fns <- c(colnames(load_mat), "c", "nu")
  load_mat <- cbind(load_mat,1,1)

  init_mat <- cbind(matrix(rep(1/d, each=ncol(data)*d), nrow=ncol(data)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE), center = TRUE, scale = TRUE)/2),
                    2
  )* load_mat

  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index,
    d = d
  ))
}


group_helper_matrices <- function(formula_string, data, eq_interval = FALSE, ngroup, group_label=NULL){
  if(is.null(colnames(data))){
    cns0 <- paste0("v", 1:ncol(data))
  } else {
    cns0 <- colnames(data)
  }
  cns <- c()
  if(is.null(group_label)){
    for(i in 1:ngroup){
      cns <- append(cns, paste0(cns0, paste0("g", i)))
    }
  } else{
    for(i in 1:ngroup){
      cns <- append(cns, paste0(cns0, group_label[i]))
    }

  }

  load_mat <- build_loading_matrix(formula_string, variable_order = cns)
  d <- ncol(load_mat)

  # initial matrix for thresholds
  data <- reorder_mat(data)
  category <- apply(data, 2, max, na.rm = TRUE)
  category <- rep(category, ngroup)
  fns <- c(colnames(load_mat), "c", "nu", paste0("t", 1:max(category)))

  init_mat2 <- matrix(nrow = nrow(load_mat), ncol = max(category)+1)
  init_mat2[,1] <- 2
  for(i in 1:nrow(init_mat2)){
    init_mat2[i, 2:(category[i]+1)] <- 0
  }


  init_mat <- cbind(matrix(rep(1/d, each=nrow(load_mat)*d), nrow=nrow(load_mat)),
                    as.numeric(scale(colMeans(data, na.rm = TRUE)/category, center = TRUE, scale = TRUE)/2)
  )* cbind(load_mat, 1)
  init_mat <- cbind(init_mat, init_mat2)

  load_mat <- cbind(load_mat,1,1,1*(!is.na(init_mat[,-(1:(d+2))])))
  par_id <- matrix(1:length(as.vector(init_mat)), ncol = ncol(init_mat))
  dimnames(init_mat) <- list(cns, fns)
  dimnames(load_mat) <- list(cns, fns)
  dimnames(par_id) <- list(cns, fns)

  eq_constraint <- build_constraint_matrix(formula_string, cns, fns)

  fxd <- apply_fixed_values(formula_string,init_matrix = init_mat, load_mat)
  init_mat <- fxd$initial
  load_mat <- fxd$loading

  if(eq_interval == TRUE) load_mat[, (d+3):ncol(load_mat)] <- 0

  par_id <- par_id * load_mat
  par_id[par_id == 0] <- NA

  par_id <- apply_equal_constraints_seq(par_id, eq_constraint)

  grouping_index <- group_parameters(par_id)

  return(list(
    init_mat = init_mat,
    load_mat = load_mat,
    eq_constraint = eq_constraint,
    par_id = par_id,
    grouping_index = grouping_index,
    d = d
  ))
}
