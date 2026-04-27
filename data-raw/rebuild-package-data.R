#!/usr/bin/env Rscript

dataset_files <- Sys.glob("data/*.rdata")

if (!length(dataset_files)) {
  stop("No .rdata files found in data/.")
}

for (path in dataset_files) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  object_names <- ls(env, all.names = TRUE)

  if (length(object_names) != 1L) {
    stop("Expected exactly one object in ", path, ".")
  }

  object_name <- object_names[[1]]
  object_value <- tibble::as_tibble(env[[object_name]])
  assign(object_name, object_value)

  output_path <- file.path("data", paste0(object_name, ".rda"))
  save(list = object_name, file = output_path, compress = "bzip2", version = 2)

  rm(list = object_name)
}

unlink(dataset_files)
