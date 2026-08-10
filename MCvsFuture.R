#!/usr/bin/env Rscript

library(dplyr)
library(cluster)
library(parallel)
library(future)
library(future.apply)

## ------------------------------------------------------------
## Load data and objects
## ------------------------------------------------------------

load("ptest.RData")

## Erwartet werden insbesondere:
## cs_p.df
## kmeans.id
## ggf. k_max / nb_iter
##
## Falls kmeans.id NICHT in ptest.RData gespeichert ist,
## hier die Funktion einfügen.

nb.cores <- 18
kmax <- 4
nstart <- 25

## ------------------------------------------------------------
## Prepare data once
## ------------------------------------------------------------

cat("\nPreparing split.df ...\n")

t_split <- system.time({

  split.df <- split(
    cs_p.df,
    interaction(
      cs_p.df$case_ID,
      cs_p.df$id,
      drop = TRUE
    )
  )

})

print(t_split)

cat(
  "\nNumber of individual clustering problems:",
  length(split.df),
  "\n"
)

## ------------------------------------------------------------
## Helper function for timing
## ------------------------------------------------------------

run_benchmark <- function(label, expr) {

  cat("\n")
  cat("====================================================\n")
  cat(label, "\n")
  cat("====================================================\n")

  gc()

  t0 <- Sys.time()

  timing <- system.time({
    result <- eval.parent(substitute(expr))
  })

  t1 <- Sys.time()

  cat("\nElapsed wall-clock time:\n")
  print(t1 - t0)

  cat("\nsystem.time():\n")
  print(timing)

  cat("\nNumber of returned elements:", length(result), "\n")

  invisible(list(
    result = result,
    timing = timing,
    walltime = as.numeric(
      difftime(t1, t0, units = "secs")
    )
  ))
}

## ------------------------------------------------------------
## 1. mclapply
## ------------------------------------------------------------

bench_mclapply <- run_benchmark(
  "1. parallel::mclapply",
  parallel::mclapply(
    split.df,
    kmeans.id,
    kmax = kmax,
    nstart = nstart,
    mc.cores = nb.cores,
    mc.preschedule = TRUE
  )
)

## Optional: check that result can be combined
cluster.mcl <- dplyr::bind_rows(bench_mclapply$result)

rm(cluster.mcl)
gc()

## ------------------------------------------------------------
## 2. future_lapply: default scheduling
## ------------------------------------------------------------

future::plan(
  future::multisession,
  workers = nb.cores
)

bench_future_default <- run_benchmark(
  "2. future_lapply - default scheduling",
  future.apply::future_lapply(
    split.df,
    kmeans.id,
    kmax = kmax,
    nstart = nstart,
    future.seed = TRUE
  )
)

cluster.future.default <-
  dplyr::bind_rows(bench_future_default$result)

rm(cluster.future.default)
gc()

## ------------------------------------------------------------
## 3. future_lapply: explicit scheduling
## ------------------------------------------------------------

bench_future_sched <- run_benchmark(
  "3. future_lapply - future.scheduling = 5",
  future.apply::future_lapply(
    split.df,
    kmeans.id,
    kmax = kmax,
    nstart = nstart,
    future.seed = TRUE,
    future.scheduling = 5
  )
)

cluster.future.sched <-
  dplyr::bind_rows(bench_future_sched$result)

rm(cluster.future.sched)
gc()

## Stop multisession workers
future::plan(future::sequential)

## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------

benchmark.results <- data.frame(
  method = c(
    "mclapply",
    "future_default",
    "future_scheduling_5"
  ),
  elapsed_seconds = c(
    bench_mclapply$walltime,
    bench_future_default$walltime,
    bench_future_sched$walltime
  )
)

benchmark.results$elapsed_minutes <-
  benchmark.results$elapsed_seconds / 60

cat("\n")
cat("====================================================\n")
cat("SUMMARY\n")
cat("====================================================\n")

print(benchmark.results)

write.csv(
  benchmark.results,
  "parallel_benchmark_results.csv",
  row.names = FALSE
)

cat(
  "\nResults saved to parallel_benchmark_results.csv\n"
)
