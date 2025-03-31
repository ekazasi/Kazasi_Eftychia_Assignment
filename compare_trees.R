#!/usr/bin/env Rscript

# Load necessary library
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape", repos = "https://cloud.r-project.org/")
}

library(ape)

# Function to display usage
usage <- function() {
  cat("Usage: compare_trees.R <tree_file1> <tree_file2>\n", file=stderr())
  quit(save = "no", status = 1)
}

# Check for command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  usage()
}

# Read tree files from command-line arguments
tree_file1 <- args[1]
tree_file2 <- args[2]

# Try reading the trees
#cat("Reading trees...\n")
tree1 <- tryCatch(read.tree(tree_file1), error = function(e) {
  cat("Error reading tree file 1:", e$message, "\n", file=stderr())
  quit(save = "no", status = 1)
})

tree2 <- tryCatch(read.tree(tree_file2), error = function(e) {
  cat("Error reading tree file 2:", e$message, "\n", file=stderr())
  quit(save = "no", status = 1)
})

# Check if both trees are rooted
#cat("Checking tree properties...\n")
#if (!is.rooted(tree1) || !is.rooted(tree2)) {
  #cat("Warning: One or both trees are not rooted. Results may be less reliable.\n")
#}

# Compare trees: Robinson-Foulds distance
#cat("Calculating Robinson-Foulds distance...\n")
rf_distance <- dist.topo(tree1, tree2)

# Topological equivalence
#cat("Checking topological equivalence...\n")
#are_equivalent <- all.equal.phylo(tree1, tree2, use.edge.length = FALSE)

# Output results
#cat("\nComparison Results:\n")
cat(paste("RF_", tree_file2, "\t",rf_distance, "\n", sep=""))
#cat("- Topological Equivalence:", ifelse(isTRUE(are_equivalent), "Yes", "No"), "\n")

#cat("Comparison complete.\n")
