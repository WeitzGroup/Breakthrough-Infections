# R Package Installation Script
# Run this script once to install all required packages

# CRAN packages
cran_packages <- c(
  "tidyverse",      # Data manipulation and visualization (includes ggplot2, dplyr, tidyr, etc.)
  "readxl",         # Reading Excel files
  "magrittr",       # Pipe operators
  "cowplot",        # Plot composition
  "usmap",          # US map visualization
  "latex2exp",      # LaTeX expressions in plots
  "ggrepel",        # Text labels that avoid overlap
  "viridis",        # Color palettes
  "scales",         # Scale functions for ggplot2
  "ggridges",       # Ridge plots
  "knitr",          # Report generation
  "kableExtra",     # Enhanced tables
  "deSolve"         # ODE solvers
)

# Install CRAN packages
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Special packages that need specific installation

# ggh4x - May need specific version for compatibility
# If you get errors, try: remotes::install_version("ggh4x", version = "0.2.8")
install.packages("ggh4x", repos = "https://cloud.r-project.org")

# ggmagnify - Install from R-universe
install.packages(
  "ggmagnify",
  repos = c("https://hughjonesd.r-universe.dev", "https://cloud.r-project.org")
)

# Verify installation
cat("\n========================================\n")
cat("Checking installed packages...\n")
cat("========================================\n")

all_packages <- c(cran_packages, "ggh4x", "ggmagnify")
for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s\n", pkg))
  } else {
    cat(sprintf("  ✗ %s (FAILED)\n", pkg))
  }
}

cat("\nInstallation complete!\n")
