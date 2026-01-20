# Two-Dose Vaccination Model
# Converted from TwoDoseModel.ipynb (Python) to R
# This script simulates a two-dose vaccination model and analyzes breakthrough infections

# ============================================================================
# Set working directory to scripts/R folder
# ============================================================================
find_project_root <- function() {
  current <- getwd()
  while (current != dirname(current)) {
    if (file.exists(file.path(current, "README.md")) && 
        dir.exists(file.path(current, "scripts", "R"))) {
      return(current)
    }
    current <- dirname(current)
  }
  return(NULL)
}

project_root <- find_project_root()
if (!is.null(project_root)) {
  setwd(file.path(project_root, "scripts", "R"))
  cat("Working directory set to:", getwd(), "\n")
}

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scales)


# Create output directory if it doesn't exist
if (!dir.exists("../../output/figures")) {
  dir.create("../../output/figures", recursive = TRUE)
}

# =============================================================================
# Two-Dose ODE Model
# Compartments: S, S1, S2, I, B1, B2, R
# S  = Susceptible (unvaccinated)
# S1 = Susceptible after 1st dose (vaccine failure)
# S2 = Susceptible after 2nd dose (vaccine failure)
# I  = Infected (unvaccinated)
# B1 = Breakthrough infection (1st dose only)
# B2 = Breakthrough infection (2nd dose)
# R  = Recovered/Immune
# =============================================================================

Two_dose <- function(t, y, pars) {
  S  <- y[1]
  S1 <- y[2]
  S2 <- y[3]
  I  <- y[4]
  B1 <- y[5]
  B2 <- y[6]
  R  <- y[7]
  
  N <- S + S1 + S2 + I + B1 + B2 + R
  I_tot <- I + B1 + B2
  
  dS  <- pars$b - (pars$beta * S / N) * I_tot - pars$v1 * pars$p1 * S - pars$m * S
  dS1 <- pars$v1 * pars$p1 * pars$eps1 * S - (pars$beta * S1 / N) * I_tot - pars$v2 * pars$p2 * S1 - pars$m * S1
  dS2 <- pars$v2 * pars$p2 * pars$eps2 * S1 - (pars$beta * S2 / N) * I_tot - pars$m * S2
  dI  <- (pars$beta * S / N) * I_tot - pars$gamma * I - pars$m * I
  dB1 <- (pars$beta * S1 / N) * I_tot - pars$gamma * B1 - pars$m * B1
  dB2 <- (pars$beta * S2 / N) * I_tot - pars$gamma * B2 - pars$m * B2
  dR  <- pars$v1 * pars$p1 * (1 - pars$eps1) * S + pars$v2 * pars$p2 * (1 - pars$eps2) * S1 + pars$gamma * I_tot - pars$m * R
  
  list(c(dS, dS1, dS2, dI, dB1, dB2, dR))
}

# =============================================================================
# Main Simulation: Baseline scenario
# =============================================================================

p1_vals <- seq(0, 1, length.out = 100)
N <- 1e5

# Storage for results
results <- data.frame(
  p1 = numeric(),
  fv1 = numeric(),
  fv2 = numeric(),
  fv3 = numeric(),
  rv1 = numeric(),
  rv2 = numeric(),
  rv3 = numeric()
)

cat("Running baseline simulations...\n")

for (i in seq_along(p1_vals)) {
  # Model parameters
  pars <- list(
    eps1 = 0.07,
    eps2 = 0.03 / 0.07,  # conditional failure rate
    gamma = 1/10,
    beta = 1.515,
    p1 = p1_vals[i],
    p2 = p1_vals[i],
    v1 = 1 / (365/2),      # rate of first dose (6 months)
    v2 = 1 / (365 * 3.5),  # rate of second dose (3.5 years)
    m = 1 / (365 * 4.5),   # mortality/birth rate
    b = 1 / (365 * 4.5)    # birth rate = mortality for constant pop
  )
  
  # Time span: 50 years
  times <- seq(0, 365 * 50, by = 1)
  
  # Initial conditions (normalized)
  y0 <- c(
    S = (N - 1) / N,
    S1 = 0,
    S2 = 0,
    I = 1 / N,
    B1 = 0,
    B2 = 0,
    R = 0
  )
  
  # Solve ODE
  out <- ode(y = y0, times = times, func = Two_dose, parms = pars, 
             method = "lsoda", atol = 1e-8, rtol = 1e-8)
  
  # Get final (equilibrium) values
  final <- out[nrow(out), ]
  S_f  <- final["S"]
  S1_f <- final["S1"]
  S2_f <- final["S2"]
  I_f  <- final["I"]
  B1_f <- final["B1"]
  B2_f <- final["B2"]
  R_f  <- final["R"]
  
  I_tot_f <- I_f + B1_f + B2_f
  
  # Calculate breakthrough fractions
  fv1 <- ifelse(I_tot_f > 0, B1_f / I_tot_f, 0)
  fv2 <- ifelse(I_tot_f > 0, B2_f / I_tot_f, 0)
  fv3 <- ifelse(I_tot_f > 0, I_f / I_tot_f, 0)
  
  results <- rbind(results, data.frame(
    p1 = p1_vals[i],
    fv1 = as.numeric(fv1),
    fv2 = as.numeric(fv2),
    fv3 = as.numeric(fv3),
    rv1 = as.numeric(B1_f),
    rv2 = as.numeric(B2_f),
    rv3 = as.numeric(I_f)
  ))
  
  if (i %% 20 == 0) cat(sprintf("  Progress: %d%%\n", i))
}

cat("Baseline simulations complete.\n")

# =============================================================================
# Plot 1: Main figure (2-panel)
# =============================================================================

# Panel A: Breakthrough fractions
p1 <- ggplot(results) +
  geom_line(aes(x = p1, y = fv1, color = "First-dose vaccinated"), linewidth = 1) +
  geom_line(aes(x = p1, y = fv2, color = "Second-dose vaccinated"), linewidth = 1) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1e-6),
    breaks = c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = scales::scientific
  ) +
  scale_color_manual(values = c("First-dose vaccinated" = "red", 
                                 "Second-dose vaccinated" = "blue"),
                     labels = c(expression("First-dose vaccinated" (frac(I[V1]^"*", I[U]^"*" + I[V1]^"*" + I[V2]^"*"))),
                               expression("Second-dose vaccinated" (frac(I[V2]^"*", I[U]^"*" + I[V1]^"*" + I[V2]^"*"))))) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = "Breakthrough Fraction (log-scaled)",
       color = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.7, 0.3),
        legend.background = element_blank())

# Panel B: Infection counts at equilibrium
p2 <- ggplot(results) +
  geom_line(aes(x = p1, y = rv1 * N, color = "First-dose vaccinated"), linewidth = 1) +
  geom_line(aes(x = p1, y = rv2 * N, color = "Second-dose vaccinated"), linewidth = 1) +
  geom_line(aes(x = p1, y = rv3 * N, color = "Unvaccinated"), linewidth = 1) +
  scale_y_log10() +
  scale_color_manual(values = c("First-dose vaccinated" = "red", 
                                 "Second-dose vaccinated" = "blue",
                                 "Unvaccinated" = "black"),
                     labels = c(expression("First-dose vaccinated"  (I[V1]^"*")),
                               expression("Second-dose vaccinated"  (I[V2]^"*")),
                               expression("Unvaccinated" (I[U]^"*")))) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = "Infections at Equilibrium (log-scaled)",
       color = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.7, 0.3),
        legend.background = element_blank())

# Combine panels
pdf("../../output/figures/TwoDose_Breakthroughs.pdf", height = 6, width = 14)
print(plot_grid(p1, p2, nrow = 1, labels = "AUTO", label_size = 20))
dev.off()

cat("Saved: ../../output/figures/TwoDose_Breakthroughs.pdf\n")

# =============================================================================
# Sensitivity Analysis: Early second dose with higher failure rate
# =============================================================================

cat("\nRunning early vaccination scenario...\n")

a <- 3    # Second dose given 3x faster
b <- 1.5  # Conditional failure rate 1.5x higher

results_early <- data.frame(
  p1 = numeric(),
  fv1_early = numeric(),
  fv2_early = numeric(),
  rv1_early = numeric(),
  rv2_early = numeric()
)

for (i in seq_along(p1_vals)) {
  R0 <- 12
  
  pars <- list(
    eps1 = 0.07,
    eps2 = min(1, b * 0.03 / 0.07),  # capped at 1
    gamma = 1/10,
    beta = (1/10) * R0,
    p1 = p1_vals[i],
    p2 = p1_vals[i],
    v1 = 1 / (365/2),
    v2 = a * 1 / (365 * 3.5),  # faster second dose
    m = 1 / (365 * 4.5),
    b = 1 / (365 * 4.5)
  )
  
  times <- seq(0, 365 * 10, by = 1)
  
  y0 <- c(
    S = (N - 1) / N,
    S1 = 0,
    S2 = 0,
    I = 1 / N,
    B1 = 0,
    B2 = 0,
    R = 0
  )
  
  out <- ode(y = y0, times = times, func = Two_dose, parms = pars,
             method = "lsoda", atol = 1e-8, rtol = 1e-8)
  
  final <- out[nrow(out), ]
  I_f  <- final["I"]
  B1_f <- final["B1"]
  B2_f <- final["B2"]
  
  I_tot_f <- I_f + B1_f + B2_f
  
  fv1_early <- ifelse(I_tot_f > 0, B1_f / I_tot_f, 0)
  fv2_early <- ifelse(I_tot_f > 0, B2_f / I_tot_f, 0)
  
  results_early <- rbind(results_early, data.frame(
    p1 = p1_vals[i],
    fv1_early = as.numeric(fv1_early),
    fv2_early = as.numeric(fv2_early),
    rv1_early = as.numeric(B1_f),
    rv2_early = as.numeric(B2_f)
  ))
  
  if (i %% 20 == 0) cat(sprintf("  Progress: %d%%\n", i))
}

cat("Early vaccination scenario complete.\n")

# =============================================================================
# Comparison Plots: On-time vs Early vaccination
# =============================================================================

# Combine data
comparison_df <- results %>%
  left_join(results_early, by = "p1")

# Plot: Breakthrough fractions comparison
pdf("../../output/figures/TwoDose_Comparison_Fractions.pdf", height = 6, width = 8)
print(ggplot(comparison_df) +
  geom_line(aes(x = p1, y = fv1, color = "On-time", linetype = "First dose"), linewidth = 1) +
  geom_line(aes(x = p1, y = fv2, color = "On-time", linetype = "Second dose"), linewidth = 1) +
  geom_line(aes(x = p1, y = fv1_early, color = "Too early", linetype = "First dose"), linewidth = 1) +
  geom_line(aes(x = p1, y = fv2_early, color = "Too early", linetype = "Second dose"), linewidth = 1) +
  scale_color_manual(values = c("On-time" = "black", "Too early" = "grey50")) +
  scale_linetype_manual(values = c("First dose" = "solid", "Second dose" = "dashed")) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = "Breakthrough Fraction",
       color = "Timing",
       linetype = "Dose") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right"))
dev.off()

# Plot: Total breakthrough fraction comparison
pdf("../../output/figures/TwoDose_Comparison_Total.pdf", height = 6, width = 8)
print(ggplot(comparison_df) +
  geom_line(aes(x = p1, y = fv1 + fv2, linetype = "On-time"), linewidth = 1.5) +
  geom_line(aes(x = p1, y = fv1_early + fv2_early, linetype = "Too early"), linewidth = 1.5) +
  scale_linetype_manual(values = c("On-time" = "solid", "Too early" = "dashed")) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = expression(paste("Total Breakthrough Fraction (", f[V1] + f[V2], ")")),
       linetype = "Timing") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 0.8)))
dev.off()

# Plot: Breakthrough cases comparison
pdf("../../output/figures/TwoDose_Comparison_Cases.pdf", height = 6, width = 8)
print(ggplot(comparison_df) +
  geom_line(aes(x = p1, y = rv1 * N, color = "First dose", linetype = "On-time"), linewidth = 1) +
  geom_line(aes(x = p1, y = rv2 * N, color = "Second dose", linetype = "On-time"), linewidth = 1) +
  geom_line(aes(x = p1, y = rv1_early * N, color = "First dose", linetype = "Too early"), linewidth = 1) +
  geom_line(aes(x = p1, y = rv2_early * N, color = "Second dose", linetype = "Too early"), linewidth = 1) +
  scale_color_manual(values = c("First dose" = "red", "Second dose" = "blue")) +
  scale_linetype_manual(values = c("On-time" = "solid", "Too early" = "dashed")) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = "Breakthrough Cases",
       color = "Dose",
       linetype = "Timing") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right"))
dev.off()

# Plot: Total breakthrough cases comparison
pdf("../../output/figures/TwoDose_Comparison_TotalCases.pdf", height = 6, width = 8)
print(ggplot(comparison_df) +
  geom_line(aes(x = p1, y = (rv1 + rv2) * N, linetype = "On-time"), linewidth = 1.5) +
  geom_line(aes(x = p1, y = (rv1_early + rv2_early) * N, linetype = "Too early"), linewidth = 1.5) +
  scale_linetype_manual(values = c("On-time" = "solid", "Too early" = "dashed")) +
  labs(x = expression(paste("First-dose vaccine coverage (", p[1], ")")),
       y = expression(paste("Total Breakthrough Cases (", I[V1]^"*" + I[V2]^"*", ")")),
       linetype = "Timing") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 0.8)))
dev.off()
