# on first run, install ggmagnify for easy inset
# install.packages(
#   "ggmagnify",
#   repos = c("https://hughjonesd.r-universe.dev", "https://cloud.r-project.org")
# )

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

require(deSolve)
require(magrittr)
require(tidyverse)
require(cowplot)
require(usmap)
require(latex2exp)
library(ggrepel)
library(scales)
library(viridis)
library(ggmagnify)

# Create output directories if they don't exist
if (!dir.exists("../../output/figures")) {
  dir.create("../../output/figures", recursive = TRUE)
}
if (!dir.exists("../../output/tables")) {
  dir.create("../../output/tables", recursive = TRUE)
}


# read in and transform the outbreak data
outbreak_df <- data.frame(State = state.name, Abb = state.abb) %>%
  right_join(read.csv("../../data/input/recent-outbreaks.csv")) %>%
  # if unknown and unvaccinated are not reported separately, assume half of combined category are unvaccinated
  mutate(Unknown = ifelse(is.na(Unknown), .5*Unknown.or.unvaccinated, Unknown),
         Unvaccinated = ifelse(is.na(Unvaccinated), .5*Unknown.or.unvaccinated, Unvaccinated)) %>%
  mutate(fv_mid = At.least.one.dose/(At.least.one.dose+Unvaccinated)) %>%
  # uncertainty: add upper and lower bound based on assumption of vaccination rate in known vs unknown populations
  mutate(fv_lower = (At.least.one.dose+fv_mid/3*Unknown)/(Unknown.or.unvaccinated+At.least.one.dose),
         fv_upper = (At.least.one.dose+fv_mid*3*Unknown)/(Unknown.or.unvaccinated+At.least.one.dose)) %>%
  rename(Coverage = State.Coverage) %>%
  mutate(Name = paste0(County, ", ", Abb))


# read in simulation output
df <- read.csv("../../data/generated/ode_out.csv") %>%
  data.frame() %>%
  rename(p = coverage,
         annual_CU = incidence_U,
         annual_CV = incidence_V,
         phi = assortativity) %>%
  # Round phi values to 2 decimal places for cleaner legend display
  mutate(phi = round(phi, 2)) %>%
  mutate(Underreported = "Neither")

df %<>% 
  rbind(df %>% 
          mutate(IU = (1/2)*IU) %>%
          mutate(Underreported="Unvaccinated Underreporting")) %>%
  rbind(rbind(df %>% 
                  mutate(IV = (1/2)*IV) %>%
                  mutate(Underreported = "Breakthrough Underreporting"))) %>%
  mutate(fV = IV/(IU+IV)) %>%
  filter(IU+IV >= 1)

# solve for HIT values (will use to draw caps on p vs f_V plot)
df_HIT <- df %>%
  filter(IU+IV >= 1) %>%
  group_by(phi) %>%
  slice_max(p) %>%
  rename(pHIT = p) %>%
  rename(fvHIT = fV) 

# Underreporting sensitivity
pdf("../../output/figures/supp-underrep.pdf", height=6, width=10)
print(ggplot(outbreak_df)+
        #predict relationship between p and fV by phi
        geom_line(data = df %>% filter(phi %in% c(0,.3, .6, .9, .98) & Underreported != "Neither"), aes(x=p, y=fV, color=as.factor(phi)))+
        facet_wrap(~Underreported)+
        # caps at herd immunity
        geom_point(data = df_HIT %>% filter(phi %in% c(0,.3, .6, .9, .98) & Underreported != "Neither"), aes(x=pHIT, y=fvHIT, color=as.factor(phi)))+
        # outbreak data as points
        geom_linerange(aes(x=Coverage, ymin=fv_lower, ymax=fv_upper))+
        geom_point(aes(x=Coverage, y=fv_mid, shape=Abb), size=2)+
        # state labels
        geom_text(data=outbreak_df %>% filter(Abb %in% c("TX", "NM")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x=.02, size=3)+
        geom_text(data=outbreak_df %>% filter(!Abb %in% c("TX", "NM", "SC", "MI")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x = -.06, nudge_y=.006, size=3)+
        geom_text(data=outbreak_df %>% filter(Abb %in% c("SC", "MI")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x=.02, nudge_y=.006, size=3)+
        # limits, labels, and theme
        scale_x_continuous(limits=range(c(0,1.1)), expand=c(0,0))+
        #scale_y_continuous(limits=range(c(-.02,.43)), expand=c(0,0))+
        scale_color_viridis_d()+
        theme_classic(base_size=20)+
        guides(color = guide_legend(override.aes=list(shape=NA),
                                    expression(paste("Assortativity (", italic(phi), ")")),
                                    title.position = "top"),
               shape="none")+
        theme(legend.position = "bottom", 
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 15, hjust=.5))+
        ylab(expression(paste("Breakthrough Fraction (", italic(f[V]), ")")))+
        xlab(expression(paste("Vaccine coverage (", italic(p),")")))+
        # manually set shapes for different states
        scale_shape_manual(values=c("MI"=15, "ND"=19,"UT"=17, "TX"=18, "NM"=8, 
                                    "AZ"=7, "SC"=18)))
dev.off()


# read in simulation output for different disease parameters
diffdis_df <- read.csv("../../data/generated/diffdisease-ode-output.csv") %>%
  # Round phi values to 2 decimal places for cleaner legend display
  mutate(assortativity = round(assortativity, 2))  %>%
  filter(assortativity %in% c(0, .6, .9, .98))

# supplemental figure - fV with different disease parameters
pdf("../../output/figures/supp-fv.pdf", height=10, width=10)

print(diffdis_df %>%
  data.frame() %>%
  rename(phi = assortativity) %>%
  # exclude small/NA values (indicative of elimination)
  filter(!is.nan(IV) & !is.nan(IU)) %>%
  filter(IU + IV >= 1) %>%
ggplot(aes(x=coverage, y=fV, color=as.factor(phi))) +
  geom_line()+
  scale_x_continuous(limits=range(c(0,1)), expand=c(0,0))+
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab(expression(paste("Breakthrough Fraction (", italic(f[V]), ")")))+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  scale_color_viridis_d()+
  # dynamic labelling of facets based on variable values
  facet_grid(rows = vars(vaccine_failure),
             cols = vars(R0),
             labeller = label_bquote(rows = epsilon ==.(vaccine_failure),
                                     cols = R[0] == .(R0)))+
  theme(panel.spacing.x = unit(3, "lines")))
dev.off()


# supplemental figure - fV with different disease parameters
pdf("../../output/figures/supp-Iv.pdf", height=10, width=10)
print(diffdis_df %>%
  data.frame() %>%
  rename(phi = assortativity) %>%
  ggplot(aes(x=coverage, y=incidence_V, color=as.factor(phi))) +
  geom_line()+
  # set limits to avoid space below 0 infections
  scale_x_continuous(limits=range(c(0,1)), expand=c(0,0))+
  coord_cartesian(ylim = c(0, NA), expand=FALSE)+ 
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab("Breakthrough Infections (annual per 100K)")+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  # dynamic labelling of facets based on variable values
  facet_grid(rows = vars(vaccine_failure),
             cols = vars(R0),
             labeller = label_bquote(rows = epsilon ==.(vaccine_failure),
                                     cols = R[0] == .(R0)))+
  scale_color_viridis_d()+
  theme(panel.spacing.x = unit(3, "lines")))
  
dev.off()

# supplemental figure - herd immunity threshold with mild assortativity

# read in and combine two dataframes here (inset gives smaller steps in coverage at high values)
# Use round(phi, 4) to preserve small phi values like 0.0001, 0.001
lowphi_df <- read.csv("../../data/generated/lowphi-ode-output.csv") %>% 
              data.frame() %>%
              mutate(phi = round(phi, 4))
lowphi_df <- read.csv("../../data/generated/lowphi-inset-ode-output.csv") %>% 
              data.frame() %>%
              mutate(phi = round(phi, 4)) %>%
              rbind(lowphi_df)


# plot larger figure with coverage from 0 to 1
lowphi_main <- lowphi_df %>%
  data.frame() %>%
  ggplot(aes(x=coverage, y=incidence_V, color=as.factor(phi))) +
  geom_line()+
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab("Breakthrough Infections (annual per 100K)")+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  coord_cartesian(ylim = c(0, max(lowphi_df$incidence_V)+1), xlim=c(0, 1.1), expand=FALSE)+ 
  scale_color_viridis_d()


pdf("../../output/figures/supp-lowphiHIT.pdf", height=8, width=8)

# use geo_magnify to add an inset
print(lowphi_main +
  geom_magnify(from=c(xmin=.95, xmax=1, ymin=0, ymax=10),
               to=c(xmin=.3, xmax=.85, ymin=1, ymax=8),
               axes="x"))
dev.off()

# ========================================
# Read estimated assortativity values from assortativity_estim.R output
# ========================================
phi_estimates <- read.csv("../../output/tables/phi_estimates.csv")
phi_mean <- phi_estimates$phi_mean
phi_lower <- phi_estimates$phi_lower
phi_upper <- phi_estimates$phi_upper

cat("\n========================================\n")
cat("Using estimated assortativity values:\n")
cat(sprintf("  phi_mean:  %.3f\n", phi_mean))
cat(sprintf("  phi_lower: %.3f\n", phi_lower))
cat(sprintf("  phi_upper: %.3f\n", phi_upper))
cat("========================================\n")

# ========================================
# Load ODE output (from Figure4.m with estimated phi values)
# ========================================
ode_df <- read.csv("../../data/generated/ode_out.csv") %>%
  data.frame() %>%
  rename(p = coverage,
         annual_CU = incidence_U,
         annual_CV = incidence_V,
         phi = assortativity) %>%
  # Round phi values to 2 decimal places for cleaner display
  mutate(phi = round(phi, 2))

# Find closest phi values in ODE output to estimated values
available_phi <- unique(ode_df$phi)
cat("\nAvailable phi values in ODE output:", paste(available_phi, collapse=", "), "\n")

# Function to find closest phi
find_closest_phi <- function(target, available) {
  available[which.min(abs(available - target))]
}

phi_mean_matched <- find_closest_phi(phi_mean, available_phi)
phi_lower_matched <- find_closest_phi(phi_lower, available_phi)
phi_upper_matched <- find_closest_phi(phi_upper, available_phi)

cat(sprintf("Matched phi values: mean=%.2f, lower=%.2f, upper=%.2f\n", 
            phi_mean_matched, phi_lower_matched, phi_upper_matched))

# Compute p_hat (coverage at peak breakthrough infections) for each phi
hatp_df <- ode_df %>%
  group_by(phi) %>%
  slice_max(annual_CV) %>%
  rename(p_hat = p) %>%
  select(phi, p_hat)

# Compute p_c (critical coverage for disease elimination) for each phi
# Use interpolation to find precise crossing point where total incidence = 1
compute_pc_interpolated <- function(df_phi) {
  df_phi <- df_phi %>% arrange(p)
  total_inc <- df_phi$annual_CU + df_phi$annual_CV
  
  # Check if incidence ever drops below 1 in simulated range
  if (min(total_inc) >= 1) {
    # Incidence never drops below 1 in simulated range
    # Interpolate between last simulated point and p=1 (where incidence=0)
    last_p <- max(df_phi$p)
    last_inc <- total_inc[which.max(df_phi$p)]
    if (last_inc > 1 && last_p < 1) {
      # Linear interpolation: find p where incidence crosses 1
      # At p=1, incidence=0; at last_p, incidence=last_inc
      # Solve: 1 = last_inc + (0 - last_inc) * (p_c - last_p) / (1 - last_p)
      p_c <- last_p + (1 - last_inc) / (0 - last_inc) * (1 - last_p)
      return(p_c)
    }
    return(NA)  # Truly never crosses 1
  }
  
  # Find first crossing point within simulated range
  idx <- which(total_inc < 1)[1]
  if (idx == 1) return(df_phi$p[1])
  
  # Linear interpolation between adjacent points where crossing occurs
  p1 <- df_phi$p[idx - 1]
  p2 <- df_phi$p[idx]
  inc1 <- total_inc[idx - 1]
  inc2 <- total_inc[idx]
  
  # Solve: 1 = inc1 + (inc2 - inc1) * (p_c - p1) / (p2 - p1)
  p_c <- p1 + (1 - inc1) / (inc2 - inc1) * (p2 - p1)
  return(p_c)
}

pc_df <- ode_df %>%
  group_by(phi) %>%
  group_modify(~ data.frame(p_c = compute_pc_interpolated(.x))) %>%
  ungroup()

# Combine p_hat and p_c
thresholds_df <- hatp_df %>%
  left_join(pc_df, by = "phi")

cat("\np_hat and p_c for each simulated phi:\n")
print(thresholds_df)

# ========================================
# Extract p_hat and p_c for estimated phi values (from ODE output)
# ========================================
# p_hat at mean estimated phi (or closest match)
p_hat_mean <- thresholds_df %>% filter(phi == phi_mean_matched) %>% pull(p_hat)
p_hat_lower <- thresholds_df %>% filter(phi == phi_lower_matched) %>% pull(p_hat)
p_hat_upper <- thresholds_df %>% filter(phi == phi_upper_matched) %>% pull(p_hat)

# Extract p_c values (now computed via interpolation)
p_c_baseline <- thresholds_df %>% filter(phi == 0) %>% pull(p_c)
p_c_at_mean <- thresholds_df %>% filter(phi == phi_mean_matched) %>% pull(p_c)
p_c_at_lower <- thresholds_df %>% filter(phi == phi_lower_matched) %>% pull(p_c)
p_c_at_upper <- thresholds_df %>% filter(phi == phi_upper_matched) %>% pull(p_c)

cat("\n========================================\n")
cat("Extracted threshold values from ODE output (p_c via interpolation):\n")
cat(sprintf("  p_hat at phi=%.2f (mean):  %.3f\n", phi_mean_matched, p_hat_mean))
cat(sprintf("  p_hat at phi=%.2f (lower): %.3f\n", phi_lower_matched, p_hat_lower))
cat(sprintf("  p_hat at phi=%.2f (upper): %.3f\n", phi_upper_matched, p_hat_upper))
cat(sprintf("  p_c at phi=0 (baseline):   %.4f\n", p_c_baseline))
cat(sprintf("  p_c at phi=%.2f (mean):    %.4f\n", phi_mean_matched, p_c_at_mean))
cat(sprintf("  p_c at phi=%.2f (lower):   %.4f\n", phi_lower_matched, p_c_at_lower))
cat(sprintf("  p_c at phi=%.2f (upper):   %.4f\n", phi_upper_matched, p_c_at_upper))
cat("========================================\n")

# Save threshold values (now with interpolated p_c for all phi values)
threshold_estimates <- data.frame(
  description = c("phi=0 baseline", "phi_mean", "phi_lower", "phi_upper"),
  phi = c(0, phi_mean_matched, phi_lower_matched, phi_upper_matched),
  p_hat = c(thresholds_df %>% filter(phi == 0) %>% pull(p_hat), p_hat_mean, p_hat_lower, p_hat_upper),
  p_c = c(p_c_baseline, p_c_at_mean, p_c_at_lower, p_c_at_upper)
)
write.csv(threshold_estimates, file = "../../output/tables/threshold_estimates.csv", row.names = FALSE)
cat("\nSaved threshold estimates to: ../../output/tables/threshold_estimates.csv\n")

# Also save full thresholds_df with all phi values
write.csv(thresholds_df, file = "../../output/tables/thresholds_all_phi.csv", row.names = FALSE)
cat("Saved full thresholds table to: ../../output/tables/thresholds_all_phi.csv\n")

# ========================================
# Figure 5: County map with data-driven thresholds
# Use p_hat and p_c at mean estimated phi (both from interpolated ODE output)
# ========================================
p_hat_map <- round(p_hat_mean, 3)
p_c_map <- round(p_c_at_mean, 4)  # Use 4 decimal places for p_c (very close to 1)

cat(sprintf("\nFor county map: p_hat = %.3f, p_c = %.4f (both at phi=%.2f)\n", 
            p_hat_map, p_c_map, phi_mean_matched))

# read in csv of county-level MMR vaccine coverage
pdf("../../output/figures/county_map.pdf", height=10, width=8)

map1 <- read.csv("../../data/input/mmr_data_us_counties.csv") %>%
  rename(fips = FIPS) %>%
  # bin values based on data-driven p_c and p_hat from estimated assortativity
  mutate(category = case_when(
    SY2022_23 <= p_hat_map ~ "1",
    SY2022_23 > p_hat_map & SY2022_23 <= p_c_map ~ "2",
    SY2022_23 > p_c_map ~ "3"
  )) %>%
  plot_usmap(data=., values="category", color="grey", size=.1)+
  theme_void(base_size=20)+
  theme(legend.position="bottom")+
  scale_fill_viridis_d(name="",
                      option="A",
                      # allow latex special characters in map
                      labels = c("3" = TeX("$p > p_{c}"), 
                                 "2" = TeX("$\\hat{p} < p \\leq p_{c}$"), 
                                 "1" = TeX("$p \\leq \\hat{p}$")),
                      guide=guide_legend(reverse=FALSE),
                      na.value="white")

county_data <- read.csv("../../data/input/mmr_data_us_counties.csv")

# bottom panel, histogram with data-driven threshold lines
map2 <- ggplot()+
  geom_histogram(data=county_data, aes(x=SY2022_23), fill="grey70", closed="left")+
  ylab("Count (Counties)")+
  xlab("Vaccine Coverage (p)")+
  # solid line for p_hat at mean estimated phi
  geom_segment(aes(x=p_hat_map, y=0, yend=470))+
  # solid line for p_c (baseline at phi=0)
  geom_segment(aes(x=p_c_map, y=0, yend=470))+
  # dotted lines show p_hat range for phi_lower and phi_upper
  geom_segment(aes(x=p_hat_lower, y=0, yend=400), linetype="dotted", alpha=0.7)+
  geom_segment(aes(x=p_hat_upper, y=0, yend=400), linetype="dotted", alpha=0.7)+
  geom_segment(aes(x=p_c_at_lower, y=0, yend=400), linetype="dotted", alpha=0.7)+
  geom_segment(aes(x=p_c_at_upper, y=0, yend=400), linetype="dotted", alpha=0.7)+
  scale_x_continuous(expand = expansion(mult = c(0, .1)),
                     breaks=seq(from=0, to=1, by=.1))+
  scale_y_continuous(limits=c(0, 650),expand = expansion(mult = c(0, 0)))+
  theme_classic(base_size=20)+
  annotate("text", x=p_c_map, y=550, size=8, label="p[c]", parse=TRUE)+
  annotate("text", x=p_hat_map, y=550, size=8, label="hat(p)", parse=TRUE)+
  theme(axis.text.x = element_text(hjust = 0),
        axis.ticks.x  = element_blank())

print(plot_grid(map1, map2, rel_heights=c(6, 2), nrow =2))
dev.off()


# get values for number of counties above/below threshold value for coverage
n_above_pc <- read.csv("../../data/input/mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 >= p_c_map) %>%
  nrow()

n_between <- read.csv("../../data/input/mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 < p_c_map & SY2022_23 > p_hat_map) %>%
  nrow()

n_below_phat <- read.csv("../../data/input/mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 <= p_hat_map) %>%
  nrow()

n_below_phat_lower <- read.csv("../../data/input/mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 <= p_hat_lower) %>%
  nrow()

n_below_phat_upper <- read.csv("../../data/input/mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 <= p_hat_upper) %>%
  nrow()

cat("\n========================================\n")
cat("County counts by coverage category:\n")
cat(sprintf("  p > p_c (%.4f):           %d counties\n", p_c_map, n_above_pc))
cat(sprintf("  p_hat < p <= p_c:         %d counties\n", n_between))
cat(sprintf("  p <= p_hat (%.3f):        %d counties\n", p_hat_map, n_below_phat))
cat(sprintf("  p <= p_hat_lower (%.3f):        %d counties\n", p_hat_lower, n_below_phat_lower))
cat(sprintf("  p <= p_hat_upper (%.3f):        %d counties\n", p_hat_upper, n_below_phat_upper))
cat("========================================\n")
