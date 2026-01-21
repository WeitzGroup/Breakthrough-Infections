# ============================================================================
# Set working directory to scripts/R folder
# This ensures relative paths work correctly regardless of how script is called
# ============================================================================
# Find the project root by looking for README.md
find_project_root <- function() {
  current <- getwd()
  while (current != dirname(current)) {  # Not at filesystem root
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
} else {
  cat("Warning: Could not find project root. Please run from scripts/R directory.\n")
  cat("Current directory:", getwd(), "\n")
}

library(tidyverse)
library(readxl)
library(ggh4x)
library(knitr)
library(kableExtra)
library(magrittr)
library(ggridges)

# Create output directories if they don't exist
if (!dir.exists("../../output/tables")) {
  dir.create("../../output/tables", recursive = TRUE)
}
if (!dir.exists("../../output/figures")) {
  dir.create("../../output/figures", recursive = TRUE)
}


# Initialize Functions ------------------------------------------------------------
# this file helps streamline reading in .xlsx files with school-level vaccination data
read_xl_flex <- function(file, sheet, skip, year){
  df <- read_xlsx(file, sheet=sheet, skip=skip) %>%
        mutate(School_year=year)
  return(df)
}

# Function to calculate population-weighted assortativity from school-level data
# Takes a dataframe with K_enrollment, school_p, School_year columns
# Returns the estimated assortativity (phi) for each School_year
calc_assortativity <- function(df) {
  df %>%
    group_by(School_year) %>%
    mutate(total_K = sum(K_enrollment)) %>% 
    # proportion of students in each school (will use for population-weighted mean)
    mutate(school_weight = K_enrollment / total_K) %>%
    # calculated population-weighted mean coverage
    mutate(weighted_p = school_p * school_weight) %>%
    group_by(School_year) %>%
    mutate(state_p = sum(weighted_p)) %>%
    # estimate school-level assortativity
    mutate(school_phi = ifelse(school_p > state_p, 
                                (school_p - state_p) / (1 - state_p), 
                                (state_p - school_p) / state_p)) %>%
    # population-weighted average estimated assortativity
    mutate(weighted_phi = school_phi * school_weight) %>%
    group_by(School_year) %>%
    summarize(phi = sum(weighted_phi), n_school = n(), .groups = "drop")
}

# function that takes in a dataframe with K_enrollment, school_p, School_year and a state name
# returns year, state, estim_phi (main estimate), across 1000 rows with different bootstrapped estimates (bootstrap_phi)
process_df <- function(df, state){
   # Calculate main estimate using the assortativity function
   calc_df <- calc_assortativity(df) %>%
     rename(estim_phi = phi) %>%
     mutate(year = as.integer(stringr::str_extract(as.character(School_year), "\\d{4}")),
            state = state)
   
   # bootstrapped estimates
   resampled_estim <- data.frame()
   
   # for each of 1000 runs, resample (with replacement) from schools & calculate phi
   set.seed(0115)
   for(i in 1:1000){
     resampled_estim <- df %>%
       group_by(School_year) %>%
       # resample with replacement to get same number of schools
       slice_sample(prop = 1, replace = TRUE) %>%
       ungroup() %>%
       calc_assortativity() %>%
       rename(bootstrap_phi = phi) %>%
       select(School_year, bootstrap_phi) %>%
       rbind(resampled_estim)
   }

   df <- calc_df %>% 
     left_join(resampled_estim, by = "School_year") %>%
     select(-School_year)
     
   return(df)
}

# initialize a dataframe that we will keep appending to with state-level dataframes below
assort_df <- data.frame()

# snippets below read in school-level vaccination data, ensure p_school is a numeric decimal
# select the desired column & standardize names, and filter out NA values
# then apply the process_df function to estimate phi (including bootstraps) and rbind to assort_df

# California ------------------------------------------------------------
assort_df <- read.csv("../../data/input/school-reports/CA-schools.csv") %>% 
              filter(CATEGORY == "MMR2") %>%
              mutate(School_year = substr(SCHOOL_YEAR, 1, 4)) %>%
              mutate(school_p = PERCENT/100) %>%
              rename(K_enrollment = ENROLLMENT) %>%
              filter(!is.na(school_p)) %>%
              process_df("CA") %>%
              rbind(assort_df)

# Colorado ------------------------------------------------------------
assort_df <- read.csv("../../data/input/school-reports/CO-schools.csv") %>%
              filter(Survey_Type=="Kindergarten" & Vaccine == "MMR" & Metric == "Fully Immunized" ) %>%
              mutate(School_year = substr(Year_, 1, 4),
                     school_p = Value_Percent/100) %>%
              rename(K_enrollment = Enrollment) %>%
              process_df("CO") %>%
              rbind(assort_df)

# Iowa ------------------------------------------------------------
assort_df <- lapply(2:10, function(n) read_xls("../../data/input/school-reports/IA-schools.xls", sheet=n) %>%
         rename(K_enrollment = `Total Enrolled`,
                school_p = `Immunization Certificates Percent`) %>%
         mutate(School_year = excel_sheets("../../data/input/school-reports/IA-schools.xls")[n] %>% substr(., 1, 4)) %>% # get schoolyears from sheet names
         select(K_enrollment, school_p, School_year)) %>%
  do.call(rbind,.) %>%
  mutate(school_p = school_p/100) %>%
  process_df("IA") %>%
  rbind(assort_df)
  
  
# Kentucky ------------------------------------------------------------
assort_df <- read_xlsx("../../data/input/school-reports/KY-schools.xlsx", sheet=2) %>%
  rename(school_p = `2+ MMR (%)`,
         K_enrollment = `Total Kindergaren Students Enrolled`) %>%
  mutate(School_year = substr(`School Year`, 1, 4),
         school_p = school_p/100) %>% 
  process_df("KY") %>%
  rbind(assort_df)

  #process_df("KY")
# Maryland ------------------------------------------------------------
assort_df <- mapply(read_xl_flex, 
       dir("../../data/input/school-reports/MD/", full.names=TRUE), 
       c(rep(2, 8), rep(1, 1)), # number of heading rows varies between files
       c(rep(5, 8), rep(4, 1)),
       2015:2023) %>%
  lapply(function(df) df %>% 
           rename_with(toupper) %>%
           select(`TOTAL K STUDENTS`, `% MMR`, `SCHOOL_YEAR`, `% MEDICAL EXEMPTION`, `% RELIGIOUS EXEMPTION`)) %>%
  do.call(rbind, .) %>%
  rename(School_year = SCHOOL_YEAR,
         K_enrollment = `TOTAL K STUDENTS`,
         p_MMR = 2,
         p_med = 4,
         p_rel = 5) %>%
  #calculate p assuming exemptions are not vaccinated
  mutate(school_p = (1-as.numeric(p_rel)-as.numeric(p_med))*as.numeric(p_MMR), 
         K_enrollment = as.numeric(K_enrollment)) %>%
  filter(!is.nan(school_p) & !is.nan(K_enrollment) & K_enrollment!=0 & school_p >= 0) %>%
  process_df("MD") %>%
  rbind(assort_df)

# Massachusetts ------------------------------------------------------------
assort_df <- lapply(2:6, function(page) read_xlsx("../../data/input/school-reports/MA-schools.xlsx", sheet=page) %>%
         rename(MMR_count = `Count of MMR Req Met`,
                K_enrollment = `Child Count`) %>%
         mutate(School_year = substr(excel_sheets("../../data/input/school-reports/MA-schools.xlsx")[page],1,4)) %>%
         select(School_year, MMR_count, K_enrollment)) %>%
      do.call(rbind, .) %>%
      # convert count to percent
      mutate(school_p = MMR_count/K_enrollment) %>%
      filter(!is.na(school_p)) %>%
      process_df("MA") %>%
      rbind(assort_df)

# Michigan ------------------------------------------------------------
assort_df <- lapply(dir("../../data/input/school-reports/MI/", full.names=TRUE), function(file) read_xlsx(file, skip=7) %>%
                      mutate(School_year = str_extract(file, "\\b\\d{4}(?=\\.xlsx)")) %>% # get school year from file names
                      rename(school_p = `%COMP`,
                             K_enrollment = N) %>%
                      mutate(school_p = school_p/100) %>%
                      select(NAME, DISTRICT, School_year, K_enrollment, school_p)) %>%
  do.call(rbind, .) %>%
  process_df("MI") %>%
  rbind(assort_df)

# Minnesota ------------------------------------------------------------
assort_df <- mapply(read_xl_flex, 
       dir("../../data/input/school-reports/MN/", full.names=TRUE), 
       c(2,2,4,3,4,4), # inconsistent headers
       c(rep(1,3), 7, rep(1, 2)),
       c(2023, 2024, 2019:2022)) %>%
  lapply(function(df) df %>% select(`Kindergarten Enrollment`, `School District`, `MMR % Vaccinated`, School_year)) %>%
  do.call(rbind, .) %>%
  rename(school_p = 3,
         K_enrollment = 1) %>%
  filter(`School District` != "Statewide") %>%
  mutate(school_p = as.numeric(school_p), 
         K_enrollment = as.numeric(K_enrollment)) %>% 
  filter(!is.na(school_p) & !is.na(K_enrollment)) %>%
  process_df("MN") %>%
  rbind(assort_df)

# Missouri ------------------------------------------------------------
assort_df <- read_xlsx("../../data/input/school-reports/MO-schools.xlsx") %>%
              filter(K_Imm!="Suppressed") %>% 
              mutate(K_enrollment = as.numeric(`K Enroll`),
                     school_p=as.numeric(K_Rate)) %>%
              filter(!is.na(school_p)) %>%
              mutate(School_year=substr(`School Year`, 1, 4)) %>%
              process_df("MO") %>%
              rbind(assort_df)

# New York State ------------------------------------------------------------
assort_df <- lapply(dir("../../data/input/school-reports/NY/", full.names=TRUE),
                    function(file) read_xlsx(file) %>%
                      rename(K_enrollment = 1, num_MMR = 2) %>%
                      select(K_enrollment, num_MMR) %>%
                      mutate(
                        School_year = stringr::str_extract(basename(file), "\\d{4}"),
                        K_enrollment = as.numeric(K_enrollment),
                        num_MMR = as.numeric(num_MMR)
                      )
) %>%
  do.call(rbind, .) %>%
  mutate(school_p = num_MMR / K_enrollment) %>%
  process_df("NY*") %>%
  rbind(assort_df)
  
# North Carolina ------------------------------------------------------------
assort_df <-lapply(dir("../../data/input/school-reports/NC/", full.names=TRUE), 
              function(file) read_xlsx(file) %>%
                            select(`Total Enrollment`, `Up to Date (%)`) %>%
                            rename(K_enrollment = 1,
                                   school_p = 2) %>%
                            # get school year from file name
                            mutate(School_year = str_extract(file, "(?<=[-_])\\d{4}(?=\\.)"))) %>%
              do.call(rbind, .) %>%
              filter(!is.na(school_p)) %>%
              process_df("NC") %>%
              rbind(assort_df)

# North Dakota ------------------------------------------------------------
assort_df <- read_xlsx("../../data/input/school-reports/ND-schools.xlsx") %>%
  rename(School_year = 1,
         school_p = `% UTD MMR`,
         K_enrollment = Enrolled) %>%
  filter(!is.na(school_p)) %>%
  # convert percent to decimal
  mutate(school_p = school_p/100) %>%
  mutate(School_year = substr(School_year, 1, 4)) %>%
  process_df("ND") %>%
  rbind(assort_df)

# Oregon ------------------------------------------------------------
assort_df <- lapply(dir("../../data/input/school-reports/OR/", full.names=TRUE), function(file) read_xlsx(file) %>%
         # get school year from file name
         mutate(School_year = str_extract(file, "\\b(19|20)\\d{2}\\b")) %>% 
         rename(school_p = 4,
                K_enrollment = 2) %>%
         select(School_year, K_enrollment, school_p)) %>%
  do.call(rbind, .) %>%
  # convert percent to decimal
  mutate(school_p = ifelse(str_ends(school_p, "%"), as.numeric(gsub("%", "", school_p))/100, as.numeric(school_p))) %>%
  mutate(school_p = as.numeric(school_p), 
         K_enrollment = as.numeric(K_enrollment)) %>% 
  filter(!is.na(school_p) & !is.na(K_enrollment)) %>%
  process_df("OR")  %>%
  rbind(assort_df)


# South Carolina ------------------------------------------------------------
sc_path <- "../../data/input/school-reports/SC-schools.xlsx"
assort_df <- lapply(1:6, function(n) read_xlsx(sc_path, sheet=n, skip=2) %>%
         mutate(School_year = substr(excel_sheets(sc_path)[n], 1, 4)) %>%
         rename(school_p = 5,
                K_enrollment = 3) %>%
         select(school_p, K_enrollment, School_year, Type)) %>%
  do.call(rbind,.) %>%
  filter(!is.na(Type)) %>%
  process_df("SC") %>%
  rbind(assort_df)

# Utah ------------------------------------------------------------
assort_df <- read_xlsx("../../data/input/school-reports/UT-schools.xlsx") %>%
  rename(school_p = 6) %>%
  mutate(School_year = substr(School_year, 1, 4)) %>%
  process_df("UT") %>%
  rbind(assort_df)


# Washington ------------------------------------------------------------
assort_df <- lapply(dir("../../data/input/school-reports/WA/", full.names=TRUE), 
       function(file) read_xlsx(file) %>% 
         select(Grade, `School Year`, `Disease or Vaccine`, Enrollment, Percent, "Immunization Status") %>%
         filter(Grade=="Kindergarten" & `Immunization Status`=="Complete" & `Disease or Vaccine` =="Measles")) %>%
  do.call(rbind, .) %>%
  rename(school_p = Percent,
         K_enrollment = Enrollment) %>%
  mutate(School_year = substr(`School Year`,1, 4)) %>%
  process_df("WA") %>%
  rbind(assort_df)
         
# write assort_df to file (can take some time to run)
write.csv(assort_df, file="../../output/tables/estim_assort_df.csv", row.names=FALSE)
 

# Output plots and table ------------------------------------------------------------        
assort_df <- read.csv("../../output/tables/estim_assort_df.csv") %>%
  group_by(year, state) %>%
  summarize(
    estim_phi = unique(estim_phi)[1],  # or mean(estim_phi)
    lower_phi = quantile(bootstrap_phi, 0.025, na.rm=TRUE),
    upper_phi = quantile(bootstrap_phi, 0.975, na.rm=TRUE),
    .groups = "drop"
  )

# set shapes for states (as in Fig 1 - 2)
assort_df_plot <- data.frame("state"=c("MI","ND","UT","SC"),
                             "state_shape"=c(15,19,17,18)) %>%
  right_join(assort_df) %>%
  mutate(state_shape = ifelse(is.na(state_shape), 1, state_shape),
         # set colors (black for states included in Figures 1 and 2 with 2025 outbreak data)
         state_color = ifelse(state %in% c("MI","ND","UT","SC"),"black",  "grey50"))

# evaluate trends since 2019
assort_df_plot %<>%
  filter(year >= 2019) %>%
  group_by(state) %>%
  summarize(slope_min = confint(lm(estim_phi ~ year), level=.95, parm="year")[1],
            slope_max = confint(lm(estim_phi ~ year), level=.95, parm="year")[2],
            slope = coef(lm(estim_phi ~ year))[2],
            pval = summary(lm(estim_phi ~ year))$coefficients["year", "Pr(>|t|)"],
            mean_phi = mean(estim_phi, na.rm=TRUE)) %>%
  right_join(assort_df_plot) 

# trends table
assort_table <- assort_df_plot %>% 
  select(state, mean_phi, slope, slope_min, slope_max, pval) %>%
  distinct() %>%
  # format values as numbers with two significant figures
  mutate(`Slope (95% CI)` = paste0(sprintf("%.2f", slope), " (", sprintf("%.2f", slope_min), ", ", sprintf("%.2f", slope_max), ")")) %>%
  mutate(`p-value` = sprintf("%.2f", pval),
         `Mean phi` = sprintf("%.2f", mean_phi)) %>%
  select(1, `Mean phi`, `Slope (95% CI)`, `p-value`) %>%
  kable(format="latex", row.names=FALSE)

writeLines(assort_table, con = "../../output/tables/slope_table.txt")


pdf("../../output/figures/supp-phi-trend.pdf", height=12, width=12)

print(ggplot(assort_df_plot)+
  # show 95% CI and trendline (separated out for different opacities)
  geom_smooth(data=assort_df_plot %>% filter(year>=2019), aes(x=year, y=estim_phi, color=NA, fill=state_color), alpha=.2, method="lm") +
  stat_smooth(data=assort_df_plot %>% filter(year>=2019), aes(x=year, y=estim_phi, color=state_color), geom="line", se=FALSE, alpha=.5, method="lm") +
  # estimated phi & band for each year
  geom_point(aes(x=year, y=estim_phi, color=state_color, shape=state_shape))+
  geom_segment(aes(x=year, y=lower_phi, yend=upper_phi, color=state_color))+
  # theme, labels, legend
  theme_bw(base_size=20)+
  ylab(expression(paste("Estimated Assortativity (", italic(tilde(phi)), ")")))+
  theme(panel.grid.major = element_line(color = "grey50", linetype = "dotted", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", linewidth = 0.2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits=c(0, NA))+
  scale_color_identity()+
  scale_fill_identity()+
  scale_shape_identity()+
  xlab("Year")+
  facet_wrap(~state, nrow=4))
dev.off()

# Plot Figure 3
# filter to 2019 and beyond
assort_hist <- read.csv("../../output/tables/estim_assort_df.csv") %>%
  filter(year >= 2019) %>%
  # color darker if in outbreak analysis
  mutate(state_color = ifelse(state %in% c("MI","ND","UT","SC"),"black",  "grey50")) %>%
  group_by(state) 

pdf("../../output/figures/fig3.pdf", height=8, width=8)
print(assort_hist %>%
  ggplot()+
  geom_density_ridges(aes(x=bootstrap_phi,
                                    y=reorder(state, estim_phi), 
                                    fill=state_color),
                      color=NA, alpha=.5)+
  scale_fill_identity()+
  theme_classic(base_size=20)+
  scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, by=.1))+
  ylab("State")+
  xlab(expression(paste("Estimated Assortativity (", italic(tilde(phi)), ")"))))
dev.off()

# get summary statistics for main text and store as constants
phi_summary <- read.csv("../../output/tables/estim_assort_df.csv") %>%
  filter(year >= 2019) %>% 
  summarize(
    phi_mean = mean(estim_phi, na.rm = TRUE),
    phi_median = median(estim_phi, na.rm = TRUE),
    phi_lower = min(estim_phi, na.rm = TRUE),
    phi_upper = max(estim_phi, na.rm = TRUE)
  )

# Print summary
cat("\n========================================\n")
cat("Estimated Assortativity Summary (2019+):\n")
cat("========================================\n")
cat(sprintf("  Mean:   %.3f\n", phi_summary$phi_mean))
cat(sprintf("  Median: %.3f\n", phi_summary$phi_median))
cat(sprintf("  Range:  [%.3f, %.3f]\n", phi_summary$phi_lower, phi_summary$phi_upper))
cat("========================================\n")

# Save assortativity estimates to file for use by other scripts
phi_estimates <- data.frame(
  phi_mean = phi_summary$phi_mean,
  phi_lower = phi_summary$phi_lower,
  phi_upper = phi_summary$phi_upper
)
write.csv(phi_estimates, file = "../../output/tables/phi_estimates.csv", row.names = FALSE)
cat("\nSaved phi estimates to: ../../output/tables/phi_estimates.csv\n")


