# Breakthrough Infections

Code and data to accompany "Anticipating and Interpreting Breakthrough Infections Amid Expanded Vaccine Refusal" by Mallory J Harris, Akash Arani, Tapan Goel, Kejia Zhang, Stephen J Beckett, Nathan C Lo, Jonathan Dushoff, and Joshua S Weitz.

---

## Project Structure

```
Breakthrough-Infections/
├── data/
│   ├── input/                      # Raw input data
│   │   ├── school-reports/         # School-level vaccination data (NOT INCLUDED)
│   │   ├── mmr_data_us_counties.csv
│   │   └── recent-outbreaks.csv
│   └── generated/                  # ODE simulation outputs
│       ├── ode_out.csv
│       ├── diffdisease-ode-output.csv
│       └── lowphi-*.csv
├── scripts/
│   ├── R/                          # R analysis scripts
│   │   ├── install_packages.R      # Package installation
│   │   ├── assortativity_estim.R   # Figure 3, phi estimation
│   │   ├── main_ode_plots.R        # Figures 1, 2, 4
│   │   ├── supplemental_figs.R     # Supplemental figures
│   │   └── 2-dose.R                # Two-dose model
│   └── matlab/                     # MATLAB ODE solvers
│       ├── main.m                  # Main ODE simulation -> ode_out.csv
│       ├── diffdisease.m           # Different disease params -> diffdisease-ode-output.csv
│       ├── lowphi.m                # Low assortativity sims -> lowphi-*.csv
│       ├── model_parameters.m
│       └── SIR_vaccinated_assortativity.m
├── output/
│   ├── figures/                    # Generated PDF figures
│   └── tables/                     # Generated CSV/TXT tables (NOT INCLUDED)
├── notebooks/
│   └── TwoDoseModel.ipynb
└── README.md
```

---

## Data Availability

### Public Data (Included)
- `data/input/mmr_data_us_counties.csv` - County-level MMR vaccination data
- `data/input/recent-outbreaks.csv` - Recent outbreak data
- `data/generated/` - ODE simulation outputs (can be regenerated)

### Restricted Data (Not Included)
- `data/input/school-reports/` - School-level vaccination data from state health departments
  - This data is not publicly redistributable due to privacy/licensing restrictions
  - Contact the authors for data access or obtain directly from state health departments
- `output/tables/` - Contains intermediate results derived from restricted school data

### Running Without School Data
If you don't have access to school-level data:
1. Skip Step 1 (`assortativity_estim.R`)
2. The MATLAB scripts will use default phi values automatically
3. Or manually create `output/tables/phi_estimates.csv` with your own estimates:
   ```csv
   "phi_mean","phi_lower","phi_upper"
   0.39,0.27,0.60
   ```

---

## Environment Setup

### R Requirements

**R Version:** 4.2.0 or higher recommended

**Required Packages:**
- tidyverse (data manipulation)
- readxl (Excel files)
- magrittr (pipe operators)
- cowplot (plot composition)
- usmap (US map visualization)
- latex2exp (LaTeX in plots)
- ggrepel (text labels)
- viridis (color palettes)
- scales (scale functions)
- ggridges (ridge plots)
- knitr, kableExtra (tables)
- deSolve (ODE solvers)
- ggh4x (extended facets)
- ggmagnify (plot insets)

**Installation:**
```r
# Run the installation script
source("scripts/R/install_packages.R")

# Or manually:
install.packages(c("tidyverse", "readxl", "magrittr", "cowplot", 
  "usmap", "latex2exp", "ggrepel", "viridis", "scales", 
  "ggridges", "knitr", "kableExtra", "deSolve", "ggh4x"))

# ggmagnify (special repo)
install.packages("ggmagnify", 
  repos = c("https://hughjonesd.r-universe.dev", 
            "https://cloud.r-project.org"))
```

**Troubleshooting:**
- ggh4x errors: `remotes::install_version("ggh4x", version = "0.2.8")`
- usmap errors: install sf first with `install.packages("sf")`

### MATLAB Requirements

**MATLAB Version:** R2020a or higher

**Required Toolboxes:**
- Parallel Computing Toolbox (for parpool/parfor)
- *Optional: modify scripts to use regular for loops*

**Check installation:**
```matlab
license('test', 'Distrib_Computing_Toolbox')  % Returns 1 if available
ver  % List all toolboxes
```

**Without Parallel Toolbox:** Edit main.m and diffdisease.m:
1. Comment out: poolobj = parpool(8);
2. Change: parfor -> for
3. Comment out: delete(poolobj);

---

## Running the Analysis

### Step 1: Estimate Assortativity (R)
```r
setwd("scripts/R")
source("assortativity_estim.R")
```

### Step 2: Run ODE Simulations (MATLAB)
```matlab
cd scripts/matlab
main        % -> data/generated/ode_out.csv
diffdisease % -> data/generated/diffdisease-ode-output.csv
lowphi      % -> data/generated/lowphi-ode-output.csv, lowphi-inset-ode-output.csv
```
(Timestamped backups are also saved in scripts/matlab/ for reproducibility)

### Step 3: Generate Figures (R)
```r
setwd("scripts/R")
source("main_ode_plots.R")
source("supplemental_figs.R")
source("2-dose.R")
```

---

## Data Flow

```
assortativity_estim.R -> phi_estimates.csv
                              |
        +---------------------+---------------------+
        v                     v                     v
    main.m              diffdisease.m          lowphi.m
        |                     |                     |
        v                     v                     v
  ode_out.csv     diffdisease-ode-output.csv   lowphi-*.csv
        |                     |                     |
        +---------------------+---------------------+
                              v
              main_ode_plots.R + supplemental_figs.R
                              |
                              v
                      output/figures/*.pdf
```

---

## Citation

Harris MJ, Arani A, Goel T, Zhang K, Beckett SJ, Lo NC, Dushoff J, Weitz JS. 
"Anticipating and Interpreting Breakthrough Infections Amid Expanded Vaccine Refusal"
