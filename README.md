# Breakthrough Infections

Code and data to accompany "Interpreting Breakthrough Infections Given Assortative Mixing of Partially Vaccinated Populations" by Mallory J Harris, Akash Arani, Tapan Goel, Kejia Zhang, Stephen J Beckett, Nathan C Lo, Jonathan Dushoff, and Joshua S Weitz.

---

## Project Structure

```
Breakthrough-Infections/
├── data/
│   ├── input/                      # Raw input data
│   │   ├── school-reports/         # School-level vaccination data (included for select states)
│   │   ├── mmr_data_us_counties.csv
│   │   └── recent-outbreaks.csv
│   └── generated/                  # ODE simulation outputs
│       ├── ode_out.csv
│       ├── diffdisease-ode-output.csv
│       ├── lowphi-*.csv
│       └── manyphi-ode-output.csv
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
│       ├── manyphi.m               # Sims across multiple phi values for heatmap -> manyphi-ode-output.csv
│       ├── model_parameters.m
│       └── SIR_vaccinated_assortativity.m
├── output/
│   ├── figures/                    # Generated PDF figures
│   └── tables/                     # Generated CSV/TXT tables (NOT INCLUDED)
├── notebooks/
│   ├── underrep-heatmap.ipynb                  
│   └── TwoDoseModel.ipynb
└── README.md
```

---

## Data Availability

### Public Data (Included)
- `data/input/mmr_data_us_counties.csv` - County-level MMR vaccination data
- `data/input/recent-outbreaks.csv` - Recent outbreak data
- `data/generated/` - ODE simulation outputs (can be regenerated)
- `data/input/school-reports/` - School-level vaccination data from state health departments where publicly available: CA,CO, MD, MN, NC, and WA

### Restricted Data (Not Included)
- `data/input/school-reports/` - School-level vaccination data from state health departments in IA, KY, MA, MI, MO, ND, NY, OR, SC, and UT
  - This data is not publicly redistributable due to privacy/licensing restrictions
  - Contact the authors for data access information or obtain directly from state health departments
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

**R Version:** 4.2.3

**Required Packages:**
| Package | Version | Description |
|---------|---------|-------------|
| tidyverse | 2.0.0 | Data manipulation and visualization |
| readxl | 1.4.5 | Excel file reading |
| magrittr | 2.0.4 | Pipe operators |
| cowplot | 1.1.1 | Plot composition |
| usmap | 0.6.4 | US map visualization |
| latex2exp | 0.9.8 | LaTeX expressions in plots |
| ggrepel | 0.9.6 | Non-overlapping text labels |
| viridis | 0.6.5 | Color palettes |
| scales | 1.4.0 | Scale functions |
| ggridges | 0.5.7 | Ridge plots |
| knitr | 1.51 | Report generation |
| kableExtra | 1.4.0 | Table formatting |
| deSolve | 1.40 | ODE solvers |
| ggh4x | 0.2.8 | Extended ggplot2 facets |
| ggmagnify | 0.4.2 | Plot insets/magnification |
| stringr | 1.5.1 | String manipulation |

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
- For exact version reproducibility: use `renv` or install specific versions with `remotes::install_version()`

### MATLAB Requirements

**MATLAB Version:** R2024a

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
manyphi     % -> data/generated/manyphi-ode-output.csv
```
(Timestamped backups are also saved in scripts/matlab/ for reproducibility)

### Step 3: Generate Figures (R, Python)
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
        +---------------------+---------------------+---------------------+
        v                     v                     v                     v
    main.m              diffdisease.m          lowphi.m              manyphi.m
        |                     |                     |                     |
        v                     v                     v                     v
  ode_out.csv     diffdisease-ode-output.csv   lowphi-*.csv     manyphi-ode-output.csv
        |                     |                     |                     |
        +---------------------+---------------------+---------------------+ 
                              v
                              |
                              v
                      output/figures/*.pdf

```

---

## Citation

Harris MJ, Arani A, Goel T, Zhang K, Beckett SJ, Lo NC, Dushoff J, Weitz JS. 
"Interpreting Breakthrough Infections Given Assortative Mixing of Partially Vaccinated Populations"
