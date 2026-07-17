Sepsis-Associated Hypoxemia Trajectory Study
============================================

R analysis scripts for a multicohort retrospective study of sepsis-associated
hypoxemia. The workflow fits group-based multi-trajectory models in MIMIC-IV and
evaluates in-hospital mortality with competing-risk models in MIMIC-IV, eICU-CRD,
and MIMIC-III.

Scripts
-------

| File | Purpose |
|---|---|
| `R/01_trajectory/01_gbmtm_fit_six_scenarios.R` | Fits GBMTM models across six variable scenarios. |
| `R/01_trajectory/02_trajectory_figures.R` | Generates trajectory figures. |
| `R/01_trajectory/03_model_selection_tables.R` | Builds model-selection tables. |
| `R/01_trajectory/04_bootstrap_stability_K2_7.R` | Runs bootstrap stability analysis. |
| `R/02_competing_risk/05_competing_risk_main_cohort.R` | Runs competing-risk analysis in the main cohort. |
| `R/02_competing_risk/06_external_validation_eicu.R` | Applies external validation in eICU-CRD. |
| `R/02_competing_risk/07_external_validation_mimic3.R` | Applies external validation in MIMIC-III. |

Requirements
------------

- R >= 4.4
- R packages: `dplyr`, `tibble`, `purrr`, `tidyr`, `gbmt`, `ggplot2`,
  `patchwork`, `ggridges`, `scales`, `cmprsk`, `survival`, `prodlim`,
  `riskRegression`, `forestploter`, `mice`, `gtsummary`, `gt`

Run Order
---------

1. `R/01_trajectory/01_gbmtm_fit_six_scenarios.R`
2. `R/01_trajectory/02_trajectory_figures.R`
3. `R/01_trajectory/03_model_selection_tables.R`
4. `R/01_trajectory/04_bootstrap_stability_K2_7.R`
5. `R/02_competing_risk/05_competing_risk_main_cohort.R`
6. `R/02_competing_risk/06_external_validation_eicu.R`
7. `R/02_competing_risk/07_external_validation_mimic3.R`

Data
----

Patient-level data are not included. Prepare derived per-stay and day-level input
tables from MIMIC-IV, eICU-CRD, and MIMIC-III before running the scripts.

License
-------

MIT License.
