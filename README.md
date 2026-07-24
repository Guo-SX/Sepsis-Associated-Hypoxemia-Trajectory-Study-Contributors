Sepsis-Associated Hypoxemia Trajectory Study
============================================

R analysis scripts for a multicohort retrospective study of sepsis-associated
hypoxemia. The workflow fits group-based multi-trajectory models in MIMIC-IV.
Competing-risk scripts are retained for reproducibility and endpoint auditing,
but run only when coherent hospital death and live-discharge times are supplied.

Scripts
-------

| File | Purpose |
|---|---|
| `R/01_trajectory/01_gbmtm_fit_six_scenarios.R` | Fits GBMTM models across six variable scenarios. |
| `R/01_trajectory/02_trajectory_figures.R` | Generates trajectory figures. |
| `R/01_trajectory/03_model_selection_tables.R` | Builds model-selection tables. |
| `R/01_trajectory/04_bootstrap_stability_K2_7.R` | Runs bootstrap stability analysis. |
| `R/02_competing_risk/05_competing_risk_main_cohort.R` | Runs the main-cohort competing-risk analysis after verifying hospital terminal-event times. |
| `R/02_competing_risk/06_external_validation_eicu.R` | Runs the eICU competing-risk workflow only when coherent hospital event times are available. |
| `R/02_competing_risk/07_external_validation_mimic3.R` | Runs the MIMIC-III competing-risk workflow only when coherent hospital event times are available. |

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

Competing-risk endpoint contract
--------------------------------

The Fine-Gray scripts require `hospital_event_time_hours`, measured from one
explicit analysis origin:

- for hospital deaths, the value must correspond to the verified death time;
- for survivors, the value must correspond to the live hospital-discharge time;
- `hosp_mortality` supplies the matching event type.

The scripts intentionally stop if `icu_los_hours` is selected as the time scale
for `hosp_mortality`. ICU exit is not interchangeable with hospital death or
live discharge. If a database export does not contain coherent hospital event
times, use binary eventual in-hospital mortality or descriptive event
proportions instead of running the Fine-Gray workflow.

License
-------

MIT License.
