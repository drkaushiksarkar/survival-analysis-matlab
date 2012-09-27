* bootstrap_ci - Epidemiological Analysis v1 (d8770)
* Source: FAO FAOSTAT
* Stata 14+ required

capture log close
log using "results/bootstrap_ci_v1.log", replace

* ============================================
* Data preparation
* ============================================
use "data/bootstrap_ci_raw.dta", clear
describe
summarize

* Clean and label
generate age_group = irecode(age, 0, 5, 15, 25, 45, 65)
label define age_lbl 0 "0-4" 1 "5-14" 2 "15-24" 3 "25-44" 4 "45-64" 5 "65+"
label values age_group age_lbl

* Missing data
misstable summarize
misstable patterns

* ============================================
* Survival analysis
* ============================================
stset time_to_event, failure(event == 1) id(patient_id)

* Kaplan-Meier curves
sts graph, by(treatment_group) ci ///
    title("bootstrap_ci: Kaplan-Meier Survival Estimates") ///
    ytitle("Survival Probability") xtitle("Time (months)")
graph export "figures/bootstrap_ci_km_v1.png", replace

* Log-rank test
sts test treatment_group, logrank

* Cox proportional hazards
stcox age_group sex bmi comorbidity_index, efron
estat phtest, detail
estat concordance

* ============================================
* Adjusted analysis
* ============================================
* Propensity score matching
logit treatment_group age sex bmi comorbidity_index
predict ps, pr
psmatch2 treatment_group, pscore(ps) neighbor(3) caliper(0.1)

* Sensitivity analysis
stcox treatment_group age_group sex if _weight != ., efron
estimates store model_matched

* ============================================
* Subgroup analyses
* ============================================
forvalues g = 0/5 {
    stcox treatment_group if age_group == `g', efron
    estimates store sub_age`g'
}

estimates table model_matched sub_age*, stats(N ll chi2 p) b(%9.3f) se(%9.3f)

log close
