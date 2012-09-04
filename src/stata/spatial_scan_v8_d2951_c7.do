* spatial_scan - Disease Burden Analysis v8 (d2951)
* Source: World Bank Health Nutrition Population
* Population-level mortality and morbidity analysis

capture log close
log using "results/spatial_scan_burden_v8.log", replace

* ============================================
* Mortality rate computation
* ============================================
use "data/spatial_scan_mortality.dta", clear

* Age-standardized rates (WHO standard population)
generate who_weight = .
replace who_weight = 0.0886 if age_group == 1
replace who_weight = 0.0869 if age_group == 2
replace who_weight = 0.0860 if age_group == 3
replace who_weight = 0.0858 if age_group == 4
replace who_weight = 0.0822 if age_group == 5
replace who_weight = 0.0635 if age_group == 6
replace who_weight = 0.0370 if age_group == 7

* Crude rate
generate crude_rate = (deaths / population) * 100000

* Standardized rate per country-year
bysort country year: egen asr = total(crude_rate * who_weight)

* ============================================
* Trend analysis
* ============================================
* Join-point regression approximation
bysort country (year): generate time = _n
xtset country_id time

* Fixed effects panel model for trends
xtreg asr time, fe
predict asr_hat, xb
predict residual, e

* Test for non-linearity
generate time2 = time^2
xtreg asr time time2, fe
testparm time2

* ============================================
* Years of Life Lost (YLL)
* ============================================
generate yll = deaths * (life_expectancy - age_midpoint)
replace yll = 0 if yll < 0

* Discounting (3% rate, as per GBD methodology)
generate yll_discounted = yll * (1 / (1.03^time))

bysort country year: egen total_yll = total(yll_discounted)

* Export results
export delimited using "output/spatial_scan_burden_v8.csv", replace

log close
