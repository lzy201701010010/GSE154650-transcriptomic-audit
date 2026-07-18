# Table 1 | Cohort, covariates and model contexts

Values are patient counts unless stated otherwise. The comparator is anal fissure surgical controls, not healthy tissue. The deposited expression-derived histological-zone classification was inferred in the original study from mRNA expression signatures, ssGSEA and clustering and is not independently pathology confirmed. BMI was modeled as a continuous centered covariate; descriptive P values were not used for covariate selection. Age and batch were unavailable in the recovered public metadata. The reduced-zone sensitivity analysis excludes anoderm and retains 18 hemorrhoidal-disease and 8 control patients.

| variable | level_or_statistic | hemorrhoidal_disease | anal_fissure_surgical_controls | missingness | interpretive_note |
| --- | --- | --- | --- | --- | --- |
| Patients | Total n | 20 | 18 | 0 | Patient/library is the statistical unit. |
| Sex | Female, n | 8 | 9 | 0 | Descriptive; no P value used for covariate selection. |
| Sex | Male, n | 12 | 9 | 0 | Descriptive; no P value used for covariate selection. |
| BMI | Mean ± SD | 27.54 ± 6.80 | 28.93 ± 5.08 | 0 | Continuous centered BMI was used in adjusted models. |
| Expression-derived histological-zone classification | Anoderm, n | 2 | 10 | 0 | Expression-derived; not independently pathology confirmed. |
| Expression-derived histological-zone classification | Transition, n | 4 | 5 | 0 | Expression-derived; not independently pathology confirmed. |
| Expression-derived histological-zone classification | Intestinal mucosa, n | 14 | 3 | 0 | Expression-derived; not independently pathology confirmed. |
| Paired analyte completeness | mRNA available, n | 20 | 18 | 0 | Same patient Combo-Seq library. |
| Paired analyte completeness | miRNA available, n | 20 | 18 | 0 | Same patient Combo-Seq library. |
| Age | Not reported | Unavailable | Unavailable | 38 | Not available in recovered public metadata. |
| Batch | Not reported | Unavailable | Unavailable | 38 | Not available in recovered public metadata. |
| Model hierarchy | Fully unadjusted model | n = 38 | n = 38 | Not applicable | expression ~ disease group |
| Model hierarchy | Sex- and BMI-adjusted model without zone | n = 38 | n = 38 | Not applicable | expression ~ disease group + sex + centered BMI |
| Model hierarchy | Full model including the expression-derived histological-zone classification | n = 38 | n = 38 | Not applicable | expression ~ disease group + sex + centered BMI + expression-derived histological-zone classification |
| Model hierarchy | Reduced-zone sensitivity analysis | n = 18 after excluding anoderm | n = 8 after excluding anoderm | Not applicable | Anoderm excluded; transition and intestinal-mucosa classes retained with zone adjustment. |
