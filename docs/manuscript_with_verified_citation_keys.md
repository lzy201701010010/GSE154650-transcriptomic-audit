# Sensitivity of transcriptomic inference to expression-derived histological-zone adjustment in hemorrhoidal disease versus anal fissure surgical controls: a patient-level mRNA and miRNA reanalysis

**Article type:** Research article

**Running title:** Zone-adjustment sensitivity in hemorrhoid transcriptomics

**Authors:** Zhiyao Li¹, Lili Wang², Qinghua Luo³, and Yunxiang Wu³*

**Affiliations:**

1. Department of General Surgery, Zhanjiang First Hospital of Traditional Chinese Medicine, Zhanjiang 524033, Guangdong Province, China.
2. Clinical Medical College, Jiangxi University of Chinese Medicine, Nanchang 330006, Jiangxi Province, China.
3. Department of Anorectal Surgery, Affiliated Hospital of Jiangxi University of Chinese Medicine, Nanchang 330006, Jiangxi Province, China.

**Corresponding author:** Yunxiang Wu, Department of Anorectal Surgery, Affiliated Hospital of Jiangxi University of Chinese Medicine, Nanchang 330006, Jiangxi Province, China. Email: Wyx841106@163.com. ORCID: 0009-0008-7343-5503.

**Author emails and ORCID identifiers:**

- Zhiyao Li: Lizhiyao159@outlook.com; ORCID 0009-0001-2696-6730
- Lili Wang: 2962018117@qq.com; ORCID 0009-0006-3120-1026
- Qinghua Luo: luohafo2818252021@163.com; ORCID 0000-0001-9767-9117
- Yunxiang Wu: Wyx841106@163.com; ORCID 0009-0008-7343-5503


## Abstract

### Background

Disease-group transcriptomic comparisons can depend on how tissue heterogeneity is represented in the statistical model. We examined this dependence in a public patient-level paired mRNA-miRNA Combo-Seq cohort comparing hemorrhoidal disease with anal fissure surgical controls.

### Methods

We reanalyzed paired mRNA and miRNA measurements from 38 patient-level libraries in the National Center for Biotechnology Information Gene Expression Omnibus dataset GSE154650: 20 hemorrhoidal disease and 18 anal fissure surgical controls. The original study inferred the expression-derived histological-zone classification from mRNA expression signatures, single-sample gene-set enrichment analysis and clustering. It was not an independent pathological or anatomical assessment. Limma-voom models with quality weights compared the fully unadjusted model, the sex- and BMI-adjusted model without zone, the full model including the expression-derived histological-zone classification, and the reduced-zone sensitivity analysis.

### Results

Among 12,161 genes, 6,111 met the false-discovery-rate threshold in the fully unadjusted model. In the sex- and BMI-adjusted model without zone, 5,389 met the threshold; none remained in the full model. Of these 5,389 signals, 5,255 attenuated and 1,016 reversed direction. Among 670 mature miRNAs, the corresponding model counts were 46, 45 and 0. Of the 45 no-zone signals, 41 attenuated and 5 reversed direction. The reduced-zone sensitivity analysis yielded 0 significant features in both analytes. No robust gene or adjusted-stable mature miRNA was identified. Directionally stable but low-precision effects remained for 4,210 genes and 229 mature miRNAs.

### Conclusions

Disease-effect estimates in both analytes were highly sensitive to the deposited expression-derived histological-zone classification. This sensitivity limits disease-specific inference but does not establish anatomical causation or exclude genuine low-precision disease effects.

**Keywords:** hemorrhoidal disease; anal fissure; transcriptomics; Combo-Seq; expression-derived histological-zone classification; sensitivity analysis


## Background

GSE154650 contains patient-level total-RNA Combo-Seq measurements from hemorrhoidal-disease tissue and anal fissure surgical controls. The deposited groups are distributed asymmetrically across three histological-zone classes inferred from mRNA expression. This structure makes disease-effect estimates potentially sensitive to whether the expression-derived classification is included in the model [1, 2].

The original study used keratinocyte and sebocyte expression signatures, single-sample gene-set enrichment analysis and clustering to derive the zone classes. It also adjusted a restricted candidate-gene mRNA analysis for BMI and zone. The present study is therefore not the first zone-adjusted analysis of these data. Later reanalyses emphasized differentially expressed genes, network or pathway summaries, immune-infiltration estimates and prediction-oriented models [3–5]. A genome-wide, patient-level audit of model-dependent migration across both mRNA and mature miRNA remained unreported in the reviewed sources.

We compared disease coefficients from the fully unadjusted model, the sex- and BMI-adjusted model without zone, the full model including the expression-derived histological-zone classification, and the reduced-zone sensitivity analysis. The two analytes were evaluated separately, without an mRNA-miRNA regulatory model. We quantified attenuation, direction reversal, low-precision directional stability and patient influence to define the limits of disease-specific inference.

## Methods

### Study design and public dataset

This was a secondary reanalysis of GSE154650, a public human total-RNA NEXTFLEX Combo-Seq dataset [1, 2]. The series contains 38 independent surgical-tissue libraries: 20 from patients with hemorrhoidal disease and 18 from anal fissure surgical controls. The patient/library was the statistical unit, and all 38 patients were retained in the primary analyses.

### Patient-level sample definition and common-feature universe

Official GEO records and recovered per-patient processed files were harmonized to one row per patient. Each total-RNA Combo-Seq library yielded both mRNA and miRNA measurements, enabling paired analyte analysis within the same patient. The exact valid-row intersection across all 38 patients contained 33,959 records: 32,153 ENST records and 1,806 miRNA records.

GSM4676457 had a source-truncated table with absent tail features. Those features were treated as unobserved and were neither imputed nor filled with zero. The shared 33,959-record universe therefore provided the same observed feature set for every patient. Subsequent mRNA aggregation yielded 12,161 formally tested genes. After family-specific filtering, 697 miRNA records remained, including 670 mature miRNAs in the primary formal family.

### Transcript-to-gene aggregation

Valid ENST records were mapped with the Ensembl GRCh38 release 100 annotation specified before modeling [6]. Transcript abundances mapping unambiguously to the same gene were summed before filtering. Deposited raw-abundance estimates were non-negative but often fractional. They were retained without rounding or integer coercion, so DESeq2 was not used. Genes were filtered with edgeR `filterByExpr` using the full-model design and normalized by the trimmed mean of M values method [7, 8].

### miRNA classification and filtering

The 1,806 miRNA records were classified with the miRBase 22.1 mapping specified before modeling [9]. Mature miRNAs formed the primary family; precursor and unresolved records remained separate. Mature 5p and 3p products were not merged. Filtering and TMM normalization were performed by family. The 697 retained records comprised 670 mature, 14 precursor and 13 unresolved records.

### Expression-derived histological-zone classification

The deposited `tissue_zone` variable comprised anoderm, transition and intestinal mucosa. The original study inferred it from mRNA expression signatures, single-sample gene-set enrichment analysis and clustering. Case/control counts were 2/10, 4/5 and 14/3 across these levels. Because the variable was derived from mRNA, its inclusion can introduce endogeneity or overadjustment. We evaluated sensitivity to this classification rather than treating it as an independently measured anatomical confounder.

### Statistical modeling

Gene and miRNA families were modeled separately with limma `voomWithQualityWeights`, followed by linear modeling and empirical-Bayes moderation [10–12]. The contrast was hemorrhoidal disease minus anal fissure surgical control. Sex and continuous centered BMI were prespecified covariates. Age and sequencing batch were unavailable in the recovered public metadata and were not imputed.

Four model contexts were retained:

```text
fully unadjusted model:
expression ~ disease_group

sex- and BMI-adjusted model without zone:
expression ~ disease_group + sex + BMI_centered

full model including the expression-derived histological-zone classification:
expression ~ disease_group + sex + BMI_centered + tissue_zone

reduced-zone sensitivity analysis:
exclude anoderm; retain disease_group + sex + BMI_centered + tissue_zone
```

The full design was estimable at rank 6/6, but sparse overlap limited biological comparability. The primary migration comparison was the no-zone coefficient versus the full-model coefficient. Fully unadjusted results were descriptive context. The reduced-zone sensitivity analysis included 26 patients: 18 with hemorrhoidal disease and 8 anal fissure surgical controls.

### Stability and patient influence

Migration summaries included FDR-count changes, absolute-effect attenuation, direction reversal and median absolute log2 fold change. Prespecified sensitivity analyses covered source-warning exclusions, alternative filtering, quality weights, reduced- and within-zone estimates, and leave-one-patient-out analyses. No patient was removed from the primary analysis to improve significance.

Mutually exclusive stability classes incorporated model-direction agreement, reduced-zone behavior, warning-sample sensitivity and leave-one-out influence. Robust or adjusted-stable status required full-model FDR support plus the prespecified stability criteria. Features retaining direction without full-model FDR support were labeled directionally stable low-precision. Detailed thresholds and all class counts are reported in the Supplementary Methods and tables.

### Exploratory competitive gene-set testing and multiple testing

Predefined mRNA modules were evaluated with camera as an exploratory competitive gene-set test that accounts for inter-gene correlation [13]. Patient-level module scores were analyzed separately. Camera results and patient-level score results were not treated as equivalent and are reported only in the Supplementary Results. Benjamini-Hochberg correction was applied within prespecified feature or module families, with FDR < 0.05 as the formal threshold [14].

### Software and reproducibility

Analyses were run under R 4.4.2 (2024-10-31 ucrt) [15] on Windows 11 x64 (build 26200), platform `x86_64-w64-mingw32/x64`. The locale was Chinese (Simplified), China UTF-8 for collation, character type, monetary formatting and time; `LC_NUMERIC` was `C`. The time zone was Asia/Shanghai.

Package versions were limma 3.62.2, edgeR 4.4.2, statmod 1.5.2, data.table 1.18.4, ggplot2 4.0.2 and digest 0.6.37. Scripts, session information and reproducibility records are maintained in the author-controlled repository. No models were refitted, FDR values recalculated, camera tests rerun, mRNA and miRNA integrated, targets predicted, networks constructed, or external or single-cell analyses added during manuscript revision.

### Use of AI-assisted tools

ChatGPT and OpenAI Codex were used to assist with language refinement, code checking and debugging, and workflow and document organization. Study conceptualization, methodological decisions, statistical interpretation, and scientific conclusions were determined by the authors. All AI-assisted outputs were independently reviewed and verified by the authors, who take full responsibility for the final work.

## Results

### Cohort structure, paired Combo-Seq profiles, and expression-derived histological-zone classification

The analysis included 38 patient-level Combo-Seq libraries: 20 hemorrhoidal disease and 18 anal fissure surgical controls. The exact common-feature universe contained 33,959 records across all patients, comprising 32,153 ENST and 1,806 miRNA records. GSM4676457 lacked tail features because of source truncation; absent records were not zero-filled. The common universe ensured that the primary analysis compared the same observed features across patients.

The ENST branch was aggregated to 12,161 formally tested genes. MiRNA filtering retained 697 records, including 670 mature miRNAs in the primary family, 14 precursors and 13 unresolved records (Figure 1). The expression-derived zone distribution was asymmetric: case/control counts were 2/10 for anoderm, 4/5 for transition and 14/3 for intestinal mucosa. The classification came from mRNA signatures, single-sample gene-set enrichment and clustering rather than independent pathological or anatomical assessment. Table 1 summarizes the cohort, covariates and model contexts.

### mRNA disease-effect estimates were highly sensitive to adjustment for the expression-derived histological-zone classification

All 12,161 genes entered formal testing. In the fully unadjusted model, 6,111 genes met Benjamini-Hochberg FDR < 0.05. After adjustment for sex and centered BMI without zone, 5,389 met the threshold. Adding the expression-derived classification reduced this direct comparison from 5,389 to 0. The reduced-zone sensitivity analysis also yielded 0 FDR-significant genes (Figure 2a).

Among the 5,389 no-zone signals, 5,255 (97.5%) attenuated in absolute magnitude. Median absolute log2 fold change declined from 0.513 to 0.149, and 1,016 (18.9%) reversed direction (Figure 2b,c). No gene met the prespecified robust definition. A further 4,210 genes were directionally stable but low precision. Sample-sensitive and zone-dependent classes, together with leave-one-patient-out changes, limited feature-level interpretation (Figure 2d). Exploratory module-level analyses are reported in the Supplementary Results.

### Mature-miRNA estimates showed a parallel sensitivity pattern

Of 1,806 deposited miRNA records, 697 passed the prespecified filters. The formal mature-miRNA denominator was 670. In this family, 46 miRNAs met FDR < 0.05 in the fully unadjusted model. Adjustment for sex and centered BMI without zone yielded 45 signals. Adding the expression-derived classification reduced the direct comparison from 45 to 0, and the reduced-zone sensitivity analysis also yielded 0 (Figure 3a).

Of the 45 no-zone signals, 41 (91.1%) attenuated. Median absolute log2 fold change declined from 1.767 to 0.537, and 5 (11.1%) reversed direction (Figure 3b,c). No adjusted-stable mature miRNA was identified. The mature-family classification contained 229 directionally stable low-precision, 333 sample-sensitive and 93 zone-dependent miRNAs (Figure 3d). The all-record supplementary counts were 348 sample-sensitive and 95 zone-dependent among 697 records. These were not substituted for the mature-family counts.

### Cross-analyte auditing defined the limits of disease-specific molecular inference

Both analytes showed broad FDR signal sets in models without zone and no single-feature FDR signal after the full model or reduced-zone sensitivity analysis (Figure 4a). Most no-zone signals attenuated, some reversed direction, and both analytes retained low-precision directional subsets with material patient sensitivity (Figure 4b-d).

Table 2 provides the cross-model summary, and Supplementary Table S1 records the claim-evidence boundary. The data support inferential sensitivity to the expression-derived classification. They do not establish a causal anatomical-zone effect, prove that zone caused all no-zone differences, or exclude genuine low-precision disease effects.

## Discussion

This reanalysis shows that disease-effect estimates in GSE154650 depend strongly on whether the deposited expression-derived histological-zone classification is included in the model. Both analytes lost all FDR-supported single-feature results in the full model and reduced-zone sensitivity analysis. Most no-zone estimates attenuated, some reversed direction, and residual directional subsets remained low precision. The contribution is an audit of inference boundaries, not a mechanism or molecular signature.

The asymmetric disease-by-zone distribution provides one plausible source of the broad no-zone contrasts. A disease coefficient may incorporate both disease-associated biology and differences in sampled tissue context. This does not mean that zone caused every apparent difference. Disease status, tissue composition, technical variation and unreported factors may contribute together, and the cohort cannot fully separate them.

The zone classification may capture meaningful histological or compositional variation, but it was derived from mRNA expression. Adjusting mRNA for that variable can absorb tissue heterogeneity and genuine disease-associated signal. The full model should be interpreted as one sensitivity specification rather than as a uniquely correct biological model. The analysis establishes model sensitivity, not an independently measured anatomical cause or complete removal of confounding.

Mature miRNAs showed a similar migration even though they were not the primary expression layer used to construct the zone classes. This parallel result strengthens the sensitivity observation. It is not fully independent validation because both analytes came from the same patients, specimens and libraries. They share the same sampling structure and influential observations.

The original Gut study derived the zone classes and adjusted a restricted candidate-gene mRNA analysis for BMI and zone [2]. Our incremental contribution is the genome-wide audit of mRNA and mature-miRNA migration, attenuation, reversal, low-precision stability and patient influence. Later GSE154650 reanalyses addressed CALM3, mitochondria-associated membrane genes or glycolysis-related programs [3–5]. These studies asked different questions. The present findings call for caution with disease-specific interpretations from no-zone models but do not invalidate every observation in those reports.

Future studies should record sampling subregions independently and, where feasible, confirm them by pathology. Case and comparator tissues should be balanced or matched at the sampling level. Adjustment variables should be prespecified rather than derived from the expression layer under test. Larger cohorts could support adequately powered within-zone estimates, interactions and independent replication. Age, batch and other technical variables should also be recorded explicitly.

Several limitations define the conclusion. This was one public cohort of 38 patients, with anal fissure surgical controls rather than healthy tissue. The zone variable was mRNA-derived, no independent pathological assessment of the classification was available, and mRNA models may be overadjusted. The miRNA layer was not fully independent. Age and batch were unavailable. Zone overlap was sparse, and the reduced-zone sensitivity analysis included 26 patients. Source values included fractional abundances and file-format warnings. Finally, 0 adjusted FDR features do not imply biological equivalence or absence of true disease effects. No sampling-matched external transcriptomic cohort was available.

## Conclusions

Disease-effect estimates in this patient-level Combo-Seq cohort were highly sensitive to adjustment for the expression-derived histological-zone classification. The cross-analyte migration and patient sensitivity limit disease-specific inference while leaving open genuine effects that the available design cannot separate from tissue context.

## List of abbreviations

- BMI: body mass index
- ENST: Ensembl transcript identifier
- FDR: false discovery rate
- GEO: Gene Expression Omnibus
- mRNA: messenger RNA
- miRNA: microRNA
- ssGSEA: single-sample gene-set enrichment analysis
- TMM: trimmed mean of M values

## Declarations

### Ethics approval and consent to participate

This study exclusively involved the secondary analysis of publicly available, de-identified data obtained from the Gene Expression Omnibus under accession number GSE154650. The authors did not recruit participants, access identifiable personal information, or collect additional biological specimens. Under the institutional requirements applicable to this secondary analysis, additional ethics committee review was not required. No ethics approval or waiver reference number was issued.

Consent to participate: Not applicable because this study involved only the secondary analysis of publicly available, de-identified data.

### Consent for publication

Not applicable.

### Availability of data and materials

The dataset analyzed in this study is publicly available in the National Center for Biotechnology Information Gene Expression Omnibus under accession number GSE154650. The analysis scripts, figure source data, derived result tables, and reproducibility documentation are currently maintained in a private author-controlled repository and will be made publicly accessible before manuscript submission. Permanent archive information will be added after repository release and verification.

### Competing interests

The authors declare that they have no competing interests.

### Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

### Authors' contributions

ZL contributed to conceptualization, methodology, software, validation, formal analysis, data curation, visualization, writing of the original draft, and review and editing of the manuscript. LW and QL contributed to validation and review and editing of the manuscript. YW contributed to supervision, project administration, and review and editing of the manuscript. All authors read and approved the final manuscript.

### Acknowledgements

The authors acknowledge the investigators who generated and deposited the GSE154650 dataset in the Gene Expression Omnibus.


## References

1. Juzenas S, Franke A. Transcriptome Analysis of Human Hemorrhoids Tissue. *NCBI Gene Expression Omnibus*. 2021. GSE154650; SRP272387. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154650. Accessed 17 Jul 2026.

2. Zheng T, Ellinghaus D, Juzenas S, et al. Genome-wide analysis of 944 133 individuals provides insights into the etiology of haemorrhoidal disease. *Gut*. 2021;70(8):1538-1549. doi:10.1136/gutjnl-2020-323868; PMID:33888516.

3. He J, Ni Z, Li Z. CALM3 affects the prognosis of leukemia and hemorrhoids. *Medicine (Baltimore)*. 2023;102(44):e36027. doi:10.1097/MD.0000000000036027; PMID:37932969.

4. Mao L, Rao Z, Wang Y, Yang J, He J, Zheng Z, Chen L. Identification and Validation of Key Genes Involved in the Coupling of Mitochondria-Associated Endoplasmic Reticulum Membrane in Hemorrhoidal Disease. *Int J Gen Med*. 2025;18:2781-2798. doi:10.2147/IJGM.S511281; PMID:40469970.

5. Li P, Hou Q, Yang X, Han W, Wang H. Role of glycolysis related genes in the pathogenesis of hemorrhoids and immune cell infiltration analysis. *Sci Rep*. 2025;15(1):32912. doi:10.1038/s41598-025-18382-3; PMID:40998906.

6. Yates AD, Achuthan P, Akanni W, et al. Ensembl 2020. *Nucleic Acids Res*. 2020;48(D1):D682-D688. doi:10.1093/nar/gkz966; PMID:31691826.

7. Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. *Nucleic Acids Res*. 2025;53(2):gkaf018. doi:10.1093/nar/gkaf018; PMID:39844453.

8. Chen Y, Lun AT, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Res*. 2016;5:1438. doi:10.12688/f1000research.8987.2; PMID:27508061.

9. Kozomara A, Birgaoanu M, Griffiths-Jones S. miRBase: from microRNA sequences to function. *Nucleic Acids Res*. 2019;47(D1):D155-D162. doi:10.1093/nar/gky1141; PMID:30423142.

10. Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol*. 2014;15(2):R29. doi:10.1186/gb-2014-15-2-r29; PMID:24485249.

11. Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat ML, Smyth GK, Ritchie ME. Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. *Nucleic Acids Res*. 2015;43(15):e97. doi:10.1093/nar/gkv412; PMID:25925576.

12. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47. doi:10.1093/nar/gkv007; PMID:25605792.

13. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res*. 2012;40(17):e133. doi:10.1093/nar/gks461; PMID:22638577.

14. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Stat Methodol*. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.

15. R Core Team. R: A Language and Environment for Statistical Computing. *R Foundation for Statistical Computing*. 2024. R version 4.4.2. https://www.r-project.org/.

## Figure legends

### Figure 1 | Cohort, expression-derived histological-zone classification provenance and model framework

**a,** Design for 38 libraries: 20 hemorrhoidal disease and 18 anal fissure surgical controls. **b,** Classification provenance from expression signatures, single-sample gene-set enrichment and clustering. **c,** Case/control classification counts: 2/10, 4/5 and 14/3. **d,** Four prespecified model contexts, including the reduced-zone sensitivity analysis with n = 26. The classification is expression-derived, not independently measured anatomy.

### Figure 2 | Migration and stability of mRNA disease-effect estimates

**a,** FDR counts across models; the isolated comparison after adding the expression-derived classification is 5,389 to 0, while 6,111 is the fully unadjusted context. **b,** Coefficients from the sex- and BMI-adjusted model without zone versus the full model for 12,161 genes. **c,** Among 5,389 signals, 5,255 attenuated and 1,016 reversed; median absolute log2 fold change changed from 0.513 to 0.149. **d,** Stability and patient influence. No robust gene was identified; 4,210 directional signals were low precision.

### Figure 3 | Migration and stability of mature-miRNA disease-effect estimates

**a,** Filtering retained 697 of 1,806 records, including 670 mature miRNAs; the mature-family ladder was 46, 45, 0 and 0. **b,** Coefficients from the sex- and BMI-adjusted model without zone versus the full model. **c,** Among 45 signals, 41 attenuated and 5 reversed; median absolute log2 fold change changed from 1.767 to 0.537. **d,** Mature-family stability and patient influence. Low-precision directional signals numbered 229; adjusted-stable signals numbered 0.

### Figure 4 | Cross-analyte inference boundary

**a,** Parallel model migration. **b,** Within-analyte attenuation and reversal proportions. **c,** Prespecified evidence-class composition. **d,** A non-causal inference boundary. The paired analytes are not independent replications, and null FDR results in the full model do not prove equivalence.


## Table legends

### Table 1 | Cohort, covariates and model contexts

Values are patient counts unless stated otherwise. The comparator is anal fissure surgical controls, not healthy tissue. The deposited expression-derived histological-zone classification was inferred in the original study from mRNA expression signatures, ssGSEA and clustering and is not independently pathology confirmed. BMI was modeled as a continuous centered covariate; descriptive P values were not used for covariate selection. Age and batch were unavailable in the recovered public metadata. The reduced-zone sensitivity analysis excludes anoderm and retains 18 hemorrhoidal-disease and 8 control patients.

### Table 2 | Cross-model mRNA and mature-miRNA summary

FDR counts use Benjamini-Hochberg correction within the prespecified gene or mature-miRNA family. Fully unadjusted model counts (6,111 genes; 46 mature miRNAs) are descriptive context. The isolated addition of the expression-derived histological-zone classification compares the sex- and BMI-adjusted model without zone with the full model (5,389 to 0 genes; 45 to 0 mature miRNAs). Stability counts are mutually exclusive within each analyte; the mature-miRNA row excludes precursor and unresolved records. Directionally stable low-precision signals are not robust biomarkers, and parallel analyte sensitivity does not imply regulatory integration.

## Additional file

### Additional file 1

**Format:** Multi-file supplementary package to be assembled as editable tables plus PDF figures during the submission-package stage.

**Title:** Supplementary methods, results, figures and tables for the GSE154650 transcriptomic sensitivity audit.

**Description:** The package contains Supplementary Methods and Results, Supplementary Figures S1–S12, Supplementary Tables S1–S10b, patient/sample records, complete model results, filtering and mapping records, stability and influence audits, software versions and reproducibility records. Supplementary Table S1 is the prespecified claim-evidence boundary. Camera competitive results are supplementary only; patient-level module-score comparisons did not pass FDR correction. Unadjusted feature lists are not biomarker tables. No target network, mechanism figure, diagnostic model or single-cell result is included.
