# A Note on Guidelines for Structural and Copy-Number Variants

Last updated: 2026-11-18

## Summary

There is no ACMG guideline for structural variants in general. ACMG covers them through two copy-number standards in *Genetics in Medicine*

- **Germline (constitutional) CNVs:** Riggs et al. 2020 (ACMG + ClinGen) [1].
- **Tumor (acquired) CNAs, CN-LOH, and acquired SVs:** Mikhail et al. 2019 (ACMG + Cancer Genomics Consortium) [3], with AMP/ASCO/CAP-style clinical tiers [4]. Mikhail 2019 states that its recommendations should also apply to acquired structural variants detected by sequencing.

## Relevance to the annotation pipeline

The SV/CNA filter and Tier 3 in `2026-04-30_WGS_Nanopore_Annotation_T-ALL` use AnnotSV's `ACMG_class`. That class applies Riggs 2020, the germline dosage-pathogenicity framework [5].

For the tumor-only somatic calls in this cohort:

- cite Mikhail 2019 [3] as the somatic standard;
- describe the AnnotSV `ACMG_class` as a pathogenicity score, not an ACMG somatic tier.

## References

1. Riggs ER, Andersen EF, Cherry AM, et al. Technical standards for the interpretation and reporting of constitutional copy-number variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics (ACMG) and the Clinical Genome Resource (ClinGen). *Genet Med*. 2020;22(2):245–257. doi:[10.1038/s41436-019-0686-8](https://doi.org/10.1038/s41436-019-0686-8). [PMC7313390](https://pmc.ncbi.nlm.nih.gov/articles/PMC7313390)
2. Correction to [1]. *Genet Med*. 2021;23(11):2230. doi:[10.1038/s41436-021-01150-9](https://doi.org/10.1038/s41436-021-01150-9)
3. Mikhail FM, Biegel JA, Cooley LD, et al. Technical laboratory standards for interpretation and reporting of acquired copy-number abnormalities and copy-neutral loss of heterozygosity in neoplasia. *Genet Med*. 2019;21(9):1903–1916. doi:[10.1038/s41436-019-0545-7](https://doi.org/10.1038/s41436-019-0545-7)
4. Li MM, Datto M, Duncavage EJ, et al. Standards and guidelines for the interpretation and reporting of sequence variants in cancer: a joint consensus recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. *J Mol Diagn*. 2017;19(1):4–23. doi:[10.1016/j.jmoldx.2016.10.002](https://doi.org/10.1016/j.jmoldx.2016.10.002)
5. AnnotSV SV ranking and classification documentation (adaptation of Riggs et al. 2020): <https://lbgi.fr/AnnotSV/ranking>
