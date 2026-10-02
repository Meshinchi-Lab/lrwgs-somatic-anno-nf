# Meshinchi-Lab/lrwgs-somatic-anno-nf: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

> Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [Severus](https://doi.org/10.1038/s41587-025-02618-8) — long-read somatic SV caller

  > Keskus, A.G., Bryant, A., Ahmad, T. et al. Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing. Nat Biotechnol 44, 247–257 (2026). doi: 10.1038/s41587-025-02618-8.

- [SAVANA](https://doi.org/10.1038/s41592-025-02708-0) — long-read somatic SV / CNA caller (secondary SV caller in this pipeline)

  > Elrick, H., Sauer, C.M., Espejo Valle-Inclan, J. et al. SAVANA: reliable analysis of somatic structural variants and copy number aberrations using long-read sequencing. Nat Methods 22, 1436–1446 (2025). doi: 10.1038/s41592-025-02708-0.

- [Wakhan](https://doi.org/10.64898/2025.12.11.25342098) — chromosome-scale CNA caller for long-read tumor genomes

  > Ahmad, T., Keskus, A.G., Aganezov, S., Goretsky, A., Rodriguez, I., Yoo, B., Lansdon, L.A., Repnikova, E.A., Zhang, L., Liu, Y., Donmez, A., Bryant, A., Tulsyan, S., Park, J., Gardner, J., McNulty, B., Sacco, S., Shetty, J., Zhao, Y., Tran, B., Malikic, S., Day, C.-P., Miga, K., Paten, B., Sahinalp, C., Farooqi, M.S., Dean, M., Kolmogorov, M. Wakhan: reconstruction of chromosome-scale copy number profiles of tumor genomes with long-read sequencing. medRxiv 2025.12.11.25342098 (2025). doi: 10.64898/2025.12.11.25342098.

- [ClairS](https://doi.org/10.1101/2023.08.17.553778) — long-read deep-learning somatic small-variant caller (ClairS-TO is the tumor-only branch used in this pipeline)

  > Zheng, Z., Su, J., Chen, L., Lee, Y.-L., Lam, T.-W., Luo, R. ClairS: a deep-learning method for long-read somatic small variant calling. bioRxiv 2023.08.17.553778 (2023). doi: 10.1101/2023.08.17.553778.

- [DeepSomatic](https://doi.org/10.1038/s41587-025-02839-x) — multi-technology deep-learning somatic small-variant caller

  > Park, J., Cook, D.E., Chang, P.-C. et al. Accurate somatic small variant discovery for multiple sequencing technologies with DeepSomatic. Nat Biotechnol (2025). doi: 10.1038/s41587-025-02839-x.

- [OncoKB](https://doi.org/10.1158/2159-8290.CD-23-0467) — clinical actionability knowledge base (source of SNV oncogenicity + therapy annotations)

  > Suehnholz, S.P., Nissan, M.H., Zhang, H., Kundra, R., Nandakumar, S., Lu, C., Carrero, S., Dhaneshwar, A., Fernandez, N., Xu, B.W., Arcila, M.E., Zehir, A., Syed, A., Brannon, A.R., Rudolph, J.E., Paraiso, E., Sabbatini, P.J., Levine, R.L., Dogan, A., Gao, J., Ladanyi, M., Drilon, A., Berger, M.F., Solit, D.B., Schultz, N., Chakravarty, D. Quantifying the Expanding Landscape of Clinical Actionability for Patients with Cancer. Cancer Discov 14(1), 49–65 (2024). doi: 10.1158/2159-8290.CD-23-0467.

- [T-ALL genomic reference](https://doi.org/10.1038/s41586-024-07807-0) — source cohort for the ST17 recurrent SV BED and ST19 recurrent CNV BED used in driver-tier assignment (`assign_driver_tier()` / `assign_driver_tier_cna()`)

  > Pölönen, P., Di Giacomo, D., Seffernick, A.E. et al. The genomic basis of childhood T-lineage acute lymphoblastic leukaemia. Nature 632, 1082–1091 (2024). doi: 10.1038/s41586-024-07807-0.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
