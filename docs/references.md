# References

Curated literature, consortium SOPs, and tool documentation that justify the
design choices in [`pipeline_design.md`](pipeline_design.md). Citation keys
match the `[Ref N]` markers used there.

## Reference panels and imputation servers

1. Taliun, D., Harris, D. N., Kessler, M. D., *et al.* **Sequencing of
   53,831 diverse genomes from the NHLBI TOPMed Program.** *Nature* 590,
   290–299 (2021). https://doi.org/10.1038/s41586-021-03205-y
2. Das, S., Forer, L., Schönherr, S., *et al.* **Next-generation genotype
   imputation service and methods.** *Nature Genetics* 48, 1284–1287
   (2016). https://doi.org/10.1038/ng.3656
3. TOPMed Imputation Server v2 documentation, "Input data preparation" and
   "Frequently encountered errors." https://statgen.github.io/tis-v2-docs/
   (accessed 2026-05).

## QC best practices

4. Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R.,
   Morris, A. P., & Zondervan, K. T. **Data quality control in genetic
   case-control association studies.** *Nature Protocols* 5, 1564–1573
   (2010). https://doi.org/10.1038/nprot.2010.116
5. Marees, A. T., de Kluiver, H., Stringer, S., *et al.* **A tutorial on
   conducting genome-wide association studies: Quality control and
   statistical analysis.** *International Journal of Methods in Psychiatric
   Research* 27, e1608 (2018). https://doi.org/10.1002/mpr.1608
6. Verma, S. S., de Andrade, M., Tromp, G., *et al.* **Imputation and
   quality control steps for combining multiple genome-wide datasets.**
   *Frontiers in Genetics* 5, 370 (2014).
   https://doi.org/10.3389/fgene.2014.00370
7. Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., &
   Chen, W.-M. **Robust relationship inference in genome-wide association
   studies.** *Bioinformatics* 26, 2867–2873 (2010).
   https://doi.org/10.1093/bioinformatics/btq559
8. Karczewski, K. J., Francioli, L. C., Tiao, G., *et al.* **The mutational
   constraint spectrum quantified from variation in 141,456 humans.**
   *Nature* 581, 434–443 (2020). (gnomAD ancestry-assignment pipeline.)
   https://doi.org/10.1038/s41586-020-2308-7

## Harmonization and tooling

9. Rayner, W. **HRC or 1000G Imputation Preparation and Checking Tool.**
   https://www.chg.ox.ac.uk/~wrayner/tools/ (TOPMed-adapted variant of
   `HRC-1000G-check-bim.pl`).
10. Zhan, X. **checkVCF.py — sanity-check VCFs prior to imputation.**
    https://github.com/zhanxw/checkVCF
11. McCarthy, S., Das, S., Kretzschmar, W., *et al.* **A reference panel of
    64,976 haplotypes for genotype imputation.** *Nature Genetics* 48,
    1279–1283 (2016). https://doi.org/10.1038/ng.3643 (Methods for
    INFO/R² interpretation and leave-one-out concordance.)

## Tools (versions to be pinned in the container)

- **PLINK 2.0** — Chang, C. C., Chow, C. C., Tellier, L. C. A. M.,
  Vattikuti, S., Purcell, S. M., & Lee, J. J. *Second-generation PLINK:
  rising to the challenge of larger and richer datasets.* *GigaScience* 4,
  7 (2015). https://doi.org/10.1186/s13742-015-0047-8
- **bcftools / samtools** — Danecek, P., Bonfield, J. K., Liddle, J.,
  *et al.* *Twelve years of SAMtools and BCFtools.* *GigaScience* 10,
  giab008 (2021). https://doi.org/10.1093/gigascience/giab008
- **Eagle v2** — Loh, P.-R., Danecek, P., Palamara, P. F., *et al.*
  *Reference-based phasing using the Haplotype Reference Consortium panel.*
  *Nature Genetics* 48, 1443–1448 (2016).
  https://doi.org/10.1038/ng.3679
- **Minimac4** — https://github.com/statgen/Minimac4
- **CrossMap** — Zhao, H., Sun, Z., Wang, J., Huang, H., Kocher, J.-P., &
  Wang, L. *CrossMap: a versatile tool for coordinate conversion between
  genome assemblies.* *Bioinformatics* 30, 1006–1007 (2014).
  https://doi.org/10.1093/bioinformatics/btt730

## Consortium SOPs consulted

- NHGRI-EBI GWAS Catalog imputation guidance.
- Pan-UK Biobank QC and imputation pipeline notes.
- TOPMed Analysis Working Group recommendations on multi-ancestry QC.
