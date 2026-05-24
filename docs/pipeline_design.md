# GWAS Array Imputation Preparation Pipeline — Design

This document describes the design of a state-of-the-art preparation pipeline
for genome-wide genotyping array data prior to imputation against the
**TOPMed r3** reference panel on the
[TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/)
([TIS v2 docs](https://statgen.github.io/tis-v2-docs/),
[server code](https://github.com/statgen/imputationserver2)).

It is the design deliverable for the issue
*"Design and Document State-of-the-Art GWAS Array Imputation Preparation
Pipeline"* and is intended to be the authoritative reference that the
forthcoming scripts in this repository implement.

References cited as `[Ref N]` are listed in
[`references.md`](references.md).

---

## 1. Scope, assumptions, and design principles

**In scope.** Sample- and variant-level QC, harmonization, and packaging of
**genotyping-array** data (typically 0.3–2 M markers) into per-chromosome
VCFs suitable for the TOPMed Imputation Server (GRCh38). Multi-ancestry
cohorts are explicitly supported.

**Out of scope.** Whole-genome / whole-exome sequencing input; local
imputation (we delegate phasing and imputation to the TOPMed server, which
uses Eagle v2.4 + Minimac4 against ~133 k TOPMed samples
[Ref 1, Ref 2]); post-imputation association testing.

**Design principles.**

1. **Reproducibility.** Every step runs inside a versioned Docker image; all
   tool versions and parameters are recorded in `logs/provenance.json`.
2. **Maximize variant yield, minimize wasted server runs.** The TOPMed
   server rejects entire chromosomes for any of: wrong build, non-ACGT
   alleles, unsorted records, ref-mismatches above its threshold, duplicate
   positions, or `chr` prefix mismatch with the reference [Ref 3]. We
   eliminate all of these *before* upload.
3. **Ancestry-agnostic QC.** Per-population thresholds (HWE, MAF) are
   applied **within inferred ancestry groups**, not globally, to avoid
   discarding variants that are merely differentiated across populations
   [Ref 4, Ref 5].
4. **Conservative palindromic handling.** A/T and C/G SNPs with
   intermediate allele frequencies are flagged and excluded from
   strand-resolution heuristics; the server's reference-based alignment is
   the source of truth [Ref 3, Ref 6].
5. **Auditability.** Every variant or sample dropped is logged with the
   step, rule, and value that triggered exclusion.

---

## 2. Inputs

| Input | Format | Notes |
| --- | --- | --- |
| Genotypes | PLINK1 `bed/bim/fam`, PLINK2 `pgen/pvar/psam`, or VCF/BCF | **Must be on GRCh38.** If on GRCh37, run the documented `liftOver` recipe (Section 10) first. |
| Reference FASTA | `GRCh38_full_analysis_set_plus_decoy_hla.fa(.fai)` | Same contigs/naming as the TOPMed panel; `chr`-prefixed contigs. |
| Sample metadata (optional) | TSV with `IID`, `sex`, `self_reported_ancestry`, `phenotype`, `batch` | Used for stratified QC and `TableOne` reporting. |
| Array manifest (optional) | Illumina/Affymetrix CSV | Used to flag known problematic probes and confirm strand. |

Build is auto-detected by comparing a random sample of variant
positions/alleles against the supplied FASTA; mismatched build aborts the
run with a clear error.

---

## 3. High-level workflow

```
            ┌───────────────────────────────────────────────────┐
            │ 0. Provenance & input inventory                   │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 1. Convert to VCF (PLINK2) and normalize          │
            │    bcftools norm -f REF -c ws --multiallelics -any│
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 2. Variant QC                                     │
            │    - drop non-ACGT, indels (arrays), dups         │
            │    - missingness > 0.05                           │
            │    - MAF < 0.0001 (server threshold) [Ref 3]      │
            │    - HWE p < 1e-6 within ancestry [Ref 4]         │
            │    - flag palindromic A/T, C/G w/ 0.40<MAF<0.60   │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 3. Sample QC                                      │
            │    - call-rate < 0.97                             │
            │    - |F_het| > 0.20 outliers                      │
            │    - sex check (chrX heterozygosity)              │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 4. Relatedness (KING) + Ancestry (PCA)            │
            │    - flag/remove 2nd-degree+ relatives (kinship   │
            │      > 0.0884) per analysis plan [Ref 7]          │
            │    - project onto 1000G+HGDP PCs for ancestry     │
            │      assignment [Ref 8]                           │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 5. TOPMed harmonization                           │
            │    - bcftools +fixref or Will Rayner's            │
            │      HRC-1000G-check-bim adapted for TOPMed       │
            │      [Ref 9]                                      │
            │    - flip strand where unambiguous                │
            │    - drop unresolvable ref-mismatches             │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 6. Per-chromosome VCFs                            │
            │    - split by chromosome (chr1..chr22, chrX)      │
            │    - bcftools sort | bgzip | tabix                │
            │    - assert <chr prefix> matches reference        │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 7. Pre-submission validation                      │
            │    - checkVCF.py vs GRCh38 reference [Ref 10]     │
            │    - imputationbot validate (optional)            │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 8. Submit to TOPMed Imputation Server             │
            │    panel=apps@topmed-r3, build=hg38,              │
            │    phasing=eagle, population=mixed,               │
            │    r2-filter=off (filter post-hoc)                │
            └───────────────┬───────────────────────────────────┘
                            ▼
            ┌───────────────────────────────────────────────────┐
            │ 9. Post-imputation QC (separate workflow)         │
            │    - INFO/R2 distribution, MAF concordance,       │
            │      genotype concordance on held-out variants    │
            └───────────────────────────────────────────────────┘
```

The corresponding Mermaid diagram in the top-level `README.md` is the
canonical rendered version.

---

## 4. Step-by-step rationale and commands

The commands below illustrate the operations the pipeline will perform.
Concrete script wrappers, logging, and parallelization will be added in
follow-up PRs.

### Step 0 — Provenance and input inventory

- Record tool versions (`plink2 --version`, `bcftools --version`,
  `king --version`) and SHA256 of all input files.
- Detect build by spot-checking 1 000 random variant
  `chrom:pos:ref:alt` tuples against the FASTA.

### Step 1 — Convert and normalize

```bash
plink2 --bfile cohort_hg38 \
       --export vcf-4.2 bgz id-paste=iid \
       --output-chr chrM \
       --out cohort.raw

bcftools norm \
       --fasta-ref GRCh38_full_analysis_set_plus_decoy_hla.fa \
       --check-ref ws \
       --multiallelics -any \
       --do-not-normalize \
       -Oz -o cohort.norm.vcf.gz cohort.raw.vcf.gz
bcftools index -t cohort.norm.vcf.gz
```

Rationale: PLINK2 emits a VCF with consistent `chr` contigs; `bcftools norm
-c ws` *warns and sets* `REF` to the FASTA-supported allele where the
genotype data carry the alternate allele in the REF slot, which is the most
common cause of TOPMed-server rejections [Ref 3, Ref 9].

### Step 2 — Variant QC

| Filter | Threshold | Justification |
| --- | --- | --- |
| Non-ACGT / symbolic alleles | drop | Server requires `[ACGT]+` [Ref 3]. |
| Indels on a SNP array | drop | Array indels are rarely well-calibrated; TOPMed will impute indels from SNPs. |
| Duplicate `chrom:pos:ref:alt` | keep best call-rate | Server rejects duplicates [Ref 3]. |
| Variant missingness | `> 0.05` drop | Standard, e.g., Anderson et al. 2010 [Ref 4]. |
| MAF | `< 0.0001` drop | Below server's minimum useful frequency; lower is dominated by genotyping error [Ref 3]. |
| HWE | `p < 1e-6` within inferred ancestry | Avoid penalizing population-differentiated loci [Ref 4, Ref 5]. |
| Palindromic A/T or C/G with `0.40 < MAF < 0.60` | flag, exclude from strand inference | Strand ambiguous; cannot be resolved by AF alone [Ref 6]. |

PLINK2 implementation:

```bash
plink2 --vcf cohort.norm.vcf.gz \
       --geno 0.05 \
       --maf 0.0001 \
       --hwe 1e-6 midp \
       --snps-only just-acgt \
       --rm-dup exclude-all \
       --make-pgen --out cohort.varqc
```

HWE is re-run **per inferred ancestry group** after Step 4 and the union of
failed variants is excluded.

### Step 3 — Sample QC

| Filter | Threshold | Justification |
| --- | --- | --- |
| Sample call rate | `< 0.97` drop | Anderson et al. 2010 [Ref 4]. |
| Autosomal heterozygosity `F` | `|F| > 0.20` flag | Contamination / inbreeding outliers [Ref 4]. |
| Reported vs. genotypic sex | mismatch flag | chrX heterozygosity / chrY call count [Ref 4]. |

```bash
plink2 --pfile cohort.varqc --missing --het --out cohort.sampleqc
plink2 --pfile cohort.varqc --check-sex --out cohort.sex
```

### Step 4 — Relatedness and ancestry

- **Relatedness:** `KING --kinship` on LD-pruned common variants. Default
  policy: flag pairs with `kinship > 0.0884` (2nd-degree); for unrelated
  analyses, drop one of each pair preferring the higher call-rate sample
  [Ref 7].
- **Ancestry:** project the cohort onto a precomputed 1000 Genomes + HGDP
  PCA reference using PLINK2 `--score` with PC loadings; assign ancestry by
  random-forest classifier (gnomAD approach) [Ref 8]. Ancestry calls feed
  back into Step 2 HWE stratification.

### Step 5 — TOPMed reference harmonization

For each variant remaining after Step 2, compare `(chrom, pos, ref, alt)`
to the TOPMed bravo sites VCF
([`bravo-dbsnp-all.vcf.gz`](https://bravo.sph.umich.edu/freeze8/hg38/)).
Resolutions, in order:

1. Exact match → keep.
2. REF/ALT swapped → swap genotypes (`bcftools +fixref -m swap`).
3. Strand-flippable to a match for non-palindromic SNPs → flip
   (`bcftools +fixref -m flip`).
4. Palindromic with concordant AF (`|MAF_cohort − MAF_TOPMed| < 0.10` and
   both < 0.40) → keep with flip if needed; otherwise drop.
5. Otherwise → drop and log.

This is the same logic as Will Rayner's
`HRC-1000G-check-bim.pl`, adapted for the TOPMed sites file [Ref 9].

### Step 6 — Per-chromosome VCFs

```bash
for chr in $(seq 1 22) X; do
    bcftools view -r chr${chr} cohort.harmonized.vcf.gz \
        | bcftools sort -Oz -o vcf/chr${chr}.vcf.gz
    bcftools index -t vcf/chr${chr}.vcf.gz
done
```

Assertions before declaring success:

- Every contig begins with `chr`.
- Every record is sorted by position.
- `bcftools view -h` shows the GRCh38 contigs only.
- `bcftools stats` reports zero non-ACGT alleles and zero duplicates.

### Step 7 — Pre-submission validation

```bash
python checkVCF.py -r GRCh38_full_analysis_set_plus_decoy_hla.fa \
                   -o checkvcf vcf/chr*.vcf.gz
imputationbot validate vcf/chr*.vcf.gz       # optional, requires API token
```

Any non-zero exit aborts the workflow with a human-readable summary.

### Step 8 — Submission parameters

When uploading to the TOPMed Imputation Server:

| Parameter | Value | Reason |
| --- | --- | --- |
| Reference panel | **TOPMed r3** | Largest, most diverse panel for hg38 [Ref 1]. |
| Array build | **hg38** | Matches our prep. |
| Phasing | **Eagle v2.4** | Server default; recommended for unrelated samples [Ref 2]. |
| Population | **vs. TOPMed Panel** (mixed) | Avoids hard ancestry assignment when cohort is multi-ancestry. |
| rsq Filter | **off** | Filter `R2` post-hoc to keep flexibility for different downstream uses [Ref 11]. |
| AES-256 encryption | **on** | Required for many IRB protocols. |

### Step 9 — Post-imputation QC (downstream)

Documented here for completeness; implemented in a follow-up workflow:

- Distribution of `INFO`/`R2` per chromosome and per MAF bin.
- MAF concordance with gnomAD v4 (matched ancestry).
- Genotype concordance on a held-out 5 % of typed variants masked before
  upload (the "leave-one-out" sanity check) [Ref 11].

---

## 5. Decision points and defaults

| Decision | Default | When to override |
| --- | --- | --- |
| HWE threshold | `1e-6` per ancestry | Use `1e-10` for very large cohorts (>50 k) to keep FDR low. |
| Relatedness cut | drop one of each ≥2nd-degree pair | Family-based analyses: keep all, record pedigree. |
| Palindromic AF window | `0.40–0.60` flagged | Tighten to `0.45–0.55` for very small (<200) cohorts where AF estimates are noisy. |
| Ancestry reference | 1000G + HGDP (gnomAD release) | Use TOPMed BRAVO summary AF when working on TOPMed-internal cohorts. |
| Chromosome scope | `chr1–22, chrX` | chrY/chrM are not imputed by TOPMed; skip. |

Each default is encoded in `config/defaults.yaml` (to be added with the
implementation PR) so they can be overridden per run.

---

## 6. Quality metrics reported

The HTML report (`report.html`) will include, with figures:

- **Pre/post-QC counts**: variants and samples removed at each step.
- **Missingness** distributions (sample- and variant-level).
- **Heterozygosity vs call-rate** scatter for sample outlier review.
- **Sex check** plot (chrX F vs. reported sex).
- **PCA** of cohort projected onto 1000G+HGDP with ancestry assignments.
- **KING kinship** histogram and related-pair table.
- **AF concordance** (cohort vs. TOPMed BRAVO) scatter, with palindromic
  SNPs colored.
- **TableOne** demographic summary by inferred ancestry.

---

## 7. Reproducibility and CI

- Pinned tool versions in the Docker image; image SHA recorded in
  `provenance.json`.
- GitHub Actions:
  - Lint shell + R + Python (`shellcheck`, `lintr`, `ruff`).
  - Run the pipeline end-to-end on a tiny synthetic GRCh38 dataset
    (~50 samples × ~1 000 markers) generated on the fly.
  - Verify the `report.html` is produced and the per-chromosome VCFs pass
    `checkVCF.py`.

---

## 8. Failure modes intentionally caught before server submission

These are the failure modes most commonly seen on the TOPMed help forum
[Ref 3] and the precise check in this pipeline that prevents each:

| Failure mode | Pipeline check |
| --- | --- |
| Wrong build (GRCh37 submitted as hg38) | Step 0 build detection. |
| Non-ACGT alleles (`I`/`D`/`0`) | Step 2 `--snps-only just-acgt`. |
| `chr` prefix mismatch | Step 1 PLINK2 `--output-chr chrM`; Step 6 assertion. |
| Duplicate variant IDs/positions | Step 2 `--rm-dup`. |
| Ref-allele mismatch above threshold | Step 1 `bcftools norm -c ws` + Step 5 harmonization. |
| Unsorted VCF | Step 6 `bcftools sort`. |
| Monomorphic / MAF=0 variants | Step 2 `--maf 0.0001`. |
| Palindromic strand errors | Step 2 flag + Step 5 AF-based resolution. |

---

## 9. Multi-ancestry considerations

- HWE filtering, AF concordance, and palindromic AF resolution are all
  performed **within ancestry strata** inferred at Step 4 [Ref 5].
- We **do not** down-sample to a single ancestry. TOPMed r3 contains ~50 %
  non-European haplotypes [Ref 1], which is the principal reason it is
  preferred over HRC for diverse cohorts.
- The submission `Population` is set to `vs. TOPMed Panel` (mixed) so
  imputation is performed against the whole reference rather than a
  super-population subset.

---

## 10. GRCh37 → GRCh38 lift-over (for users with legacy data)

When input is on GRCh37:

```bash
# 1. PLINK -> VCF on b37 with chr-prefixed contigs
plink2 --bfile cohort_b37 --recode vcf-4.2 bgz --output-chr chrM \
       --out cohort.b37
# 2. CrossMap / Picard LiftoverVcf with hg19ToHg38 chain
CrossMap.py vcf hg19ToHg38.over.chain.gz cohort.b37.vcf.gz \
            GRCh38_full_analysis_set_plus_decoy_hla.fa cohort.hg38.vcf
# 3. Re-normalize against GRCh38 FASTA and re-enter the pipeline at Step 1
bcftools norm -f GRCh38_full_analysis_set_plus_decoy_hla.fa -c ws \
              -Oz -o cohort.hg38.norm.vcf.gz cohort.hg38.vcf
```

Variants that fail to lift over (typically <1 %) are logged in
`logs/liftover_failed.tsv`.

---

## 11. Open questions / future work

- Should we offer a "fast" mode that skips Step 4 ancestry projection when
  metadata already contains validated ancestry calls? Likely yes; gated on
  a `--trust-metadata-ancestry` flag.
- Add support for the **1000G 30x** reference panel as a fallback for
  cohorts not approved to use TOPMed.
- Integrate **GLIMPSE** for low-coverage WGS inputs (currently out of
  scope).

See [`references.md`](references.md) for the full bibliography.
