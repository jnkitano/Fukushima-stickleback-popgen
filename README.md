# Population genomic analyses of inland threespine sticklebacks in Fukushima, Japan

This repository provides a Python script for calculating genome-wide nucleotide diversity (π) within populations and absolute nucleotide divergence (Dxy) between populations from whole-genome sequencing data.

The script accompanies a study of six inland populations of *Gasterosteus aculeatus* in Fukushima Prefecture, Japan.

## Script

`calculate_Dxy_pi_genomewide_supplement.py` calculates genome-wide π and Dxy from all-sites BCF files covering 20 autosomes. It also reports site filtering statistics for each chromosome.

## Requirements

- Python 3.8 or later (standard library only)
- `bcftools` available on `PATH`

## Input files

The project directory should contain:

- `00_metadata/sample_sheet_analysis.tsv` with columns `sample_id` and `population`.
- `00_metadata/coverage_QC_51.tsv` with columns `sample_id` and `mean_depth`.
- `10_Dxy/allsites/` with 20 files matching `*.allsites.bcf`.

Each BCF must include genotype (`GT`) and per-sample depth (`FORMAT/DP`) information. Sample IDs must match those in both metadata tables. The population labels used in the script are `AgaRiv`, `Hakusan`, `Inawashiro`, `Kashira`, `Kitakata`, and `Tajima`.

## Filtering and calculations

Both invariant sites and biallelic SNPs are included; indels and multiallelic sites are excluded. A genotype is callable when its depth is at least 5 and no greater than three times the mean sequencing depth of that individual. Sites are retained only when at least 80% of individuals in **each** of the six populations have callable genotypes.

At each retained site, nucleotide diversity is calculated as π = 2p(1 − p) × 2n/(2n − 1), where p is the alternative allele frequency and n is the number of diploid individuals with callable genotypes. For populations A and B, Dxy = pA(1 − pB) + pB(1 − pA). Genome-wide values are averages over the same retained sites for all populations and pairs.

## Usage

```bash
python3 calculate_Dxy_pi_genomewide_supplement.py \
  --root /path/to/project \
  --bcf-dir /path/to/project/10_Dxy/allsites \
  --output-dir /path/to/project/10_Dxy/reproducibility_check
```

`--root` is optional and defaults to the current working directory (`.`). `--bcf-dir` and `--output-dir` can also be omitted when using the directory structure above.

## Outputs

- `genomewide_pi_common.tsv`: π for each population.
- `genomewide_Dxy_common.tsv`: Dxy for each population pair.
- `chromosome_callability.tsv`: per-BCF site counts, with columns `BCF`, `Records`, `Callable_sites`, `Skipped_non_biallelic_ref`, `Skipped_non_biallelic_alt`, and `Skipped_low_callability`.

## Reproducibility

Using the study's all-sites BCF files, the script reproduced the manuscript's π and Dxy estimates to six decimal places.

Input BCF files are not included in this repository. The script reproduces π and Dxy calculations from the specified inputs, not the upstream read mapping or variant calling.

## Associated manuscript

*Population genetic structure and diversity of inland threespine sticklebacks in Fukushima, Japan* (publication details to be added).

## References

Tajima F (1983) Evolutionary relationship of DNA sequences in finite populations. *Genetics* 105:437–460. https://doi.org/10.1093/genetics/105.2.437

Nei M, Li WH (1979) Mathematical model for studying genetic variation in terms of restriction endonucleases. *Proc Natl Acad Sci USA* 76:5269–5273. https://doi.org/10.1073/pnas.76.10.5269
