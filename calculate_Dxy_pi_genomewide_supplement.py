#!/usr/bin/env python3
"""Calculate genome-wide pi and Dxy from autosomal all-sites BCFs.

Reproduces the site-level filters and estimators used in the Fukushima
stickleback analysis (calculate_Dxy_windows.py), without window calculations.

Usage:
  python calculate_Dxy_pi_genomewide_supplement.py \
      --root /path/to/project \
      --bcf-dir /path/to/allsites \
      --output-dir /path/to/output

Requirements: Python 3.8+, bcftools available on PATH. Standard library only.
The input BCFs must contain GT and FORMAT/DP for all 51 samples.
"""

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path

POPS = ("AgaRiv", "Hakusan", "Inawashiro", "Kashira", "Kitakata", "Tajima")

# Map of unphased/phased diploid genotypes to ALT dosage (0, 1, or 2).
# Any other genotype string (e.g. "./.", "./0") resolves to None via .get().
_DOSAGE = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}


def read_table(path: Path, key_field: str, value_field: str) -> dict:
    """Return a dict mapping key_field → value_field from a TSV file."""
    with path.open(newline="") as f:
        return {r[key_field]: r[value_field] for r in csv.DictReader(f, delimiter="\t")}


def _parse_dosages(samples, gt_fields, dp_fields, max_dp):
    """Convert per-sample GT/DP fields to ALT dosages (int or None).

    Returns a list of length len(samples).  None means the sample is excluded
    at this site (missing DP, depth out of range, or unrecognised genotype).
    """
    dosages = []
    for i in range(len(samples)):
        raw_dp = dp_fields[i]
        if raw_dp in ("", "."):
            dosages.append(None)
            continue
        depth = int(raw_dp)
        if depth < 5 or depth > max_dp[i]:
            dosages.append(None)
            continue
        # Normalise phased genotypes before lookup
        genotype = gt_fields[i].replace("|", "/")
        dosages.append(_DOSAGE.get(genotype))  # None for "./." etc.
    return dosages


def calculate(bcf: Path, meta: dict, depth_means: dict) -> tuple:
    """Stream one BCF through bcftools query and accumulate pi / Dxy sums.

    Returns
    -------
    records : int
        Total lines emitted by bcftools (all sites).
    retained : int
        Sites that passed all filters and contributed to the statistics.
    skipped_counts : dict
        Breakdown of filtered-out sites by reason.
    pi_sum : defaultdict(float)
        Per-population sum of site-level pi (unscaled).
    dxy_sum : defaultdict(float)
        Per-pair sum of site-level Dxy (unscaled).
    """
    # ------------------------------------------------------------------ setup
    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(bcf)], text=True
    ).splitlines()

    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate sample IDs in {bcf}")
    if not set(samples) <= meta.keys():
        raise ValueError(f"Samples missing from metadata: {bcf}")
    if not set(samples) <= depth_means.keys():
        raise ValueError(f"Samples missing from depth table: {bcf}")

    groups = {
        pop: [i for i, s in enumerate(samples) if meta[s] == pop]
        for pop in POPS
    }
    missing_pops = [pop for pop, idx in groups.items() if not idx]
    if missing_pops:
        raise ValueError(f"Populations absent from {bcf}: {missing_pops}")

    max_dp = [3.0 * float(depth_means[s]) for s in samples]

    # --------------------------------------------------------- accumulators
    pi_sum  = defaultdict(float)
    dxy_sum = defaultdict(float)
    retained = records = 0
    skipped_counts = {
        "non_biallelic_ref": 0,   # REF is not a single ACGT base
        "non_biallelic_alt": 0,   # ALT is not . or a single ACGT SNP
        "low_callability":   0,   # ≥1 population below 80 % callable threshold
    }

    # --------------------------------------------------------- bcftools stream
    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP]\n",
        str(bcf),
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    try:
        for line in process.stdout:
            records += 1
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[:4]
            values = fields[4:]

            # --- site-level filters ---
            # Keep only biallelic SNPs (REF = single ACGT base)
            if len(ref) != 1 or ref not in "ACGT":
                skipped_counts["non_biallelic_ref"] += 1
                continue

            # ALT must be missing (monomorphic) or a single ACGT SNP.
            # Multi-allelic sites have commas in ALT and are rejected here.
            if alt not in (".", "A", "C", "G", "T"):
                skipped_counts["non_biallelic_alt"] += 1
                continue

            # Fail explicitly on malformed sample fields.
            if len(values) != 2 * len(samples):
                raise RuntimeError(
                    f"Sample field count mismatch in {bcf}: {chrom}:{pos}"
                )

            # --- per-sample genotype parsing ---
            gt_fields = values[0::2]
            dp_fields = values[1::2]
            dosages = _parse_dosages(samples, gt_fields, dp_fields, max_dp)

            # --- per-population frequency / diversity ---
            frequencies = {}
            diversities = {}
            low_callability = False

            for pop in POPS:
                called = [dosages[i] for i in groups[pop] if dosages[i] is not None]
                n = len(called)
                min_called = math.ceil(0.8 * len(groups[pop]))
                if n < min_called:
                    low_callability = True
                    break
                p = sum(called) / (2 * n)
                frequencies[pop] = p
                # Unbiased estimator of nucleotide diversity (Tajima 1983):
                # π̂ = 2p(1−p) · 2n/(2n−1)
                diversities[pop] = 2.0 * p * (1.0 - p) * (2 * n) / (2 * n - 1)

            if low_callability:
                skipped_counts["low_callability"] += 1
                continue

            # All populations passed → accumulate
            retained += 1
            for pop in POPS:
                pi_sum[pop] += diversities[pop]
            for pop_a, pop_b in combinations(POPS, 2):
                pa, pb = frequencies[pop_a], frequencies[pop_b]
                # Nei & Li (1979) site-level Dxy:
                # d̂xy = pa(1−pb) + pb(1−pa)
                dxy_sum[(pop_a, pop_b)] += pa * (1.0 - pb) + pb * (1.0 - pa)

    finally:
        process.stdout.close()

    if process.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {bcf}")

    return records, retained, skipped_counts, pi_sum, dxy_sum


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--root",       type=Path, default=Path("."))
    parser.add_argument("--bcf-dir",    type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--pattern",    default="*.allsites.bcf", help="BCF filename glob")
    parser.add_argument("--expected-chromosomes", type=int, default=20)
    args = parser.parse_args()

    root       = args.root
    bcf_dir    = args.bcf_dir    or root / "10_Dxy/allsites"
    output_dir = args.output_dir or root / "10_Dxy/reproducibility_check"

    bcfs = sorted(bcf_dir.glob(args.pattern))
    if len(bcfs) != args.expected_chromosomes:
        parser.error(
            f"Expected {args.expected_chromosomes} BCFs, "
            f"found {len(bcfs)} in {bcf_dir}"
        )

    meta         = read_table(root / "00_metadata/sample_sheet_analysis.tsv", "sample_id", "population")
    depth_means  = read_table(root / "00_metadata/coverage_QC_51.tsv",        "sample_id", "mean_depth")

    output_dir.mkdir(parents=True, exist_ok=True)

    pi_totals  = defaultdict(float)
    dxy_totals = defaultdict(float)
    total_sites = 0

    # ---------------------------------------------------------- per-chromosome log
    callability_path = output_dir / "chromosome_callability.tsv"
    with callability_path.open("w") as f:
        f.write(
            "BCF\tRecords\tCallable_sites"
            "\tSkipped_non_biallelic_ref\tSkipped_non_biallelic_alt"
            "\tSkipped_low_callability\n"
        )
        for bcf in bcfs:
            records, sites, skipped, pi, dxy = calculate(bcf, meta, depth_means)
            f.write(
                f"{bcf.name}\t{records}\t{sites}"
                f"\t{skipped['non_biallelic_ref']}"
                f"\t{skipped['non_biallelic_alt']}"
                f"\t{skipped['low_callability']}\n"
            )
            total_sites += sites
            for k, v in pi.items():
                pi_totals[k] += v
            for k, v in dxy.items():
                dxy_totals[k] += v
            print(f"{bcf.name}: {sites:,} callable sites", flush=True)

    if not total_sites:
        raise ValueError("No callable sites found across all BCFs.")

    # ----------------------------------------------------------------- pi output
    pi_path = output_dir / "genomewide_pi_common.tsv"
    with pi_path.open("w") as f:
        f.write("Population\tCallable_sites\tPi\n")
        for pop in sorted(POPS):
            pi_val = pi_totals[pop] / total_sites
            f.write(f"{pop}\t{total_sites}\t{pi_val:.10f}\n")

    # --------------------------------------------------------------- Dxy output
    dxy_path = output_dir / "genomewide_Dxy_common.tsv"
    with dxy_path.open("w") as f:
        f.write("Population1\tPopulation2\tCallable_sites\tDxy\n")
        for pop_a, pop_b in combinations(POPS, 2):
            dxy_val = dxy_totals[(pop_a, pop_b)] / total_sites
            f.write(f"{pop_a}\t{pop_b}\t{total_sites}\t{dxy_val:.10f}\n")

    print(f"\nTotal callable sites : {total_sites:,}")
    print(f"Outputs written to   : {output_dir}")


if __name__ == "__main__":
    main()
