"""Pairwise two-sided Fisher tests of original Figure 4 VAF-bin N/S counts.

This script reads the validated, unchanged 56-row all_positive dN/dS table.
It does not recompute sites, VAF, PAML backgrounds, dN/dS, or the figure.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


HERE = Path(__file__).resolve().parent
TASK1 = Path(os.environ["MITOEVO_FIG4_TASK_ROOT"])
OUT = TASK1 / "Figure4_pairwise_20261105_v1"
SOURCE = TASK1 / "Tables" / "figure4_56_original_results_with_display_labels.tsv"
ORIGINAL = Path(os.environ["MITOEVO_SOURCE_CODE_ROOT"]) / "manuscript_revision_20261104_v1" / "remote_outputs" / "tables" / "old_AllSpecies_dot_dnds_results_withP.tsv"
FIGURE = TASK1 / "Figures" / "Figure4.pdf"
BINS = ["(0,0.1]", "(0.1,0.3]", "(0.3,0.5]", "(0.5,1]"]
SPECIES = ["human", "mouse", "rat", "chicken", "cobra", "sturgeon", "stingray"]
LINEAGES = ["T_B_NK", "Mono"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(p_values: list[float]) -> list[float]:
    """BH adjusted values in input order, across one species × lineage family."""
    if not p_values:
        return []
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranks = np.arange(1, len(p) + 1)
    adjusted_sorted = np.minimum.accumulate((p[order] * len(p) / ranks)[::-1])[::-1]
    adjusted = np.empty_like(p)
    adjusted[order] = np.minimum(adjusted_sorted, 1)
    return adjusted.tolist()


def validate_input(df: pd.DataFrame) -> None:
    assert len(df) == 56
    assert set(df.species) == set(SPECIES)
    assert set(df.celltype_group) == set(LINEAGES)
    assert set(df.bin) == set(BINS)
    assert not df.duplicated(["species", "celltype_group", "bin"]).any()
    assert (df.nonsynonymous >= 0).all() and (df.synonymous >= 0).all()
    assert (df.nonsynonymous + df.synonymous == df.total_mut).all()
    assert (df.groupby(["species", "celltype_group"]).size() == 4).all()
    for row in df.itertuples():
        if row.synonymous > 0:
            expected = row.nonsynonymous / row.synonymous / row.NS_ratio_expected
            assert np.isclose(expected, row.dnds, rtol=1e-12, atol=1e-12)
        else:
            assert pd.isna(row.dnds)

    # Independently archived original output must match the table used for this analysis.
    old = pd.read_csv(ORIGINAL, sep="\t")
    keys = ["species", "celltype_group", "bin"]
    merged = df.merge(old, on=keys, suffixes=("_copy", "_original"), validate="one_to_one")
    assert len(merged) == 56
    for col in ["nonsynonymous", "synonymous", "total_mut"]:
        assert (merged[f"{col}_copy"] == merged[f"{col}_original"]).all()
    for col in ["dnds", "p_binom", "NS_ratio_expected"]:
        assert np.allclose(
            merged[f"{col}_copy"], merged[f"{col}_original"],
            rtol=1e-12, atol=1e-12, equal_nan=True,
        )


def main() -> None:
    df = pd.read_csv(SOURCE, sep="\t")
    validate_input(df)
    before_figure = sha256(FIGURE)
    tests: list[dict] = []
    skipped: list[dict] = []
    summaries: list[dict] = []

    for species in SPECIES:
        for lineage in LINEAGES:
            block = df[(df.species == species) & (df.celltype_group == lineage)].set_index("bin")
            family: list[dict] = []
            family_skipped: list[dict] = []
            for bin1, bin2 in itertools.combinations(BINS, 2):
                a, b = block.loc[bin1], block.loc[bin2]
                n1, s1 = int(a.nonsynonymous), int(a.synonymous)
                n2, s2 = int(b.nonsynonymous), int(b.synonymous)
                t1, t2 = n1 + s1, n2 + s2
                base = {
                    "species": species, "lineage": lineage, "bin1": bin1, "bin2": bin2,
                    "bin1_N_obs": n1, "bin1_S_obs": s1, "bin1_coding_mutations": t1,
                    "bin1_dnds_original": a.dnds,
                    "bin2_N_obs": n2, "bin2_S_obs": s2, "bin2_coding_mutations": t2,
                    "bin2_dnds_original": b.dnds,
                }
                if t1 == 0 or t2 == 0:
                    family_skipped.append({**base, "reason": "one_or_both_bins_empty"})
                    continue
                # scipy returns the cross-product/sample odds ratio and a two-sided exact P.
                result = fisher_exact([[n1, s1], [n2, s2]], alternative="two-sided")
                odds = float(result.statistic)
                if np.isnan(odds):
                    direction = "undefined_no_N_vs_S_contrast"
                elif odds > 1:
                    direction = "bin1_higher_nonsynonymous_odds"
                elif odds < 1:
                    direction = "bin2_higher_nonsynonymous_odds"
                else:
                    direction = "equal_sample_odds"
                family.append({
                    **base,
                    "N_obs_both_bins": n1 + n2,
                    "S_obs_both_bins": s1 + s2,
                    "coding_mutations_both_bins": t1 + t2,
                    "odds_ratio_N_over_S_bin1_vs_bin2": odds,
                    "effect_direction": direction,
                    "fisher_p_two_sided": float(result.pvalue),
                    "low_count_bin_under_5": t1 < 5 or t2 < 5,
                    "zero_N_or_S_cell": min(n1, s1, n2, s2) == 0,
                })
            q_values = bh_adjust([item["fisher_p_two_sided"] for item in family])
            for item, q in zip(family, q_values, strict=True):
                item["bh_p_within_species_lineage"] = q
                item["n_tests_in_BH_family"] = len(family)
                item["bh_significant_0_05"] = q < 0.05
                tests.append(item)
            skipped.extend(family_skipped)
            significant = [x for x in family if x["bh_significant_0_05"]]
            raw_only = [x for x in family if x["fisher_p_two_sided"] < 0.05 and not x["bh_significant_0_05"]]
            sparse = [x for x in family if x["low_count_bin_under_5"]]
            label = lambda x: f'{x["bin1"]} vs {x["bin2"]}'
            summaries.append({
                "species": species, "lineage": lineage,
                "bin_counts_N_plus_S": "; ".join(f"{name}:{int(block.loc[name, 'total_mut'])}" for name in BINS),
                "nonempty_bins": int((block.total_mut > 0).sum()),
                "tested_pairs": len(family), "skipped_empty_pairs": len(family_skipped),
                "bh_significant_pairs_q_lt_0_05": "; ".join(label(x) for x in significant),
                "raw_p_lt_0_05_only_pairs": "; ".join(label(x) for x in raw_only),
                "pairs_with_bin_under_5_mutations": "; ".join(label(x) for x in sparse),
            })

    assert len(tests) + len(skipped) == 14 * 6
    assert sha256(FIGURE) == before_figure
    table_dir = OUT / "Tables"
    docs_dir = OUT / "Docs"
    table_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(tests).to_csv(table_dir / "figure4_vaf_bin_pairwise_fisher.tsv", sep="\t", index=False, float_format="%.15g")
    pd.DataFrame(skipped).to_csv(table_dir / "figure4_vaf_bin_pairs_skipped_empty.tsv", sep="\t", index=False, float_format="%.15g")
    pd.DataFrame(summaries).to_csv(table_dir / "figure4_species_lineage_pairwise_summary.tsv", sep="\t", index=False)
    validation = {
        "input_table": str(SOURCE), "input_sha256": sha256(SOURCE),
        "original_archived_table": str(ORIGINAL), "original_sha256": sha256(ORIGINAL),
        "original_figure": str(FIGURE), "original_figure_sha256_unchanged": before_figure,
        "n_original_rows": len(df), "n_possible_pairs": 84,
        "n_tested_pairs": len(tests), "n_skipped_empty_pairs": len(skipped),
        "test": "two-sided Fisher exact on [[N_bin1,S_bin1],[N_bin2,S_bin2]]",
        "odds_ratio": "sample cross-product: (N1*S2)/(S1*N2); bin1 versus bin2",
        "multiple_testing": "Benjamini-Hochberg within each species x lineage among tested nonempty pairs",
        "low_count_under_5": "descriptive sparsity flag only; does not exclude a comparison",
    }
    (docs_dir / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(validation, indent=2, ensure_ascii=False))
    print(pd.DataFrame(tests).loc[lambda x: x.bh_significant_0_05, [
        "species", "lineage", "bin1", "bin2", "bin1_N_obs", "bin1_S_obs",
        "bin2_N_obs", "bin2_S_obs", "odds_ratio_N_over_S_bin1_vs_bin2",
        "fisher_p_two_sided", "bh_p_within_species_lineage",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
