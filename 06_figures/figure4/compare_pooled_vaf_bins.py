"""Four prespecified one-sided N/S comparisons after merging Figure 4 lineages.

Reads the validated independent coding-site table. No new site or cell filter,
no alteration to original mean VAF, annotation, bins, PAML or Figure 4.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
CODE = Path(os.environ["MITOEVO_SOURCE_CODE_ROOT"])
SOURCE = CODE / "mitoevo_targeted_analysis/Tables/coding_site_all_unique.tsv"
BG = CODE / "manuscript_revision_20261104_v1/remote_outputs/tables/background_sites_summary.csv"
PRIOR = CODE / "manuscript_revision_20261104_v1/tables/fig4_dnds_CI.tsv"
LINEAGE_SITES = CODE / "manuscript_revision_20261104_v1/remote_outputs/tables/fig4_site_vaf_both_filters.tsv"
FIGURE = CODE / "manuscript_story_task1/Figures/Figure4.pdf"
SP = ["human","mouse","rat","chicken","cobra","sturgeon","stingray"]
BINS = ["(0,0.1]","(0.1,0.3]","(0.3,0.5]","(0.5,1]"]
PAIRS = [(2,0),(2,1),(3,0),(3,1)]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh_planned(pvals: pd.Series, family_n: int) -> pd.Series:
    """BH with planned family size; skipped tests stay NA, count as p=1."""
    result = pd.Series(np.nan,index=pvals.index,dtype=float)
    valid = pvals.dropna().sort_values()
    if valid.empty:
        return result
    q = np.minimum.accumulate((valid.to_numpy()*family_n/np.arange(1,len(valid)+1))[::-1])[::-1]
    result.loc[valid.index] = np.minimum(q,1)
    return result


def main() -> None:
    original_figure = sha(FIGURE)
    sites = pd.read_csv(SOURCE,sep="\t")
    bg = pd.read_csv(BG).set_index("species")
    assert len(sites)==977 and set(sites.species)==set(SP)
    assert not sites.duplicated(["species","mutation"]).any()
    assert set(sites["class"])=={"nonsynonymous","synonymous"}
    assert sites.mean_vaf.between(0,1,inclusive="right").all()
    assert (sites.n_observations>0).all()
    expected_bin = pd.cut(sites.mean_vaf,[0,.1,.3,.5,1],labels=BINS).astype(str)
    assert (expected_bin==sites.bin).all()
    original_lineage=pd.read_csv(LINEAGE_SITES,sep="\t")
    original_lineage=original_lineage[original_lineage.filter_mode.eq("legacy_all_positive")]
    original_lineage=original_lineage.merge(sites[["species","mutation"]],on=["species","mutation"],validate="many_to_one")
    assert not original_lineage.duplicated(["species","mutation","lineage"]).any()
    recalc=(original_lineage.assign(vaf_sum=lambda x:x.mean_vaf*x.n_observations)
            .groupby(["species","mutation"],as_index=False)
            .agg(vaf_sum=("vaf_sum","sum"),n_observations=("n_observations","sum")))
    recalc["mean_vaf_check"]=recalc.vaf_sum/recalc.n_observations
    check=sites.merge(recalc,on=["species","mutation"],validate="one_to_one",suffixes=("","_recalc"))
    assert len(check)==len(sites)
    assert (check.n_observations==check.n_observations_recalc).all()
    assert np.allclose(check.mean_vaf,check.mean_vaf_check,rtol=1e-12,atol=1e-12)

    rows=[]
    for sp in SP:
        z=sites[sites.species.eq(sp)]
        bg_ratio=float(bg.loc[sp,"non_synonymous_sites"]/bg.loc[sp,"synonymous_sites"])
        for i,label in enumerate(BINS,1):
            g=z[z.bin.eq(label)]
            N=int(g.is_nonsyn.sum());S=len(g)-N
            rows.append(dict(species=sp,bin_number=i,bin=label,N_obs=N,S_obs=S,
                             n_unique_coding_sites=N+S,PAML_N_over_S_background=bg_ratio,
                             dnds_original_background=N/S/bg_ratio if S>0 else np.nan,
                             dnds_status="defined" if S>0 else ("empty" if N==0 else "undefined_S_zero"),
                             n_positive_cell_variant_observations=int(g.n_observations.sum())))
    bins=pd.DataFrame(rows)
    assert len(bins)==28 and bins.n_unique_coding_sites.sum()==len(sites)

    # Confirm all 28 pooled rows against the previously validated All_unique table.
    prior=pd.read_csv(PRIOR,sep="\t")
    prior=prior[(prior.filter_mode.eq("legacy_all_positive")) &
                (prior.lineage.eq("All_unique")) & prior.bin.isin(BINS)]
    assert len(prior)==28
    merged=bins.merge(prior,on=["species","bin"],validate="one_to_one")
    assert (merged.N_obs==merged.N).all() and (merged.S_obs==merged.S).all()
    # Prior helper uses +Inf for N>0,S=0; this delivery explicitly records NA instead.
    for row in merged.itertuples():
        if row.S_obs>0:
            assert np.isclose(row.dnds_original_background,row.dnds,rtol=1e-12,atol=1e-12)

    pairs=[]
    for sp in SP:
        block=bins[bins.species.eq(sp)].set_index("bin_number")
        for high_i,low_i in PAIRS:
            high=block.loc[high_i+1];low=block.loc[low_i+1]
            high_N,high_S=int(high.N_obs),int(high.S_obs)
            low_N,low_S=int(low.N_obs),int(low.S_obs)
            nh, nl=high_N+high_S,low_N+low_S
            rec=dict(species=sp,comparison=f"Bin {high_i+1} vs Bin {low_i+1}",
                     high_bin_number=high_i+1,high_bin=high.bin,
                     low_bin_number=low_i+1,low_bin=low.bin,
                     high_N_obs=high_N,high_S_obs=high_S,high_n_sites=nh,
                     high_dnds=high.dnds_original_background,
                     high_dnds_status=high.dnds_status,
                     low_N_obs=low_N,low_S_obs=low_S,low_n_sites=nl,
                     low_dnds=low.dnds_original_background,
                     low_dnds_status=low.dnds_status,
                     PAML_N_over_S_background=high.PAML_N_over_S_background,
                     odds_ratio_high_vs_low=np.nan,p_one_sided_less=np.nan,
                     status="empty_bin_skipped" if nh==0 or nl==0 else "tested",
                     high_bin_under_5_sites=nh<5)
            if nh and nl:
                result=fisher_exact([[high_N,high_S],[low_N,low_S]],alternative="less")
                odds=float(result.statistic)
                rec["odds_ratio_high_vs_low"]=odds
                rec["p_one_sided_less"]=float(result.pvalue)
                if np.isnan(odds):
                    rec["effect_direction"]="undefined_no_N_S_contrast"
                    rec["status"]="tested_no_N_S_contrast"
                elif odds<1:
                    rec["effect_direction"]="high_bin_lower_N_over_S"
                elif odds>1:
                    rec["effect_direction"]="high_bin_higher_N_over_S"
                else:
                    rec["effect_direction"]="equal_sample_N_over_S"
            else:
                rec["effect_direction"]="not_estimable_empty_bin"
            pairs.append(rec)
    tests=pd.DataFrame(pairs)
    assert len(tests)==28
    tests["BH_q_within_species_four_planned"] = np.nan
    for sp in SP:
        idx=tests.index[tests.species.eq(sp)]
        tests.loc[idx,"BH_q_within_species_four_planned"] = bh_planned(tests.loc[idx,"p_one_sided_less"],4)
    tests["BH_q_all_species_28_planned"] = bh_planned(tests.p_one_sided_less,28)
    tests["n_planned_tests_species_family"]=4
    tests["n_planned_tests_global_family"]=28
    tests["n_tested_tests_species_family"]=tests.groupby("species").p_one_sided_less.transform("count")
    tests.to_csv(HERE/"four_prespecified_one_sided_fisher.tsv",sep="\t",index=False,float_format="%.15g")
    bins.to_csv(HERE/"pooled_lineages_four_vaf_bin_counts.tsv",sep="\t",index=False,float_format="%.15g")
    assert sha(FIGURE)==original_figure
    checks=dict(source=str(SOURCE),source_sha256=sha(SOURCE),
                original_Figure4_sha256_unchanged=original_figure,
                pooled_unique_coding_sites=len(sites),bin_rows=len(bins),
                original_lineage_site_rows=len(original_lineage),
                shared_across_lineages=len(original_lineage)-len(sites),
                pooled_mean_vaf_recomputed_from_positive_observation_weights=True,
                planned_pairs=len(tests),tested_pairs=int(tests.p_one_sided_less.notna().sum()),
                skipped_empty_pairs=int(tests.p_one_sided_less.isna().sum()),
                within_species_BH_family="4 prespecified comparisons per species; empty pairs count conservatively as p=1 for family size, output P/q remains NA",
                global_BH_family="28 prespecified comparisons across 7 species; empty pairs count conservatively as p=1 for family size, output P/q remains NA",
                test="one-sided Fisher exact, rows high/low VAF and columns N/S, alternative high N/S < low N/S",
                dnds_S_zero="NA (undefined), never set to zero")
    (HERE/"validation.json").write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding="utf-8")
    print(bins.to_string(index=False))
    print(tests[["species","comparison","high_N_obs","high_S_obs","low_N_obs","low_S_obs",
                 "odds_ratio_high_vs_low","p_one_sided_less","BH_q_within_species_four_planned",
                 "BH_q_all_species_28_planned","status"]].to_string(index=False))
    print(json.dumps(checks,ensure_ascii=False,indent=2))


if __name__=="__main__":main()
