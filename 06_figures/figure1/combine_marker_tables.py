#!/usr/bin/env python3
import csv
import pathlib
import re
import subprocess
from collections import defaultdict
ROOT=pathlib.Path('/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC')
OUT=ROOT/'annotation_visual_revision/v3'
CODE=ROOT/'code/annotation_visual_revision/v3'
SPECIES=('human','mouse','rat','chicken','cobra','sturgeon','stingray')
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t')
        w.writeheader();w.writerows(rows)
for suffix,target in [('marker_evidence','seven_species_marker_evidence.tsv'),
                      ('candidate_availability','seven_species_candidate_availability.tsv'),
                      ('subtype_coverage','seven_species_subtype_coverage.tsv'),
                      ('umap_number_key','seven_species_umap_number_key.tsv')]:
    rows=[r for sp in SPECIES for r in read(OUT/f'{sp}.{suffix}.tsv')]
    write(OUT/target,rows)
markers=read(OUT/'seven_species_marker_evidence.tsv')
avail=read(OUT/'seven_species_candidate_availability.tsv')
coverage=read(OUT/'seven_species_subtype_coverage.tsv')
bygroup=defaultdict(list)
for r in markers:bygroup[(r['species'],r['celltype'])].append(r)
rows=[]
for c in coverage:
    sp,ct=c['species'],c['celltype']
    shown=[r for r in bygroup[(sp,ct)] if r['display']=='TRUE']
    missing=[r['gene'] for r in avail if r['species']==sp and r['celltype']==ct and r['present_in_matrix']=='FALSE']
    if not shown:level='no_candidate_displayed'
    elif all(r['assessment']=='unverified_identity' for r in shown):level='unverified_ID_candidate_only'
    elif any(r['assessment']=='uncertain_subtype' for r in shown):level='lineage_only_subtype_unresolved'
    elif any(r['signal_status']!='relative_enriched' for r in shown):level='shared_or_weak_signal'
    else:level='biological_candidates_with_relative_signal'
    rows.append(dict(species=sp,sub_celltype=ct,n_cells=c['n_cells'],
                     n_displayed=len(shown),displayed_features=';'.join(r['plot_label'] for r in shown),
                     candidate_features_absent=';'.join(missing),evidence_level=level,
                     interpretation=('Unassigned category; no subtype markers required.' if ct=='Others' else
                     'Two candidates displayed; inspect assessment and signal_status for lineage specificity.' if len(shown)>=2 else
                     'Fewer than two usable candidates found; no unsupported replacement was invented.')))
write(OUT/'seven_species_subtype_evidence_summary.tsv',rows)
print('combined',len(markers),'marker rows and',len(rows),'subtype rows')
