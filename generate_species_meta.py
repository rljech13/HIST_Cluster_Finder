#!/usr/bin/env python3
"""
generate_species_meta.py

Generates species_meta.tsv, using assembly_accession from hits with organism name from summary:
 1. Reads unique assembly IDs from hits (parse first GCA_/GCF_\d+\.\d+ in column 'subject')
 2. Downloads summary.tsv from pandas, automatically опpickingределяя:
    - accession column (contains 'accession')
    - spieces name column (contains 'organism' or 'species')
 3. Filters summary by accession list from hits
 4. Writes TSV: genomic_id<tab>species_label
"""
import argparse
import re
import sys

import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--hits',    required=True, help='TSV with H2B7 hiyts')
    p.add_argument('--summary', required=True, help='Assembly summary TSV')
    p.add_argument('--out',     required=True, help='Output species_meta.tsv')
    return p.parse_args()

def extract_assemblies(hits_path):
    df = pd.read_csv(hits_path, sep='\t', header=None, dtype=str, comment='#')
    subjects = df.iloc[:,1].dropna().astype(str)
    pattern = re.compile(r'^(GCA_\d+\.\d+|GCF_\d+\.\d+)')
    accs = subjects.map(lambda s: pattern.match(s).group(1) if pattern.match(s) else None)
    return sorted({a for a in accs if a})

def detect_columns(df):
    cols = df.columns.tolist()
    acc_cols = [c for c in cols if 'accession' in c.lower()]
    name_cols = [c for c in cols if any(x in c.lower() for x in ('organism','species'))]
    if not acc_cols:
        sys.exit('[ERROR] Can`t find column with accession in summary')
    if not name_cols:
        sys.exit('[ERROR] Can`t find column with organism/species in summary')
    return acc_cols[0], name_cols[0]

def main():
    args = parse_args()

    # 1) Getting list of accession from hits
    assemblies = extract_assemblies(args.hits)
    if not assemblies:
        sys.exit('[ERROR] Couldn`t manage to find any assembly from hits')

    # 2) Reading summary
    df = pd.read_csv(args.summary, sep='\t', dtype=str, comment='#', low_memory=False)
    acc_col, name_col = detect_columns(df)

    # 3) Filter by assemblies
    df_sel = df[df[acc_col].isin(assemblies)].copy()
    if df_sel.empty:
        sys.exit('[ERROR] Couldn`t find any assembly from hits in summary')
    df_sel = df_sel[[acc_col, name_col]].drop_duplicates()
    df_sel.columns = ['genomic_id','species_label']

    missing = set(assemblies) - set(df_sel['genomic_id'])
    if missing:
        print(f"[WARN] {len(missing)} assemblies haven`t foun summary:", file=sys.stderr)
        print(','.join(sorted(missing)), file=sys.stderr)

    # 4) Сохраняем TSV
    df_sel.to_csv(args.out, sep='\t', index=False)
    print(f"✔ Записано {len(df_sel)} строк в {args.out}")

if __name__=='__main__':
    main()