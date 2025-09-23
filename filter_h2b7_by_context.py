#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import argparse
import pandas as pd
from collections import defaultdict

###############################################################################
# I/O helpers
###############################################################################

def load_meta(path):
    """Load mapping genomic_id -> species_label from a meta table.

    Supported formats:
      - With or without header:  genomic_id<TAB>species_label
      - Or:                      genomic_id<2+ spaces>species_label
    """
    import re as _re
    meta = {}
    with open(path, 'r', encoding='utf-8') as f:
        first = True
        for line in f:
            line = line.rstrip('\n').strip()
            if not line or line.startswith('#'):
                continue
            if first and ('genomic_id' in line.lower() and 'species_label' in line.lower()):
                first = False
                continue
            first = False
            parts = _re.split(r'\t+|\s{2,}', line, maxsplit=1)
            if len(parts) < 2:
                continue
            gid = parts[0].strip()
            sp  = parts[1].strip()
            if gid and sp:
                meta[gid] = sp
    return meta


def parse_sseqid(sseqid):
    """Parse sseqid of form:
         <ASM>_genomic_<cluster>_<geneToken>_...
       Return (genome_id, cluster_id, gene_idx:int or None, gene_token:str)."""
    if '_genomic_' not in sseqid:
        return None, None, None, None
    asm_part, tail = sseqid.split('_genomic_', 1)
    asm = asm_part + '_genomic'
    parts = tail.split('_')
    if len(parts) < 2:
        return None, None, None, None

    cluster = parts[0]
    gene_token = parts[1]

    m = re.match(r'^(GCA_|GCF_)\d+\.\d+', asm)
    genome_id = m.group(0) if m else asm

    num = re.findall(r'\d+', gene_token)
    gene_idx = int(num[0]) if num else None

    return genome_id, cluster, gene_idx, gene_token


def read_fasta_dict(path):
    """Read FASTA into dict {header: seq}; header is the first token after '>'."""
    seqs = {}
    header = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(buf)
                header = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if header is not None:
            seqs[header] = ''.join(buf)
    return seqs


def write_fasta(records, out_path):
    """Write [(header, seq), ...] to FASTA."""
    with open(out_path, 'w') as out:
        for h, s in records:
            out.write(f">{h}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")


def sanitize_species(name):
    return (name.replace(' ', '_')
                .replace('\t', '_')
                .replace('(', '_')
                .replace(')', '_'))


###############################################################################
# Loading labels
###############################################################################

def load_labels_from_hits(hits_tsv):
    """Load labels from TSV.

    Supports:
      - 3+ cols: qseqid sseqid label => take [1,2]
      - 2 cols : sseqid label        => take [0,1]
    Returns DataFrame with columns: sseqid, label
    """
    df = pd.read_csv(hits_tsv, sep='\t', header=None, dtype=str, comment='#')
    if df.shape[1] >= 3:
        sdf = df[[1, 2]].copy()
    elif df.shape[1] >= 2:
        sdf = df[[0, 1]].copy()
    else:
        sys.exit(f"[ERROR] {hits_tsv}: need ≥2 columns")
    sdf.columns = ['sseqid', 'label']
    sdf['sseqid'] = sdf['sseqid'].astype(str).str.strip()
    sdf['label']  = sdf['label'].astype(str).str.strip()
    sdf = sdf.dropna(subset=['sseqid', 'label'])
    return sdf


###############################################################################
# Context selection
###############################################################################

def select_h2b7_by_context(labels_df):
    """Select TRUE H2B7 by neighbor context inside (genome_id, cluster_id) blocks.

    Rule A (priority):   H4 | H2B7 | H4
    Rule B (fallback):   neighbors are H1 and H2A-IV (any order),
                         applied ONLY if cluster has no Rule-A hits.

    Returns:
        dict with keys:
          - 'A': [sseqid, ...]
          - 'B': [sseqid, ...]
          - 'report': list of dict rows (for TSV)
          - 'all_h2b7': [sseqid, ...]
          - 'other': [sseqid, ...]
    """
    parsed = labels_df['sseqid'].apply(parse_sseqid)
    labels_df[['genome_id', 'cluster_id', 'gene_idx', 'gene_token']] = pd.DataFrame(parsed.tolist(),
                                                                                    index=labels_df.index)
    all_h2b7_ids = labels_df[labels_df['label'] == 'H2B7']['sseqid'].unique().tolist()

    labels_ctx = labels_df.dropna(subset=['genome_id', 'cluster_id', 'gene_idx']).copy()
    if labels_ctx.empty:
        return {'A': [], 'B': [], 'report': [], 'all_h2b7': all_h2b7_ids, 'other': all_h2b7_ids}

    labels_ctx['gene_idx'] = labels_ctx['gene_idx'].astype(int)

    groups = defaultdict(list)
    for row in labels_ctx.itertuples(index=False):
        groups[(row.genome_id, row.cluster_id)].append((row.gene_idx, row.sseqid, row.label))
    for k in groups:
        groups[k].sort(key=lambda x: x[0])

    keepA, keepB, rows_report = [], [], []

    for (genome_id, cluster_id), genes in groups.items():
        local_A, local_B = [], []

        for i, (gidx, sid, lab) in enumerate(genes):
            if lab != 'H2B7':
                continue
            if i == 0 or i == len(genes) - 1:
                rows_report.append({
                    'genome_id': genome_id, 'cluster_id': cluster_id,
                    'gene_idx': gidx, 'sseqid': sid,
                    'rule': 'NA_no_neighbors',
                    'left_label': genes[i-1][2] if i > 0 else 'NONE',
                    'right_label': genes[i+1][2] if i < len(genes)-1 else 'NONE'
                })
                continue

            left_lab  = genes[i-1][2]
            right_lab = genes[i+1][2]

            if left_lab == 'H4' and right_lab == 'H4':
                local_A.append((sid, gidx, left_lab, right_lab))
            elif {left_lab, right_lab} == {'H1', 'H2A-IV'}:
                local_B.append((sid, gidx, left_lab, right_lab))
            else:
                rows_report.append({
                    'genome_id': genome_id, 'cluster_id': cluster_id,
                    'gene_idx': gidx, 'sseqid': sid,
                    'rule': 'NA_context_mismatch',
                    'left_label': left_lab, 'right_label': right_lab
                })

        if local_A:
            for sid, gidx, l, r in local_A:
                keepA.append(sid)
                rows_report.append({
                    'genome_id': genome_id, 'cluster_id': cluster_id,
                    'gene_idx': gidx, 'sseqid': sid,
                    'rule': 'A_H4-H2B7-H4',
                    'left_label': l, 'right_label': r
                })
        else:
            for sid, gidx, l, r in local_B:
                keepB.append(sid)
                rows_report.append({
                    'genome_id': genome_id, 'cluster_id': cluster_id,
                    'gene_idx': gidx, 'sseqid': sid,
                    'rule': 'B_H2AIV/H1',
                    'left_label': l, 'right_label': r
                })

    chosen = set(keepA) | set(keepB)
    other = [sid for sid in all_h2b7_ids if sid not in chosen]

    return {'A': keepA, 'B': keepB, 'report': rows_report, 'all_h2b7': all_h2b7_ids, 'other': other}


###############################################################################
# CLI
###############################################################################

def main():
    ap = argparse.ArgumentParser(
        description="Extract ALL H2B7 like the base script + context classification: "
                    "True(A): H4|H2B7|H4; True(B): H1/H2A-IV (fallback if no A in cluster)."
    )
    ap.add_argument('--hits_tsv', required=True,
                    help='TSV labels: either [qseqid sseqid label] (≥3 cols) or [sseqid label] (2 cols). '
                         'To classify context, labels for ALL cluster genes must be specific (H4/H1/H2A-IV/H2B7/…).')
    ap.add_argument('--hits_fa', required=True,
                    help='FASTA with sequences keyed by sseqid (e.g., hits_marked.fa).')
    ap.add_argument('--meta', required=False,
                    help='Meta table: genomic_id<TAB>species_label (or 2+ spaces). Optional.')
    ap.add_argument('--outdir', required=True, help='Output directory.')
    ap.add_argument('--prefix', default='H2B7', help='Prefix for output files (default: H2B7).')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # meta
    meta = {}
    if args.meta and os.path.isfile(args.meta):
        try:
            meta = load_meta(args.meta)
        except Exception as e:
            print(f"[WARN] meta parsing failed: {e}. Continue without species labels.", file=sys.stderr)

    # load labels and FASTA
    labels_df = load_labels_from_hits(args.hits_tsv)
    fa = read_fasta_dict(args.hits_fa)

    # === BASE funtionality
    df_h2b7 = labels_df[labels_df['label'] == 'H2B7'].copy()
    ids_all = sorted(df_h2b7['sseqid'].dropna().unique().tolist())

    combined = []
    by_species = defaultdict(list)
    missing_all = []

    for sid in ids_all:
        genome_id, cluster, gid, gene_token = parse_sseqid(sid)
        if genome_id is None:
            continue
        species = meta.get(genome_id, genome_id)
        sp_sanit = sanitize_species(species)

        if sid not in fa:
            missing_all.append(sid)
            continue

        header = f"{sp_sanit}|{genome_id}|cl{cluster}|{gid if gid is not None else gene_token}"
        seq = fa[sid]
        combined.append((header, seq))
        by_species[sp_sanit].append((header, seq))

    out_all = os.path.join(args.outdir, f"{args.prefix}_all.faa")
    write_fasta(combined, out_all)
    print(f"[OK] {out_all}: {len(combined)} seqs")

    for sp, recs in by_species.items():
        out_sp = os.path.join(args.outdir, f"{args.prefix}_{sp}.faa")
        write_fasta(recs, out_sp)
        print(f"[OK] {out_sp}: {len(recs)} seqs")

    if missing_all:
        print(f"[WARN] absent in FASTA {len(missing_all)} sseqid from main H2B7-set (first 20):", file=sys.stderr)
        for m in missing_all[:20]:
            print("   ", m, file=sys.stderr)

    # === КОНТЕКСТНАЯ КЛАССИФИКАЦИЯ TRUE vs OTHER ===
    specific_needed = {'H4', 'H1', 'H2A-IV'}
    has_specific = labels_df['label'].isin(specific_needed | {'H2B7'}).any()
    if not has_specific:
        print("[WARN] Seems, in hits_tsv for non-H2B7 there are unhelpful namings (example, 'Other'). "
              "Context classifiaction may be empty. Base files created anyway.",
              file=sys.stderr)

    ctx = select_h2b7_by_context(labels_df)

    def emit_fasta(name, id_list):
        recs, miss = [], []
        for sid in id_list:
            genome_id, cluster, gid, gene_token = parse_sseqid(sid)
            if genome_id is None:
                continue
            species = meta.get(genome_id, genome_id)
            sp_sanit = sanitize_species(species)
            if sid in fa:
                header = f"{sp_sanit}|{genome_id}|cl{cluster}|{gid if gid is not None else gene_token}"
                recs.append((header, fa[sid]))
            else:
                miss.append(sid)
        out_path = os.path.join(args.outdir, f"{args.prefix}_{name}.faa")
        if recs:
            write_fasta(recs, out_path)
            print(f"[OK] {out_path}: {len(recs)} seqs")
        else:
            print(f"[WARN] {out_path}: 0 seqs")
        if miss:
            print(f"[WARN] {name}: absent in FASTA {len(miss)} (first 10):", file=sys.stderr)
            for m in miss[:10]:
                print("   ", m, file=sys.stderr)

    # True sets
    emit_fasta("true_between_H4_H4", ctx['A'])
    emit_fasta("true_between_H2AIV_H1", ctx['B'])
    emit_fasta("true_union", list(dict.fromkeys(ctx['A'] + ctx['B'])))

    # Other (все H2B7 минус True)
    emit_fasta("other", ctx['other'])

    # Report
    report_rows = ctx['report']
    for r in report_rows:
        species = meta.get(r['genome_id'], r['genome_id'])
        r['species_label'] = species

    cols = ['species_label', 'genome_id', 'cluster_id', 'gene_idx',
            'sseqid', 'rule', 'left_label', 'right_label']
    report_df = pd.DataFrame(report_rows)
    if not report_df.empty:
        for c in cols:
            if c not in report_df.columns:
                report_df[c] = ''
        report_df = report_df[cols]
        report_df.sort_values(by=['species_label', 'genome_id', 'cluster_id', 'gene_idx', 'sseqid'],
                              inplace=True, kind='mergesort')
    report_path = os.path.join(args.outdir, f"{args.prefix}_selection_report.tsv")
    report_df.to_csv(report_path, sep='\t', index=False)
    print(f"[OK] {report_path}: {report_df.shape[0] if not report_df.empty else 0} rows")

    # === NEW: per-cluster / per-genome stats of TRUE presence =================
    true_set = set(ctx['A']) | set(ctx['B'])

    df_all_h2b7 = labels_df[labels_df['label'] == 'H2B7'].copy()
    parsed_all = df_all_h2b7['sseqid'].apply(parse_sseqid)
    df_all_h2b7[['genome_id','cluster_id','gene_idx','gene_token']] = pd.DataFrame(parsed_all.tolist(),
                                                                                   index=df_all_h2b7.index)
    df_all_h2b7 = df_all_h2b7.dropna(subset=['genome_id','cluster_id'])

    def cluster_has_true(subdf):
        sids = set(subdf['sseqid'].tolist())
        return int(len(sids & true_set) > 0)

    clust_stats = (df_all_h2b7
                   .groupby(['genome_id','cluster_id'], as_index=False)
                   .apply(lambda g: pd.Series({
                       'has_true': cluster_has_true(g),
                       'n_h2b7_in_cluster': g.shape[0]
                   }))
                   .reset_index(drop=True))

    clust_stats['species_label'] = clust_stats['genome_id'].map(lambda gid: meta.get(gid, gid))
    clust_stats = clust_stats[['species_label','genome_id','cluster_id','n_h2b7_in_cluster','has_true']]
    clust_stats.sort_values(['species_label','genome_id','cluster_id'], inplace=True, kind='mergesort')

    clust_report = os.path.join(args.outdir, f"{args.prefix}_true_stats_by_cluster.tsv")
    clust_stats.to_csv(clust_report, sep='\t', index=False)
    print(f"[OK] {clust_report}: {clust_stats.shape[0]} rows")

    genome_stats = (clust_stats
                    .groupby(['genome_id','species_label'], as_index=False)
                    .agg(n_clusters_with_h2b7=('cluster_id','nunique'),
                         n_clusters_with_true=('has_true','sum')))
    genome_stats['has_any_true'] = genome_stats['n_clusters_with_true'] > 0
    genome_stats.sort_values(['species_label','genome_id'], inplace=True, kind='mergesort')

    genome_report_all = os.path.join(args.outdir, f"{args.prefix}_true_stats_by_genome.tsv")
    genome_stats.to_csv(genome_report_all, sep='\t', index=False)
    print(f"[OK] {genome_report_all}: {genome_stats.shape[0]} rows")

    no_true = genome_stats[~genome_stats['has_any_true']].copy()
    no_true_report = os.path.join(args.outdir, f"{args.prefix}_no_true_by_genome.tsv")
    no_true.to_csv(no_true_report, sep='\t', index=False)
    print(f"[OK] {no_true_report}: {no_true.shape[0]} genomes with zero TRUE")

    # Summary
    union_true = list(dict.fromkeys(ctx['A'] + ctx['B']))
    print("\n=== SUMMARY ===")
    print(f"All H2B7 (base set): {len(ids_all)} ids (written to {os.path.basename(out_all)} + per-species)")
    print(f"Rule A (H4-H2B7-H4): {len(ctx['A'])}")
    print(f"Rule B (H2AIV/H1):   {len(ctx['B'])}   (only in clusters without Rule-A hits)")
    print(f"True union:          {len(union_true)} (written to {args.prefix}_true_union.faa)")
    print(f"Other H2B7:          {len(ctx['other'])} (written to {args.prefix}_other.faa)")
    print(f"Reports:             {os.path.basename(report_path)}, "
          f"{os.path.basename(clust_report)}, {os.path.basename(genome_report_all)}, "
          f"{os.path.basename(no_true_report)}")

if __name__ == "__main__":
    main()