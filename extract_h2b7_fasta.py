#!/usr/bin/env python3
import os, sys, re, argparse
import pandas as pd
from collections import defaultdict

def load_meta(path):
    """
    Reads meta-table with format:
        genomic_id<TAB>species_label
      or
        genomic_id<2+ пробелов>species_label

    Returns dict: {genomic_id -> species_label}.
    Don't use pandas.
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
    if '_genomic_' not in sseqid:
        return None, None, None
    asm_part, tail = sseqid.split('_genomic_', 1)
    asm = asm_part + '_genomic'
    parts = tail.split('_')
    if len(parts) < 2:
        return None, None, None
    cluster = parts[0]
    gene_id = parts[1]
    m = re.match(r'^(GCA_|GCF_)\d+\.\d+', asm)
    genome_id = m.group(0) if m else asm
    return genome_id, cluster, gene_id

def read_fasta_dict(path):
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
    with open(out_path, 'w') as out:
        for h, s in records:
            out.write(f">{h}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits_tsv', required=True, help='hits_marked.tsv')
    ap.add_argument('--hits_fa',  required=True, help='hits_marked.fa')
    ap.add_argument('--meta',     required=False, help='meta TSV: genomic_id<TAB>species_label (or 2+ spaces between columns)')
    ap.add_argument('--outdir',   required=True, help='куда писать FASTA')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # meta
    meta = {}
    if args.meta and os.path.isfile(args.meta):
        try:
            meta = load_meta(args.meta)
        except Exception as e:
            print(f"[WARN] meta parsing failed: {e}. Continue without labels.", file=sys.stderr)

    df = pd.read_csv(args.hits_tsv, sep='\t', header=None, dtype=str, comment='#')
    if df.shape[1] < 3:
        sys.exit("At least 3 columns in hits_marked.tsv is needed (qseqid, sseqid, label).")
    df = df[[0,1,2]].copy()
    df.columns = ['qseqid','sseqid','label']
    df = df[df['label'] == 'H2B7']

    ids = sorted(df['sseqid'].dropna().unique().tolist())

    annot = {}
    for sid in ids:
        genome_id, cluster, gid = parse_sseqid(sid)
        if genome_id is None:
            continue
        species = meta.get(genome_id, genome_id)
        annot[sid] = (species, genome_id, cluster, gid)

    seqs = read_fasta_dict(args.hits_fa)

    combined = []
    by_species = defaultdict(list)
    missing = []

    for sid in ids:
        if sid not in annot:
            continue
        if sid not in seqs:
            missing.append(sid); continue
        species, genome_id, cluster, gid = annot[sid]
        sp_sanit = species.replace(' ','_')
        header = f"{sp_sanit}|{genome_id}|cl{cluster}|{gid}"
        seq = seqs[sid]
        combined.append((header, seq))
        by_species[sp_sanit].append((header, seq))

    out_all = os.path.join(args.outdir, "H2B7_all.faa")
    write_fasta(combined, out_all)
    print(f"[OK] {out_all}: {len(combined)} seqs")

    for sp, recs in by_species.items():
        out_sp = os.path.join(args.outdir, f"H2B7_{sp}.faa")
        write_fasta(recs, out_sp)
        print(f"[OK] {out_sp}: {len(recs)} seqs")

    if missing:
        print(f"[WARN] missing in FASTA {len(missing)} sseqid (look up for STDERR)", file=sys.stderr)
        for m in missing[:20]:
            print(f"  missing: {m}", file=sys.stderr)

if __name__ == "__main__":
    main()