#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path

def parse_gff_for_proteins(gff_path: Path, root: Path):
    """
    From one GFF finds pairs (header, protein_seq),
    even if sequence is divided for several lines.
    Header = ASM_clusterIdx_geneID (shorter than 50 symbols)
    """
    rel = gff_path.relative_to(root)
    
    asm = rel.parts[0]
   
    m = re.match(r"cluster(\d+)_", gff_path.stem)
    cluster_idx = m.group(1) if m else "0"

    records = []
    gene_id = None
    in_prot = False
    seq_parts = []

    with open(gff_path, 'r') as fh:
        for line in fh:
            if line.startswith('# start gene'):
                m2 = re.match(r"# start gene (\S+)", line)
                gene_id = m2.group(1) if m2 else None
                in_prot = False
                seq_parts = []

            elif line.startswith('# protein sequence'):
                in_prot = True
                after = line.partition('[')[2].rstrip('\n').rstrip(']')
                seq_parts = [after]

            elif in_prot and line.startswith('#'):
                part = line.lstrip('#').strip()
                if part.endswith(']'):
                    part = part.rstrip(']')
                    seq_parts.append(part)
                    full_seq = ''.join(seq_parts)
                    header = f"{asm}_{cluster_idx}_{gene_id}"
                    records.append((header, full_seq))
                    in_prot = False
                    seq_parts = []
                    gene_id = None
                else:
                    seq_parts.append(part)
    return records

def main():
    p = argparse.ArgumentParser(
        description="Parses protein (multiline!) from Augustus-GFF"
    )
    p.add_argument('root_histones', type=Path,
                   help="Root dir with subdirs and GFF")
    p.add_argument('-o','--out', type=Path, default=Path('all_proteins.fa'),
                   help="Output FASTA")
    args = p.parse_args()

    gffs = list(args.root_histones.rglob('*.gff'))
    if not gffs:
        print(f"[ERROR] couldnt find any .gff in {args.root_histones}", file=sys.stderr)
        sys.exit(1)

    total = 0
    with args.out.open('w') as outfh:
        for gff in gffs:
            recs = parse_gff_for_proteins(gff, args.root_histones)
            if not recs:
                print(f"[DEBUG] in {gff} can`t find proteins", file=sys.stderr)
            for hdr, seq in recs:
                outfh.write(f">{hdr}\n")
                for i in range(0, len(seq), 60):
                    outfh.write(seq[i:i+60] + "\n")
                total += 1

    print(f"[OK] written {total} proteins in {args.out}")

if __name__ == '__main__':
    main()