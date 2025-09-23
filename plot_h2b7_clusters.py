#!/usr/bin/env python3
"""
plot_h2b7_clusters.py

For evry assembly:
 1. Draws PNG for every cluster:
    - horizontal line form beg…end (global coords)
    - arrows for every gene (hit-hits dyed by strict/loose/remote or H2B7/H3/…)
    - above line signed contig: beg-end (cluster coords)
 2. Merges all clusters in one final PNG
 3. For final PNG adds:
    -  headear with species_label (from metadata) 
    - clolor legen from COLOR_MAP (example, strict/loose/remote or H2B7/H3/H2A-IV/H1/H4) 
"""
import argparse, os, re, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont

COLOR_MAP = {
    'H2B7':   'red',
    'H3':     'blue',
    'H2A-IV': 'green',
    'H1':     'gold',
    'H4':     'brown',
}
DEFAULT_COLOR = 'black'
DEFAULT_ALPHA = 0.4

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--hits',     required=True,  help='TSV с hit-ами H2B7')
    p.add_argument('--meta',     required=False, help='(optional) TSV: genomic_id<TAB>species_label')
    p.add_argument('--histones', required=True,  help='Root of assembly dir')
    p.add_argument('--outdir',   required=True,  help='Where to save')
    return p.parse_args()

def load_meta(path):
    """
    Format: genomic_id<TAB>species_label
            or: genomic_id<2+ sapeces>species_label
    Returns: {genomic_id -> species_label}
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

def read_hits(path):
    df = pd.read_csv(path, sep='\t', header=None, dtype=str, comment='#')
    recs = []
    for row in df.itertuples(index=False):
        subj, grp = row[1], row[-1]
        if 'genomic_' not in subj:
            continue
        asm_part, tail = subj.split('_genomic_',1)
        asm = asm_part + '_genomic'
        parts = tail.split('_')
        if len(parts) < 2:
            continue
        recs.append((asm, parts[0], parts[1], grp))
    return recs

def find_gff(root, asm, cl):
    d = os.path.join(root, asm)
    if not os.path.isdir(d):
        return None
    for fn in os.listdir(d):
        if fn.startswith(f"cluster{cl}_") and fn.endswith('.gff'):
            return os.path.join(d, fn)
    return None

def parse_gff(path):
    hdr = re.compile(r'name\s*=\s*::([^:]+):(\d+)-(\d+)')
    contig = beg = end = None
    genes = {}
    with open(path) as fh:
        for L in fh:
            if L.startswith('#'):
                m = hdr.search(L)
                if m:
                    contig, beg, end = m.group(1), int(m.group(2)), int(m.group(3))
                continue
            cols = L.strip().split('\t')
            if len(cols) < 9 or cols[2] != 'gene':
                continue
            s_loc, e_loc, strand = int(cols[3]), int(cols[4]), cols[6]
            gid = None
            for a in cols[8].split(';'):
                if a.startswith('ID='):
                    gid = a.split('=',1)[1].split('.')[0]
                    break
            if not gid:
                continue
            s_glob = beg + s_loc - 1
            e_glob = beg + e_loc - 1
            genes[gid] = (s_glob, e_glob, strand)
    if contig is None:
        raise RuntimeError(f"Region couldn`t be parsed in {path}")
    return genes, contig, beg, end

def plot_cluster_png(path, genes, hits, cl, contig, beg, end):
    fig, ax = plt.subplots(figsize=(8,1.5))
    ax.text(beg, 1.05, f"{contig}: {beg}-{end}", ha='left', va='bottom', fontsize=8)
    y = 0.5
    ax.hlines(y, beg, end, color='gray', lw=2, alpha=0.7)
    ax.set_xlim(beg, end); ax.set_ylim(0,1); ax.set_yticks([])
    rng = end - beg
    for gid,(s,e,strand) in genes.items():
        dx = e - s; start = s
        if strand == '-':
            dx, start = -dx, e
        if gid in hits:
            col = COLOR_MAP.get(hits[gid], DEFAULT_COLOR); alp = 1.0
        else:
            col = DEFAULT_COLOR; alp = DEFAULT_ALPHA
        head_l = max(abs(dx)*0.02, rng*(0.005 if gid in hits else 0.002))
        arr = FancyArrow(start, y, dx, 0,
                         length_includes_head=True,
                         head_width=0.1, head_length=head_l,
                         width=0.02, color=col, alpha=alp)
        ax.add_patch(arr)
        if gid in hits:
            ax.text(start + dx/2, y + 0.1*1.5, gid,
                    ha='center', va='bottom', fontsize=6, rotation=90, color=col)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: str(int(x))))
    plt.tight_layout(pad=0.1)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def merge_with_header(img_paths, out_path, title):
    try:
        font_h = ImageFont.truetype("DejaVuSans-Bold.ttf", 24)
        font_l = ImageFont.truetype("DejaVuSans.ttf", 16)
    except:
        font_h = ImageFont.load_default()
        font_l = ImageFont.load_default()

    imgs = [Image.open(p) for p in img_paths]
    max_w = max(im.width for im in imgs)
    total_h = sum(im.height for im in imgs)
    header_h = 80
    merged = Image.new('RGB', (max_w, total_h + header_h), 'white')
    draw = ImageDraw.Draw(merged)

    bbox = draw.textbbox((0,0), title, font=font_h)
    tw, th = bbox[2]-bbox[0], bbox[3]-bbox[1]
    draw.text(((max_w - tw)//2, 10), title, fill='black', font=font_h)

    x, y = 20, 10 + th + 10
    for name, col in COLOR_MAP.items():
        draw.rectangle([x, y, x+20, y+20], fill=col)
        draw.text((x+25, y), name, fill='black', font=font_l)
        tb = draw.textbbox((0,0), name, font=font_l)
        text_w = tb[2]-tb[0]
        x += 25 + text_w + 20

    y_off = header_h
    for im in imgs:
        merged.paste(im, (0, y_off))
        y_off += im.height

    merged.save(out_path)

import re as _re
def _safe_slug(s: str) -> str:
    s = _re.sub(r'[^\w.-]+', '_', s, flags=_re.UNICODE)
    s = _re.sub(r'_{2,}', '_', s)
    return s.strip('_')[:200] if s else 'unknown'

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = {}
    if args.meta and os.path.isfile(args.meta):
        meta = load_meta(args.meta)

    recs = read_hits(args.hits)
    asm_cl = defaultdict(lambda: defaultdict(dict))
    for asm,cl,gid,grp in recs:
        asm_cl[asm][cl][gid] = grp

    for asm, clusters in asm_cl.items():
        acc = re.match(r'^(GCA_|GCF_)\d+\.\d+', asm)
        genome_id = acc.group(0) if acc else asm
        species = meta.get(genome_id, None)
        if not species:
        
            print(f"[WARN] meta: no species_label for {genome_id}; using accession", file=sys.stderr)
            species = genome_id

        tmp = os.path.join(args.outdir, f"{genome_id}_tmp")
        os.makedirs(tmp, exist_ok=True)
        pngs = []

        for cl in sorted(clusters, key=lambda x: int(x)):
            gff = find_gff(args.histones, asm, cl)
            if not gff:
                print(f"[WARN] no GFF for {asm} cl{cl}")
                continue
            genes, contig, beg, end = parse_gff(gff)
            out = os.path.join(tmp, f"{genome_id}_cl{cl}.png")
            plot_cluster_png(out, genes, clusters[cl], cl, contig, beg, end)
            pngs.append(out)

        if not pngs:
            continue


        final_name = _safe_slug(species)
        final = os.path.join(args.outdir, f"{final_name}.png")
        merge_with_header(pngs, final, title=species)
        print(f"[OK] {final}")

        for p in pngs:
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(tmp)
        except OSError:
            pass

if __name__=='__main__':
    main()