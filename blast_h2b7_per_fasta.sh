#!/usr/bin/env bash
set -euo pipefail

# Usage: script <fasta-file|dir> <REF_MULTIFA.faa> <out_dir> [threads]
if (( $# < 3 )); then
  echo "Usage: $0 <fasta-file|dir> <REF_MULTIFA.faa> <out_dir> [threads]" >&2
  exit 1
fi

FA_INPUT="$1"
REF="$2"          # multi - Fasta
OUTDIR="$3"
THREADS="${4:-8}"

# tool checkup
for tool in makeblastdb blastp seqtk awk; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

# Collect Fasta
if [[ -f "$FA_INPUT" ]]; then
  FASTAS=( "$FA_INPUT" )
elif [[ -d "$FA_INPUT" ]]; then
  shopt -s nullglob
  FASTAS=( "$FA_INPUT"/*.fa "$FA_INPUT"/*.faa "$FA_INPUT"/*.fasta )
  shopt -u nullglob
else
  echo "ERROR: '$FA_INPUT' not a file or dir" >&2
  exit 1
fi

mkdir -p "$OUTDIR"


QOUT="qseqid sseqid pident length mismatch gapopen evalue bitscore"
if blastp -help 2>/dev/null | grep -q "qcovhsp"; then
  OUTFMT_FIELDS="qseqid sseqid pident length mismatch gapopen qcovhsp evalue bitscore"
else
  echo "[INFO] qcovhsp не поддерживается этой версией BLAST; используем формат без него"
  OUTFMT_FIELDS="$QOUT"
fi

sanitize_fasta() {
  # $1=in.fa  $2=out.fa  $3=map.tsv(new_id \t original_header_first_token)
  local IN="$1" OUT="$2" MAP="$3"
  awk '
    BEGIN{ OFS="\t"; }
    /^>/ {
      hdr=$0; sub(/^>/,"",hdr);
      split(hdr, a, /[ \t]/); orig=a[1];
      ++n;
      newid=sprintf("seq%08d", n);   # всегда <= 11 символов
      print newid, orig >> MAP;
      print ">"newid >> OUT;
      next
    }
    { print >> OUT }
  ' MAP="$MAP" OUT="$OUT" "$IN"
}

for FA in "${FASTAS[@]}"; do
  [[ -f "$FA" ]] || continue
  NAME=$(basename "${FA%.*}")
  WORK="$OUTDIR/$NAME"
  mkdir -p "$WORK/logs"

  echo ">> Processing $NAME"

  DB_BASENAME="$WORK/proteins_db"
  DB_FASTA="$FA"
  MAP_TSV=""

  # 1) makeblastdb 
  if ! makeblastdb -in "$DB_FASTA" -dbtype prot -parse_seqids -out "$DB_BASENAME" >"$WORK/logs/makeblastdb.log" 2>&1; then
    if grep -q "local id is too long" "$WORK/logs/makeblastdb.log"; then
      echo "[WARN] long FASTA IDs detected; sanitizing headers to <=50 chars with mapping…"
      DB_FASTA_SAN="$WORK/proteins_sanit.fa"
      MAP_TSV="$WORK/id_map.tsv"
      :> "$DB_FASTA_SAN"; :> "$MAP_TSV"
      sanitize_fasta "$FA" "$DB_FASTA_SAN" "$MAP_TSV"

      DB_FASTA="$DB_FASTA_SAN"
      if ! makeblastdb -in "$DB_FASTA" -dbtype prot -parse_seqids -out "$DB_BASENAME" >"$WORK/logs/makeblastdb.sanit.log" 2>&1; then
        echo "[ERROR] makeblastdb failed even after sanitizing. See logs." >&2
        exit 1
      fi
    else
      echo "[ERROR] makeblastdb failed for $FA (see $WORK/logs/makeblastdb.log)." >&2
      exit 1
    fi
  fi

  # 2) blastp for every query from $REF
  FMT="6 $OUTFMT_FIELDS"
  blastp \
    -query "$REF" \
    -db "$DB_BASENAME" \
    -outfmt "$FMT" \
    -max_hsps 1 \
    -max_target_seqs 5000 \
    -num_threads "$THREADS" \
    >"$WORK/hits_all.tsv" \
    2>"$WORK/logs/blastp.log"

  if [[ -n "$MAP_TSV" && -s "$MAP_TSV" ]]; then
    awk -F'\t' 'BEGIN{OFS="\t"}
      FNR==NR { map[$1]=$2; next }
      { if($2 in map){ $2=map[$2] } print }
    ' "$MAP_TSV" "$WORK/hits_all.tsv" > "$WORK/hits_all.orig.tsv"
    mv "$WORK/hits_all.orig.tsv" "$WORK/hits_all.tsv"
  fi

  # 3) filter by e-value
  EV_COL=$(awk -v hdr="$OUTFMT_FIELDS" '
    BEGIN{
      n=split(hdr, a, " ");
      for(i=1;i<=n;i++) if(a[i]=="evalue") print i;
    }')
  awk -F'\t' -v c="$EV_COL" '($c+0) <= 1e-3' "$WORK/hits_all.tsv" > "$WORK/hits_filt.tsv" || true

  # 4) mapp qseqid -> histine type (H2B7, H3, H2A-IV, H1, H4)
  awk -F'\t' -v c="$EV_COL" '
    BEGIN{
      OFS="\t";
      map["NP_001073188.1"]="H2B7";   # Gallus gallus H2B7
      map["XP_001232833.1"]="H3";
      map["XP_425455.2"]="H2A-IV";
      map["XP_425456.1"]="H1";
      map["XP_425458.2"]="H4";
    }
    ($c+0) <= 1e-3 && map[$1] != "" { print $1, $2, map[$1] }
  ' "$WORK/hits_all.tsv" > "$WORK/hits_marked.tsv"

  # 5)  
  cut -f2 "$WORK/hits_marked.tsv" | sort -u > "$WORK/ids_marked.txt"
  if [[ -s "$WORK/ids_marked.txt" ]]; then
    seqtk subseq "${MAP_TSV:+$FA}" "$WORK/ids_marked.txt" > "$WORK/hits_marked.fa"
  fi

  echo ">> Done $NAME; outputs in $WORK/"
  echo "   - hits_all.tsv      : all BLAST hits (with orig sseqid)"
  echo "   - hits_filt.tsv     : e<=1e-3"
  echo "   - hits_marked.tsv   : qseqid, sseqid, {H2B7|H3|H2A-IV|H1|H4}"
  echo "   - hits_marked.fa    : (opt.) target sequences"
  [[ -n "${MAP_TSV:-}" ]] && echo "   - id_map.tsv         : if sanitized: new_id -> orig_header"
done