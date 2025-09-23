#!/usr/bin/env bash
# collect_H2B7.sh  <histones_root>
# 2025-07-14  –  версия с фильтрами: cov≥95 %, ≤4 mm, len 120-140 aa

set -euo pipefail
trap 'echo "ERR: «$BASH_COMMAND» (строка $LINENO)" >&2' ERR

echo "RUNNING collect_H2B7 v2025-07-14" >&2

###############################################################################
#tunabale params
ROOT=$(realpath "${1:?write down dir for ROOT/histones}")
QUERY="/mnt/data/i.rubinshteyn/genomes/galgal6/ncbi_dataset/data/GCF_000002315.6/H2B7.fa"

MIN_COV=95      # min coverage
MAX_MM=4        
MIN_LEN=120    
MAX_LEN=140

LOG="$ROOT/collect_H2B7.log" ; : >"$LOG"
OUT="$ROOT/all_H2B7.fa"      ; : >"$OUT"
TMP=$(mktemp -d)

###############################################################################
extract_prot(){          # $1=gff → stdout=FAA
  awk '
    BEGIN{RS="\n"; FS=""}
    /^# start gene/      {gid=$4; n=0}
    /^# protein sequence/{
        seq=""; sub(/^.*\[/,"",$0)
        while (1){
            gsub(/\302\240| /,""); gsub(/^#[[:space:]]*/,"")
            p=index($0,"]")
            if(p){seq=seq substr($0,1,p-1); break}
            seq=seq $0; getline
        }
        if(seq!=""){print ">"gid"_"++n"\n"seq}
    }' "$1"
}

###############################################################################
echo "### $(date)  start" | tee -a "$LOG"

#--- collecting all cluster.gff
mapfile -t GFFS < <(find "$ROOT" -type f -name 'cluster*.gff' | sort)
(( ${#GFFS[@]} > 0 )) || { echo "Нет *.gff – выхожу" ; exit 0; }

tot=0; found=0
for GFF in "${GFFS[@]}"; do
  (( ++tot ))
  BAS=$(basename "$GFF")
  echo ">> $BAS" | tee -a "$LOG"

  DIR=$(dirname  "$GFF")
  ASM=$(basename "$DIR")
  GROUP=$(basename "$(dirname "$DIR")")   
  CLU=${BAS%.gff}

  PROT="$TMP/$CLU.faa"
  extract_prot "$GFF" >"$PROT"
  [[ -s $PROT ]] || { echo "   …нет protein sequence" | tee -a "$LOG"; continue; }

  #--- local blast
  makeblastdb -in "$PROT" -dbtype prot -logfile /dev/null

  hit=$(blastp -query "$QUERY" -db "$PROT" \
        -evalue 1e-5 -outfmt "6 sseqid pident length qlen" -max_target_seqs 5 2>/dev/null |
        awk -v cov=$MIN_COV -v mm=$MAX_MM -v lo=$MIN_LEN -v hi=$MAX_LEN '
          {
            mism = $3 - int($3 * $2 / 100 + 0.5)
            pcov = 100 * $3 / $4
            if (pcov >= cov && mism <= mm && $3 >= lo && $3 <= hi)
              {print $1; exit}
          }')

  if [[ -z $hit ]]; then
      echo "   …H2B7 ниже порога" | tee -a "$LOG"
      continue
  fi

  #--- write down final fasta
  awk -v id=">$hit" \
      -v hdr=">${GROUP}|${ASM}|${CLU}|${hit}" \
      '$0==id{print hdr; getline; print; exit}' "$PROT" >>"$OUT"

  echo "   + added $hit  (cov≥${MIN_COV}, ≤${MAX_MM} mm)" | tee -a "$LOG"
  (( ++found ))
done

rm -rf "$TMP"

echo -e "\n=== summary ==="        | tee -a "$LOG"
echo "GFF-files : $tot"         | tee -a "$LOG"
echo "H2B7       : $found"       | tee -a "$LOG"
echo "FASTA → $OUT"
echo "LOG   → $LOG"