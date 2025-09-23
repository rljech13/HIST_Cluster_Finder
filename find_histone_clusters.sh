#!/usr/bin/env bash
set -euo pipefail

#################### НАСТРОЙКИ ####################
# Root dir - were to fina all .fna or .fa files
GENOMES_ROOT="/mnt/data/i.rubinshteyn/genomes/galgal6/ncbi_dataset/data/GCF_000002315.6/GCF_000002315.6_GRCg6a_genomic.fna"

# Anchor file with your protein
ANCHOR_PROT="/mnt/data/i.rubinshteyn/genomes/galgal6/ncbi_dataset/data/GCF_000002315.6/H2B7.fa"

# Augustus parameters
AUG_SPECIES="chicken"
THREADS=96

# Padding around blast-hit (in bp)
PAD=50000

# Max distance between two hits for merging them in one cluster 
CLUSTER_DIST=20000
###################################################

# main cycle
find "$GENOMES_ROOT" -type f -name '*.fna' | while read -r GENOME; do
  ASM=$(basename "$GENOME" .fna)
  OUTDIR="$(dirname "$GENOME")/../../histones/${ASM}"
  LOG="$OUTDIR/scan_H2B7_${ASM}.log"

  mkdir -p "$OUTDIR"
  echo ">>> $(date): START scan H2B7 in ${ASM}" | tee "$LOG"

  (
    
    set -euo pipefail
    exec >> "$LOG" 2>&1

    cd "$OUTDIR"

    # 1) faidx
    if [[ ! -f "${GENOME}.fai" ]]; then
      echo "[INFO] faidx absent - creating one..." 
      samtools faidx "$GENOME"
    fi

    # 2) makeblastdb
    if [[ ! -f db.nhr ]]; then
      echo "[INFO] making BLAST DB..."
      makeblastdb -in "$GENOME" -dbtype nucl -parse_seqids -out db
    fi

    # 3) tblastn
    echo "[INFO] launching tblastn..."
    tblastn \
      -query "$ANCHOR_PROT" \
      -db db \
      -evalue 1e-20 \
      -max_hsps 5 \
      -num_threads "$THREADS" \
      -outfmt '6 sseqid sstart send pident length evalue' \
      > raw_hits.tsv

    echo "[DEBUG] raw_hits.tsv lines: $(wc -l < raw_hits.tsv)"

    # 4) filter
    awk -vOFS='\t' '$4>=70 && $5>=80' raw_hits.tsv > hits.tsv
    echo "[DEBUG] hits.tsv (pident>=70,len>=80): $(wc -l < hits.tsv)"

    # 5) into BED ±PAD
    awk -vOFS='\t' -v PAD="$PAD" '
    {
      split($1, A, "|"); chr=A[2]
      s=$2; e=$3
      if(s>e){t=s;s=e;e=t}
      start=s-1-PAD; if(start<0) start=0
      end=e+PAD
      print chr, start, end
    }
    ' hits.tsv > hits.bed

    echo "[DEBUG] hits.bed lines: $(wc -l < hits.bed)"

    # 6) merging
    echo "[INFO] sorting and merging intervals (d<=${CLUSTER_DIST})..."
    sort -k1,1 -k2,2n hits.bed > hits.sorted.bed
    bedtools merge -i hits.sorted.bed -d "$CLUSTER_DIST" > clusters.bed

    echo "[DEBUG] clusters.bed lines: $(wc -l < clusters.bed)"
    cat clusters.bed

    # 7) cluster cycles
    i=0
    while read -r CHR BEG END; do
      i=$((i+1))
      # contig length from faidx
      CONTIG_LEN=$(awk -v chr="$CHR" '$1==chr{print $2}' "${GENOME}.fai" || echo "")
      if [[ -z "$CONTIG_LEN" ]]; then
        echo "[ERROR] contig $CHR is not found in ${GENOME}.fai — skipping"
        continue
      fi
      # correction
      (( BEG<0   )) && BEG=0
      (( END>CONTIG_LEN )) && END=$CONTIG_LEN

      ID="cluster${i}_${CHR}_${BEG}_${END}"
      echo "[INFO] → $ID: $CHR:$BEG-$END (len=$CONTIG_LEN)"

      # cut and annotate
      printf "%s\t%d\t%d\n" "$CHR" "$BEG" "$END" > "${ID}.bed"
      bedtools getfasta -fi "$GENOME" -bed "${ID}.bed" -name -fo "${ID}.fa"
      if augustus \
           --species="$AUG_SPECIES" \
           --strand=both \
           --UTR=off \
           --gff3=on \
           --genemodel=partial \
           --outfile="${ID}.gff" \
           "${ID}.fa" \
           > "${ID}.aug.log" 2>&1
      then
        echo "[OK] ${ID}.gff is ready"
      else
        echo "[ERROR] Augustus failed on ${ID}"
      fi
      echo
    done < clusters.bed

    echo ">>> $(date): DONE scan ${ASM}"
  )

done

echo ">>> ALL DONE"