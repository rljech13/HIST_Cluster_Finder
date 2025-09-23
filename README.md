# HIST_Cluster_Finder

A smal pipeline to identify and analyze histone clusters (with focus on **H2B7**) across vertebrate genomes.  
The workflow was originally tested on chicken (Gallus gallus), then extended to all avians, crocodile and other taxa.

---

You can use environment.yml file to create a conda env that i used

## How it works


**0. Initial dataset**  
- Identify histone cluster in chicken genome.  
- Extract H2B7 and ±2 neighboring proteins.  

**1. Find histone clusters**  
```bash
bash find_histone_clusters.sh
```

Edit directories inside the script before running.

2. Parse proteins from GFF

python3 parse_gff_proteins.py <genome_histone_dir> -o chick_check.fa

3. Identify H2B7 with BLAST

bash blast_h2b7_per_fasta.sh chick_check.fa cluster5_genes.faa check_dir

4. Curate metadata manually.

5. Plot clusters

python3 plot_h2b7_clusters.py --hits hits_marked.tsv --histones histones_dir --outdir check_dir

6. Filter contextual H2B7 FASTA

python3 filter_h2b7_by_context.py \
  --hits_tsv hits_marked.tsv \
  --hits_fa  hits_marked.fa \
  --meta     meta.tsv \
  --outdir   h2b7_fasta_ctx \
  --prefix   peng_ctx

! For crocodiles, synteny-based search fails (closest homolog is not H2B7).

7. Manual curation
Keep only real H2B7 sequences and double-check genomes.

8. Add query histones
Extend curated FASTA with candidate H2B7s.

9. Multiple sequence alignment (MSA)

mafft --localpair --maxiterate 1000 --thread 16 input.faa > aligned.faa

10. Trim alignment (Trimal)

trimal -in aligned.faa -out aligned.trim.faa -automated1

11. Visualize MSA
(e.g., MView).

12. Curate representatives per species
Merge into all_for_tree.fa.

13. Clean FASTA headers

awk '/^>/{gsub(/[ \t()]/,"_"); print; next} {print}' all_h2b7.faa \
  | awk '/^>/{h=$0; if (seen[h]++) {printf("%s_DUP%d\n",h,seen[h]);} else {print h}} !/^>/' \
  > all_h2b7.clean.faa

14. Final alignment & trimming

mafft --localpair --maxiterate 1000 --thread -1 all_h2b7.clean.faa > all_h2b7.mafft.faa
trimal -in all_h2b7.mafft.faa -out all_h2b7.trim.faa -automated1

15. Build phylogenetic tree

iqtree -s all_h2b7.trim.faa -m MFP -bb 1000 -alrt 1000 -nt 32 -seed 42 -pre all_h2b7

16. Root tree

nw_reroot -l all_h2b7.treefile "H2B7_Gallus_gallus" > all_h2b7.rooted.nwk

17. Visualization
Prepare ITOL metadata for annotation and beauty.

18. Synteny and exact match counts
Performed via awk one-liners (see scripts for details).

⸻

Notes
	•	Steps 4 & 6 require manual curation.
	•	Chicken genome serves as the anchor (reference H2B7).
	•	For taxa with weak synteny conservation (e.g., crocodiles), BLAST-only approaches may fail.

⸻

Requirements
	•	Bash, Python 3
	•	MAFFT
	•	Trimal
	•	IQ-TREE
	•	Newick Utilities

⸻

Citation

If you use this pipeline, please cite this repository:
https://github.com/rljech13/HIST_Cluster_Finder
