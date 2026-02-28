# =========================================
# Bioinformatics 2025/2026 - Day 2
# All-in-one pipeline (Exercises 1–8)
# make clean && make all should reproduce outputs without manual commands
# =========================================

SHELL := /bin/bash

# ---------- Environment ----------
CONDA_ENV := bioinfo_lab2
PYTHON_VERSION := 3.10
CHANNEL := conda-forge
PACKAGES := python=$(PYTHON_VERSION) emboss seqkit

# ---------- Ensembl release 115 URLs ----------
# cDNA
HUMAN_CDNA_URL := https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
DMEL_CDNA_URL  := https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.54.cdna.all.fa.gz

# TOP-LEVEL DNA (unmasked)
HUMAN_DNA_URL := https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
DMEL_DNA_URL  := https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz

# Human Chromosome Y DNA (unmasked)
HUMAN_CHRY_URL := https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz

# ---------- Local paths ----------
# Downloaded inputs
HUMAN_CDNA := data/human_cdna.fa.gz
DMEL_CDNA  := data/dmel_cdna.fa.gz
HUMAN_DNA  := data/human_dna.fa.gz
DMEL_DNA   := data/dmel_dna.fa.gz
HUMAN_CHRY := data/human_chry.fa.gz

# Exercise 5 outputs (fast exploratory analysis)
HUMAN_CDNA_NUMNUC := results/human_cdna.num_nucleotides.txt
DMEL_CDNA_NUMNUC  := results/dmel_cdna.num_nucleotides.txt

HUMAN_CDNA_NUMSEQ := results/human_cdna.num_seqs.txt
DMEL_CDNA_NUMSEQ  := results/dmel_cdna.num_seqs.txt

HUMAN_CDNA_NUMN   := results/human_cdna.num_N.txt
DMEL_CDNA_NUMN    := results/dmel_cdna.num_N.txt

# Exercise 5.4/5.5 outputs
HUMAN_TRANSFORMED := data/human_cdna.transformed.fa.gz
DMEL_TRANSFORMED  := data/dmel_cdna.transformed.fa.gz
HUMAN_MASKED      := data/human_cdna.masked_polyA.fa.gz
DMEL_MASKED       := data/dmel_cdna.masked_polyA.fa.gz

# Exercise 5.6 output
ECORI_COUNT_OUT := results/GAATTC.num.txt

# Exercise 6 outputs (Python vs EMBOSS GC profiling on cDNA)
GC_PY_HUMAN  := results/human_trans.gc.txt
GC_PY_DMEL   := results/dmel_trans.gc.txt
GC_EMB_HUMAN := results/human_trans.gc.emboss.txt
GC_EMB_DMEL  := results/dmel_trans.gc.emboss.txt
GC_COMPARE   := results/gc_compare.txt

# Exercise 7 output
BENCH_TSV := results/gc_bench_table.tsv

# Exercise 8 outputs
HUMAN_CHRY_ORFS := results/human_chry.orfs.fa
HUMAN_CDNA_ORFS := results/human_cdna.orfs.fa
ORF_STATS_TXT   := results/orf_stats.txt

# ---------- Phony ----------
.PHONY: all env data results download_data ex5 ex6 ex7 ex8 transform mask_polyA ecori gc_py gc_emboss gc_stats bench clean

# ---------- Main ----------
all: env download_data ex5 transform mask_polyA ecori ex6 ex7 ex8

# ---------- Directories ----------
data:
	@mkdir -p data

results:
	@mkdir -p results

# ---------- Environment ----------
env:
	@echo "==> Ensuring conda env exists: $(CONDA_ENV)"
	@conda env list | awk '{print $$1}' | grep -qx "$(CONDA_ENV)" || \
		conda create -n $(CONDA_ENV) -c $(CHANNEL) $(PACKAGES) -y
	@echo "==> Environment ready."
	@chmod +x gc_stats.py || true

# ---------- Downloads (Exercises 3 & 4) ----------
# NOTE: Ensembl FTP for drosophila_melanogaster/dna/ does not provide dna.chromosome.* files.
# For "chromosome-level" DNA we use DMEL_DNA (dna.toplevel).
download_data: data $(HUMAN_CDNA) $(DMEL_CDNA) $(HUMAN_DNA) $(DMEL_DNA) $(HUMAN_CHRY)
	@echo "==> All files downloaded."

$(HUMAN_CDNA):
	curl -L $(HUMAN_CDNA_URL) -o $@

$(DMEL_CDNA):
	curl -L $(DMEL_CDNA_URL) -o $@

$(HUMAN_DNA):
	curl -L $(HUMAN_DNA_URL) -o $@

$(DMEL_DNA):
	curl -L $(DMEL_DNA_URL) -o $@

$(HUMAN_CHRY):
	curl -L $(HUMAN_CHRY_URL) -o $@

# ============================================================
# Exercise 5.1 / 5.2 / 5.3 (Fast exploratory analysis on cDNA)
# ============================================================

ex5: $(HUMAN_CDNA_NUMNUC) $(DMEL_CDNA_NUMNUC) $(HUMAN_CDNA_NUMSEQ) $(DMEL_CDNA_NUMSEQ) $(HUMAN_CDNA_NUMN) $(DMEL_CDNA_NUMN)
	@echo "==> Exercise 5 stats ready."

# 5.1 total nucleotides in cDNA (exclude headers and newlines)
$(HUMAN_CDNA_NUMNUC): results $(HUMAN_CDNA)
	zcat $(HUMAN_CDNA) | grep -v "^>" | tr -d "\n" | wc -c > $@

$(DMEL_CDNA_NUMNUC): results $(DMEL_CDNA)
	zcat $(DMEL_CDNA) | grep -v "^>" | tr -d "\n" | wc -c > $@

# 5.2 total number of sequences (transcripts)
$(HUMAN_CDNA_NUMSEQ): results $(HUMAN_CDNA)
	zgrep -c "^>" $(HUMAN_CDNA) > $@

$(DMEL_CDNA_NUMSEQ): results $(DMEL_CDNA)
	zgrep -c "^>" $(DMEL_CDNA) > $@

# 5.3 total number of 'N' characters across transcriptome
$(HUMAN_CDNA_NUMN): results $(HUMAN_CDNA)
	zcat $(HUMAN_CDNA) | grep -v "^>" | tr -cd "Nn" | wc -c > $@

$(DMEL_CDNA_NUMN): results $(DMEL_CDNA)
	zcat $(DMEL_CDNA) | grep -v "^>" | tr -cd "Nn" | wc -c > $@

# ============================================================
# Exercise 5.4 Transform (headers w/o spaces + lowercase) -> .transformed.fa.gz
# ============================================================

transform: $(HUMAN_TRANSFORMED) $(DMEL_TRANSFORMED)

%.transformed.fa.gz: %.fa.gz
	zcat $< | \
	awk 'BEGIN{OFS=""} \
	/^>/ {split($$0,a," "); print a[1]; next} \
	{print tolower($$0)}' | \
	gzip > $@

# ============================================================
# Exercise 5.5 Soft mask polyA (A{5,} -> lowercase) -> .masked_polyA.fa.gz
# ============================================================

mask_polyA: $(HUMAN_MASKED) $(DMEL_MASKED)

%.masked_polyA.fa.gz: %.fa.gz
	zcat $< | \
	awk '/^>/ {print; next} \
	{ \
	  while (match($$0,/A{5,}/)) { \
	    seq=substr($$0,RSTART,RLENGTH); \
	    lower=tolower(seq); \
	    $$0=substr($$0,1,RSTART-1) lower substr($$0,RSTART+RLENGTH); \
	  } \
	  print; \
	}' | \
	gzip > $@

# ============================================================
# Exercise 5.6 EcoRI count in Human chrY (expected 7096)
# ============================================================

ecori: $(ECORI_COUNT_OUT)

$(ECORI_COUNT_OUT): results $(HUMAN_CHRY)
	zcat $(HUMAN_CHRY) | \
	awk 'BEGIN{seq=""} \
	/^>/ {next} \
	{gsub(/[ \r\n]/,"",$$0); seq=seq $$0} \
	END{motif="GAATTC"; n=0; \
	    for(i=1;i<=length(seq)-length(motif)+1;i++){ \
	      if(substr(seq,i,length(motif))==motif) n++; \
	    } \
	    print n; \
	}' > $@

# ============================================================
# Exercise 6 (GC profiling: Python + EMBOSS + comparison)
# ============================================================

ex6: gc_stats $(GC_COMPARE)
	@echo "==> Exercise 6 outputs ready."

# 6.2 Python GC on cDNA
gc_py: $(GC_PY_HUMAN) $(GC_PY_DMEL)

$(GC_PY_HUMAN): results gc_stats.py $(HUMAN_CDNA)
	zcat $(HUMAN_CDNA) | ./gc_stats.py > $@

$(GC_PY_DMEL): results gc_stats.py $(DMEL_CDNA)
	zcat $(DMEL_CDNA) | ./gc_stats.py > $@

# 6.2 EMBOSS GC on cDNA (quiet)
gc_emboss: $(GC_EMB_HUMAN) $(GC_EMB_DMEL)

$(GC_EMB_HUMAN): results $(HUMAN_CDNA)
	geecee -auto -sequence <(zcat $(HUMAN_CDNA)) -outfile $@ >/dev/null

$(GC_EMB_DMEL): results $(DMEL_CDNA)
	geecee -auto -sequence <(zcat $(DMEL_CDNA)) -outfile $@ >/dev/null

gc_stats: gc_py gc_emboss
	@echo "==> GC stats ready (Python + EMBOSS)."

# 6.3 comparison summary into a text file (short + concrete)
$(GC_COMPARE): results $(GC_PY_HUMAN) $(GC_PY_DMEL) $(GC_EMB_HUMAN) $(GC_EMB_DMEL)
	@echo "Python weighted mean (human_cdna) = $$(tail -n 1 $(GC_PY_HUMAN))" >  $@
	@echo "Python weighted mean (dmel_cdna)  = $$(tail -n 1 $(GC_PY_DMEL))"  >> $@
	@echo 'EMBOSS arithmetic mean (human_cdna) = '$$(awk 'NF==2 && $$2 ~ /^[0-9.]+$$/ {sum+=$$2;n++} END{if(n>0) print sum/n; else print "NA"}' $(GC_EMB_HUMAN)) >> $@
	@echo 'EMBOSS arithmetic mean (dmel_cdna)  = '$$(awk 'NF==2 && $$2 ~ /^[0-9.]+$$/ {sum+=$$2;n++} END{if(n>0) print sum/n; else print "NA"}' $(GC_EMB_DMEL))  >> $@
	@echo "" >> $@
	@echo "Note: Python is weighted by nucleotide count; EMBOSS mean above is unweighted across sequences." >> $@

# ============================================================
# Exercise 7 (Benchmarking table TSV, clean)
# ============================================================

ex7: bench
	@echo "==> Exercise 7 benchmark ready."

bench: results $(HUMAN_CDNA) $(DMEL_CDNA) $(HUMAN_CHRY) $(DMEL_DNA) gc_stats.py
	@echo -e "Tool\tInput\tTime(s)\tPeakMem(KB)" > $(BENCH_TSV)

	@/usr/bin/time -f "gc_stats.py\thuman_cdna.fa.gz\t%e\t%M" \
		bash -lc 'zcat "$(HUMAN_CDNA)" | ./gc_stats.py > /dev/null' 2>> $(BENCH_TSV)
	@/usr/bin/time -f "EMBOSS\thuman_cdna.fa.gz\t%e\t%M" \
		bash -lc 'geecee -auto -sequence <(zcat "$(HUMAN_CDNA)") -outfile /dev/null >/dev/null' 2>> $(BENCH_TSV)

	@/usr/bin/time -f "gc_stats.py\tdmel_cdna.fa.gz\t%e\t%M" \
		bash -lc 'zcat "$(DMEL_CDNA)" | ./gc_stats.py > /dev/null' 2>> $(BENCH_TSV)
	@/usr/bin/time -f "EMBOSS\tdmel_cdna.fa.gz\t%e\t%M" \
		bash -lc 'geecee -auto -sequence <(zcat "$(DMEL_CDNA)") -outfile /dev/null >/dev/null' 2>> $(BENCH_TSV)

	@/usr/bin/time -f "gc_stats.py\thuman_chry.fa.gz\t%e\t%M" \
		bash -lc 'zcat "$(HUMAN_CHRY)" | ./gc_stats.py > /dev/null' 2>> $(BENCH_TSV)
	@/usr/bin/time -f "EMBOSS\thuman_chry.fa.gz\t%e\t%M" \
		bash -lc 'geecee -auto -sequence <(zcat "$(HUMAN_CHRY)") -outfile /dev/null >/dev/null' 2>> $(BENCH_TSV)

	@/usr/bin/time -f "gc_stats.py\tdmel_dna.fa.gz\t%e\t%M" \
		bash -lc 'zcat "$(DMEL_DNA)" | ./gc_stats.py > /dev/null' 2>> $(BENCH_TSV)
	@/usr/bin/time -f "EMBOSS\tdmel_dna.fa.gz\t%e\t%M" \
		bash -lc 'geecee -auto -sequence <(zcat "$(DMEL_DNA)") -outfile /dev/null >/dev/null' 2>> $(BENCH_TSV)

	@echo "==> Bench written to $(BENCH_TSV)"

# ============================================================
# Exercise 8 (getorf + seqkit stats)
# ============================================================

ex8: $(HUMAN_CHRY_ORFS) $(HUMAN_CDNA_ORFS) $(ORF_STATS_TXT)
	@echo "==> Exercise 8 ORFs + stats ready."

# 8.1 getorf on human chrY (min 300 nt, find between START/STOP, output FASTA)
$(HUMAN_CHRY_ORFS): results $(HUMAN_CHRY)
	zcat $(HUMAN_CHRY) > results/human_chry.fa.tmp
	conda run -n $(CONDA_ENV) \
		getorf -sequence results/human_chry.fa.tmp -minsize 300 -find 1 -outseq $@.tmp >/dev/null
	rm -f results/human_chry.fa.tmp
	mv $@.tmp $@

$(HUMAN_CDNA_ORFS): results $(HUMAN_CDNA)
	zcat $(HUMAN_CDNA) > results/human_cdna.fa.tmp
	conda run -n $(CONDA_ENV) \
		getorf -sequence results/human_cdna.fa.tmp -minsize 300 -find 1 -outseq $@.tmp >/dev/null
	rm -f results/human_cdna.fa.tmp
	mv $@.tmp $@

# 8.2 summary stats: number of ORFs + longest ORF length (aa) and comparison
$(ORF_STATS_TXT): results $(HUMAN_CHRY_ORFS) $(HUMAN_CDNA_ORFS)
	@echo "seqkit seq -n results/human_chry.orfs.fa | wc -l" > $@
	@echo "> $$(seqkit seq -n $(HUMAN_CHRY_ORFS) | wc -l)  (total ORFs on human chrY)" >> $@
	@echo "" >> $@

	@echo "seqkit fx2tab -nl results/human_chry.orfs.fa | sort -k2,2nr | head -n 1 | cut -f2" >> $@
	@echo "> $$(seqkit fx2tab -nl $(HUMAN_CHRY_ORFS) | sort -k2,2nr | head -n 1 | cut -f2)  (longest ORF length in aa on human chrY)" >> $@
	@echo "" >> $@

	@echo "seqkit seq -n results/human_cdna.orfs.fa | wc -l" >> $@
	@echo "> $$(seqkit seq -n $(HUMAN_CDNA_ORFS) | wc -l)  (total ORFs in human cDNA)" >> $@

# ---------- Clean ----------
clean:
	rm -rf data results
	@echo "==> Project cleaned."