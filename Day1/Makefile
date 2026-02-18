SHELL := /bin/bash

FASTA_IN   := lambda_phage.fasta
PY_OUT     := lambda_phage.revcomp.fasta
EMBOSS_OUT := lambda_phage.emboss.fasta

.PHONY: all download revcomp emboss validate clean

all: validate

download: $(FASTA_IN)

$(FASTA_IN):
	esearch -db nuccore -query "NC_001416" | efetch -format fasta > $@

revcomp: $(PY_OUT)

$(PY_OUT): revcomp.py $(FASTA_IN)
	python revcomp.py $(FASTA_IN) > $@

emboss: $(EMBOSS_OUT)

$(EMBOSS_OUT): $(FASTA_IN)
	revseq -sequence $(FASTA_IN) -outseq $@ -notag

validate: $(PY_OUT) $(EMBOSS_OUT)
	@py_md5=$$(grep -v '^>' $(PY_OUT)     | tr -d '\n' | md5sum | awk '{print $$1}'); \
	em_md5=$$(grep -v '^>' $(EMBOSS_OUT) | tr -d '\n' | md5sum | awk '{print $$1}'); \
	if [ "$$py_md5" = "$$em_md5" ]; then \
		echo "Validation Passed"; \
	else \
		echo "Validation Failed"; \
		exit 1; \
	fi
