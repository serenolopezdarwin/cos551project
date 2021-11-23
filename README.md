# Cos551 Final Project

Sereno Lopez-Darwin \
November 2021

Recreates analysis originally performed by Cao et al., 2019[^1].

### cos551.py

Module file that contains object classes and shared functions for my 
python scripts. 

**Current classes:** 
- Mean: Calculates rolling means, stores count and mean value.

**Current functions:**
- log(string): Initializes a timed and dated logfile if none exists, 
then writes the string to it.

### parsesourcefiles.py

Preprocessing script that filters and formats our sparse matrix. 

Filters out:

- Cells marked as doublets in author's annotations
- Cells with expression level more than 2 standard deviations 
away from the mean expression level across all cells.
- Genes that are not members of the 2000 most variant genes across 
the full dataset.

### jobsubmission.sh

Job submission script for use with SLURM. Should be modified
based on user's cluster system.

### differentialexpressionanalysis.py

Plotting script that performs differential expression analysis on 
each gene per-cluster and plots the genes with highest percentage
expression and log-fold-change from other clusters.

[^1]: Cao, J et al. The single-cell transcriptional landscape of 
mammalian organogenesis. *Nature* **566**, 496–502 (2019).