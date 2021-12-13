# Cos551 Final Project

Sereno Lopez-Darwin \
November 2021

Recreates analysis originally performed by Cao et al., 2019[^1].

### cos551.py

Module file that contains object classes and shared functions for my 
python scripts. 

**Current classes:** 
- Mean: Calculates rolling means, stores count and mean value.
- Bidict: Bi-directional hash table with dictionary structure.

**Current functions:**
- log(string): Initializes a timed and dated logfile if none exists, 
then writes the string to it.
- log_time(string): Decorator that logs the time it takes for the
decorated function to execute, as well as a string that can be used
to denote function task.
- invert_hash(dict): Builds a Bidict class bidirectional hash given
a python dictionary object.
- aggregate_expression_level(*args): Given a sparse expression matrix,
sorts it by either cell or gene id and returns a list of tuples with 
each index corresponding to appropriate id and each tuple being the
other index and expression level.

### parsesourcefiles.py

Preprocessing script that filters and formats our sparse matrix. 

Filters out:

- Cells marked as doublets in author's annotations
- Cells with expression level more than 2 standard deviations 
away from the mean expression level across all cells.
- Genes that are not members of the 2000 most variant genes across 
the full dataset.

### job.sh

Job to be submitted to SLURM. Should be modified based on 
user's cluster system.

### differentialexpressionanalysis.py

Plotting script that performs differential expression analysis on 
each gene per-cluster and plots the genes with highest percentage
expression and log-fold-change from other clusters. Also generates
t-SNE and UMAP plots based on Cao et al.'s original annotations.

[^1]: Cao, J et al. The single-cell transcriptional landscape of 
mammalian organogenesis. *Nature* **566**, 496–502 (2019).