"""
Given a filtered cell-by-gene expression sparse matrix, calculates the most enriched genes in each of our defined
cell types and plots them in percentage-by-log-fold-change dot plots.
"""
# Imports
from cos551 import *
import numpy as np
import os
import pickle as pkl
# Globals
CELL_TYPES = ["Connective tissue progenitors", "Chondrocytes and osteoblasts", "Intermediate mesoderm",
              "Jaw and tooth progenitors", "Excitatory neurons", "Epithelial cells", "Radial glia",
              "Early mesenchyme", "Neural progenitor cells", "Postmitotic premature neurons",
              "Oligodendrocyte progenitors", "Isthmic organizer cells", "Myocytes", "Dorsal neural tube cells",
              "Inhibitory neurons", "Stromal cells", "Osteoblasts", "Inhibitory neuron progenitors",
              "Premature oligodendrocytes", "Endothelial cells", "Chondrocyte progenitors",
              "Definitive erythrocyte lineage", "Schwann cell precursors", "Sensory neurons", "Limb mesenchyme",
              "Primitive erythroid lineage", "Inhibitory interneurons", "Granule neurons", "Hepatocytes",
              "Notochord and floor plate cells", "White blood cells", "Ependymal cells", "Cholinergic neurons",
              "Cardiac muscle lineage", "Megakaryocytes", "Melanocytes", "Lens", "Neutrophils"]


def build_gene_bidict(gene_ids: list) -> Bidict:
    """Given a list of gene ids, builds a bidirectional dictionary linking the 2000 gene ids remaining after filtration
    to the indices 0:1999."""
    gene_key = Bidict({})
    for idx, gene_id in enumerate(gene_ids):
        gene_key[gene_id] = idx
    return gene_key


def swap_gene_indexing(gene_ids: list, key: Bidict, to="standard") -> list:
    """Given an array of gene ids, translates each id either to indices between 0 and 1999 ('standard') or back to the
    original gene ids to be used with our data dictionary ('original')."""
    if to == "standard":
        applied_key = key
    elif to == "original":
        applied_key = key.inverse
    else:
        raise KeyError("Incorrect 'to' argument passed to gene id swapping function.")
    new_ids = []
    for gene_id in gene_ids:
        new_id = applied_key[gene_id]
        if to == "standard":
            new_ids.append(new_id)
        # The inverted dictionary's values are lists, so we need to grab the first item.
        else:
            new_ids += new_id[0]
    return new_ids


@log_time("Gene expression percentages calculated")
def count_gene_percentages(filt_mat: list, cell_data_dict: dict, gene_ids: list) -> dict:
    """Counts the percentage of cells each gene is expressed in, per-cluster. Returns a dict of each cell type (in
    order of our CELL_TYPES) and their corresponding expression percentage for each gene."""
    expression_counts = {}
    for cell_type in CELL_TYPES:
        # This will store the expressed and unexpressed counts for each of our genes. First entry is total count,
        # second is expressed count. Needs to be a dict to be unique in memory and a list to be mutable.
        expression_counts[cell_type] = {}
        for gene_id in gene_ids:
            expression_counts[cell_type][gene_id] = [0, 0]
    filt_mat_chunk_cell_path = "intermediates/filt_mat_chunk_cell.pkl"
    if os.path.exists(filt_mat_chunk_cell_path):
        with open(filt_mat_chunk_cell_path, 'rb') as filt_mat_by_cells_in:
            matrix_by_cells = pkl.load(filt_mat_by_cells_in)
    else:
        _, matrix_by_cells = aggregate_expression_level("cells", filt_mat)
        with open(filt_mat_chunk_cell_path, 'wb') as filt_mat_by_cells_out:
            pkl.dump(matrix_by_cells, filt_mat_by_cells_out)
    for cell_id, cell_sparse_mat in enumerate(matrix_by_cells):
        # Skips empty cells
        if cell_id not in cell_data_dict or not cell_sparse_mat:
            continue
        cell_type = cell_data_dict[cell_id][3]
        # Skips cells with no cell type membership
        if cell_type == "NA":
            continue
        # Builds a list of the genes expressed in this cell and translates their indexes to our 1-1999 standard.
        expressed_genes = []
        for gene_id, _ in cell_sparse_mat:
            expressed_genes.append(gene_id)
        # expressed_gene_idxs = swap_gene_indexing(expressed_genes, gene_key)
        for gene_id in gene_ids:
            expression_counts[cell_type][gene_id][0] += 1
            if gene_id in expressed_genes:
                expression_counts[cell_type][gene_id][1] += 1
    expression_percentages = {}
    for cell_type, expression_data in expression_counts.items():
        expression_percentages[cell_type] = {}
        for gene_id, gene_exp in expression_data.items():
            expression_percentage = gene_exp[1] / gene_exp[0]
            expression_percentages[cell_type][gene_id] = expression_percentage
    return expression_percentages


@log_time("Expression levels calculated")
def calculate_fold_changes(filt_mat: list, cell_data_dict: dict):
    """Calculates the relative fold change from background of each gene in each cluster. Returns a dict keyed to
    clusters with subdicts keyed to genes, with values of the fold change of that gene in that cluster."""
    filt_mat_chunk_gene_path = "intermediates/filt_mat_chunk_gene.pkl"
    cell_ids = list(set(filt_mat[1]))
    if os.path.exists(filt_mat_chunk_gene_path):
        with open(filt_mat_chunk_gene_path, 'rb') as filt_mat_by_genes_in:
            matrix_by_genes = pkl.load(filt_mat_by_genes_in)
    else:
        _, matrix_by_genes = aggregate_expression_level("genes", filt_mat)
        with open(filt_mat_chunk_gene_path, 'wb') as filt_mat_by_genes_out:
            pkl.dump(matrix_by_genes, filt_mat_by_genes_out)
    # Will hold all legitimate cell ids and their corresponding clusters.
    good_cells = Bidict({})
    # Will hold the fold changes of each gene in each cluster.
    fold_changes = {}
    for cell_id in cell_ids:
        # Catches and skips erroneous cells
        if cell_id not in cell_data_dict:
            continue
        cell_data = cell_data_dict[cell_id]
        cell_doublet = cell_data[2]
        cell_cluster = cell_data[3]
        # Initializes the clusters for our fold_changes dictionary.
        if cell_cluster not in fold_changes:
            fold_changes[cell_cluster] = {}
        # Skips badly annotated cells.
        if cell_doublet or cell_cluster == "NA":
            continue
        good_cells[cell_id] = cell_cluster
    # Overall cell count for background calculations.
    norm_count = len(good_cells.keys())
    # Tracks progress, mostly for debugging.
    log("Calculating fold changes...")
    processed_genes = 0
    for gene_id, gene_sparse_mat in enumerate(matrix_by_genes):
        # Skips empty genes.
        if not gene_sparse_mat:
            continue
        else:
            processed_genes += 1
        # Makes a dictionary of each expressing cell keyed to its expression level.
        exp_dict = {}
        cells = []
        exprs = []
        for cell, exp in gene_sparse_mat:
            exp_dict[cell] = exp
            cells.append(cell)
            exprs.append(exp)
        # This will be a dictionary of the gene's expression level as lists across each cell.
        cluster_exp = {}
        total_exp = []
        for cell, cluster in good_cells.items():
            if cluster not in cluster_exp:
                cluster_exp[cluster] = []
            # Adds pseudocount of 0.01.
            if cell not in exp_dict:
                exp = 0.01
            else:
                exp = exp_dict[cell]
            cluster_exp[cluster].append(exp)
            total_exp.append(exp)
        exp_mean = np.mean(total_exp)
        exp_stdev = np.std(total_exp)
        # Normalizes to zero mean and unit variance, then stores the mean of the normalized expression for each cluster.
        # Also builds a full count of expression and total count to get the other piece of the mean.
        cluster_means = {}
        for cluster, exprs in cluster_exp.items():
            normalized_exp = [(exp - exp_mean) / exp_stdev for exp in exprs]
            cluster_sum = sum(normalized_exp)
            cluster_count = len(normalized_exp)
            cluster_means[cluster] = (cluster_sum, cluster_count)
        for cluster, (cluster_sum, cluster_count) in cluster_means.items():
            # The normalized sum across all data will be zero since it's mean-centered.
            bg_sum = 0 - cluster_sum
            bg_count = norm_count - cluster_count
            bg_mean = bg_sum / bg_count
            cluster_mean = cluster_sum / cluster_count
            fold_change = cluster_mean / bg_mean
            fold_changes[cluster][gene_id] = fold_change
        if processed_genes % 20 == 0:
            percent_processed = processed_genes / 20
            log(f"{percent_processed}% done.")
    return fold_changes


def main() -> None:
    """Manager function.Details of each called function are in the corresponding docstrings."""
    double_filt_mat_path = "intermediates/double_filt_exp_mat.pkl"
    if not os.path.exists(double_filt_mat_path):
        exit("Filtered matrix not found. Run parsesourcefiles.py to generate.")
    else:
        log("Log initialized...")
        log("Loading data...")
        with open(double_filt_mat_path, 'rb') as filt_mat_in:
            # This matrix is sparse, with three corresponding lists of gene id, cell id, and expression level.
            filt_mat = pkl.load(filt_mat_in)
        with open("intermediates/cell_data_dict.pkl", 'rb') as cell_data_in:
            cell_data_dict = pkl.load(cell_data_in)
        with open("intermediates/gene_data_dict.pkl", 'rb') as gene_data_in:
            gene_data_dict = pkl.load(gene_data_in)
            # genes_by_name = invert_hash(gene_data_dict, array_vals=True, identifier=0).inverse
        log("Data Loaded.")
    expression_percentage_path = "intermediates/expression_percentages.pkl"
    if os.path.exists(expression_percentage_path):
        with open(expression_percentage_path, 'rb') as exp_perc_in:
            expression_percentages = pkl.load(exp_perc_in)
    else:
        gene_ids = list(set(filt_mat[0]))
        expression_percentages = count_gene_percentages(filt_mat, cell_data_dict, gene_ids)
        with open(expression_percentage_path, 'wb') as exp_perc_out:
            pkl.dump(expression_percentages, exp_perc_out)
    fold_change_path = "intermediates/fold_changes.pkl"
    if os.path.exists(fold_change_path):
        with open(fold_change_path, 'rb') as fc_in:
            fold_changes = pkl.load(fc_in)
    else:
        fold_changes = calculate_fold_changes(filt_mat, cell_data_dict)
        with open(fold_change_path, 'wb') as fc_out:
            pkl.dump(fold_changes, fc_out)


if __name__ == "__main__":
    main()
