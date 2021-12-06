"""
Given a filtered cell-by-gene expression sparse matrix, calculates the most enriched genes in each of our defined
cell types and plots them in percentage-by-log-fold-change dot plots.
"""
# Imports
from cos551 import *
from time import perf_counter
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


def count_gene_percentages(filt_mat: list, cell_data_dict: dict, gene_ids: list) -> dict:
    """Counts the percentage of cells each gene is expressed in, per-cluster. Returns a dict of each cell type (in
    order of our CELL_TYPES) and their corresponding expression percentage for each gene."""
    start = perf_counter()
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
    log(f"Expression percentages calculated in {str(perf_counter() - start)}")
    return expression_percentages


def main() -> None:
    """Manager function.Details of each called function are in the corresponding docstrings."""
    double_filt_mat_path = "intermediates/double_filt_exp_mat.pkl"
    if not os.path.exists(double_filt_mat_path):
        exit("Filtered matrix not found. Run parsesourcefiles.py to generate.")
    else:
        log("Log initialized...")
        with open(double_filt_mat_path, 'rb') as filt_mat_in:
            # This matrix is sparse, with three corresponding lists of gene id, cell id, and expression level.
            filt_mat = pkl.load(filt_mat_in)
        with open("intermediates/cell_data_dict.pkl", 'rb') as cell_data_in:
            cell_data_dict = pkl.load(cell_data_in)
        with open("intermediates/gene_data_dict.pkl", 'rb') as gene_data_in:
            gene_data_dict = pkl.load(gene_data_in)
    expression_percentage_path = "intermediates/expression_percentages.pkl"
    if os.path.exists(expression_percentage_path):
        with open(expression_percentage_path, 'rb') as exp_perc_in:
            expression_percentages = pkl.load(exp_perc_in)
    else:
        gene_ids = list(set(filt_mat[0]))
        expression_percentages = count_gene_percentages(filt_mat, cell_data_dict, gene_ids)
        with open(expression_percentage_path, 'wb') as exp_perc_out:
            pkl.dump(expression_percentages, exp_perc_out)


if __name__ == "__main__":
    main()
