"""

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


def count_gene_percentages(filt_mat: list, cell_data_dict: dict) -> dict:
    """Counts the percentage of cells each gene is expressed in, per-cluster. Returns a list of each cell type (in
    order of our CELL_TYPES) and their corresponding expression percentage for each gene."""
    start = perf_counter()
    for cell_type in CELL_TYPES:

    _, matrix_by_cells = aggregate_expression_level("cells", filt_mat)
    with open("intermediates/filt_mat_chunk_cell.pkl", 'wb')
    for cell_id, cell_sparse_mat in enumerate(matrix_by_cells):
        if cell_id not in cell_data_dict:
            continue
        cell_type = cell_data_dict[cell_id][3]
        if cell_type == "NA":
            continue


    log(f"Expression matrix separated by cell type in {str(perf_counter() - start)}")
    return matrix_by_cell_type



def main() -> None:
    """Manager function.Details of each called function are in the corresponding docstrings."""
    double_filt_mat_path = "intermediates/double_filt_exp_mat.pkl"
    if not os.path.exists(double_filt_mat_path):
        exit("Filtered matrix not found. Run parsesourcefiles.py to generate.")
    else:
        with open(double_filt_mat_path, 'rb') as filt_mat_in:
            # This matrix is sparse, with three corresponding lists of gene id, cell id, and expression level.
            filt_mat = pkl.load(filt_mat_in)
        with open("intermediates/cell_data_dict.pkl", 'rb') as cell_data_in:
            cell_data_dict = pkl.load(cell_data_in)
        with open("intermediates/gene_data_dict.pkl", 'rb') as gene_data_in:
            gene_data_dict = pkl.load(gene_data_in)



if __name__ == "__main__":
    main()