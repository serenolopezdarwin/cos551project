"""

"""
# Imports
from cos551 import *
import os
import pickle as pkl


def analyze_cell_data_dict(cell_data_dict: dict) -> [int, list]:
    """Given a dictionary of cell data, returns the cell count and a list of all cell types."""
    cell_count = len(cell_data_dict.keys())
    cell_types = []
    for cell_data in cell_data_dict.items():
        cell_type = cell_data[3]
        if cell_type not in cell_types:
            cell_types.append(cell_type)
    return cell_count, cell_types


def aggregate_matrix_by_cell_type(matrix: list, cell_data_dict: dict) -> dict:
    """Given a sparse matrix, splits it into a separate sparse matrix for each cell type. Stores the separated matrix
    in a dictionary keyed to cell types."""
    matrix_by_cell_type = {}
    for

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