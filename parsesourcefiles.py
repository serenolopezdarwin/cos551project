"""
Given cell and gene annotations as well as a sparse expression matrix, filters out cells that are detected doublets,
cells whose expression is more than 2 stdev from mean across all cells, and picks genes with top 2000 variance in
expression levels across cells.
"""
# Imports
from cos551 import *
import csv
import numpy as np
import os
import pickle as pkl
from time import perf_counter


def generate_exp_matrix(exp_data_path: str, sparse_mat_path: str):
    """Given a sparse matrix text file, creates a list of lists corresponding to each entry (row, col, val)"""
    start = perf_counter()
    with open(exp_data_path) as exp_data_in:
        # Skips header lines
        for _ in range(2):
            next(exp_data_in)
        # The three sublists correspond to gene ids, cell ids, and expression levels.
        sparse_mat = [[], [], []]
        for exp_data in exp_data_in:
            split_exp_data = exp_data.split()
            for idx in range(3):
                sparse_mat[idx].append(int(split_exp_data[idx]))
    with open(sparse_mat_path, 'wb') as sparse_mat_out:
        pkl.dump(sparse_mat, sparse_mat_out)
    log(f"Expression matrix generated in {str(perf_counter() - start)}")
    return sparse_mat


def aggregate_expression_level(by: str, exp_mat: list, sorted_exp_mat_path: str):
    """Given a dimension and an expression matrix, aggregates the expression level along each dimension, like doing the
    row or column sums of a non-sparse matrix. Returns the sums indexed by dimension ids (not indices, since dimensions
    are 1-indexed while indices are 0-indexed due to python's basic rules). Also returns the sorted zip object so that
    we avoid re-sorting in other functions, as it's computationally expensive for large lists."""
    start = perf_counter()
    if by == "genes":
        slice_idx = 0
        other_idx = 1
    elif by == "cells":
        slice_idx = 1
        other_idx = 0
    else:
        raise KeyError("Incorrect 'by' argument passed to aggregation function")
    # Avoids sorting bottleneck, and lets us verify our chunking method more easily.
    if os.path.exists(sorted_exp_mat_path):
        with open(sorted_exp_mat_path, 'rb') as sorted_exp_mat_in:
            exp_mat_sort = pkl.load(sorted_exp_mat_in)
    else:
        id_list = exp_mat[slice_idx]
        unprocessed_list = exp_mat[other_idx]
        exp_list = exp_mat[2]
        exp_mat_sort = sorted(zip(id_list, exp_list, unprocessed_list))
        with open(sorted_exp_mat_path, 'wb') as sorted_exp_mat_out:
            pkl.dump(exp_mat_sort, sorted_exp_mat_out)
    # The exp_mat_sort call finds the max of id_list without having to iterate over the list again. Adds one because
    # Anything we initialize with it is 1-indexed.
    max_idx = exp_mat_sort[-1][0] + 1
    # We also add an extra index, since the first index (0) will be empty due to cells and genes both being 1-indexed.
    # This will simplify later references to the generated data object downstream.
    aggregate_xp_list = [0] * max_idx
    aggregate_xp, aggregate_data = 0, []
    # Maximum index we can reference in the matrix
    exp_mat_max = len(exp_mat_sort) - 1
    chunked_exp_mat = [[] for _ in range(max_idx)]
    for idx, (id_num, exp, other_id) in enumerate(exp_mat_sort):
        aggregate_xp += exp
        aggregate_data.append((other_id, exp))
        # Checks if the next entry in the list is from a different id or doesn't exist.
        if idx == exp_mat_max or exp_mat_sort[idx+1][0] != id_num:
            # That being the case, dumps our current aggregated expression level and resets it.
            aggregate_xp_list[id_num] = aggregate_xp
            # Also dumps the list of other coord, exp tuples to the index of this id to be calculated later.
            chunked_exp_mat[id_num] = aggregate_data
            aggregate_xp, aggregate_data = 0, []
    log(f"Exp matrix chunked by {by} in {str(perf_counter() - start)}")
    return aggregate_xp_list, chunked_exp_mat


def process_cell_data(cell_annotation_path: str, exp_mat: list, cell_data_path: str, filtered_mat_path: str):
    """Reads a csv file of cell annotations and extracts various data fields for each cell. Also reads an input sparse
    expression matrix, sorts it by cell ids, and calculates each cell's aggregate expression level. Then flags cells for
    removal from the dataset based on their overall expression levels or doublet membership."""
    start = perf_counter()
    chunked_exp_mat_path = "intermediates/exp_mat_chunk_cells.pkl"
    # Avoids chunking bottleneck, and lets us verify our sorting method more easily.
    if os.path.exists(chunked_exp_mat_path):
        with open(chunked_exp_mat_path, 'rb') as chunked_exp_mat_in:
            exp_by_cells = pkl.load(chunked_exp_mat_in)
            cell_exp_levels = []
            for cell_data in exp_by_cells:
                cell_exp = sum([entry[1] for entry in cell_data])
                cell_exp_levels.append(cell_exp)
    else:
        sort_exp_mat_path = "intermediates/exp_mat_sort_cells.pkl"
        cell_exp_levels, exp_by_cells = aggregate_expression_level("cells", exp_mat, sort_exp_mat_path)
        with open(chunked_exp_mat_path, 'wb') as chunked_exp_mat_out:
            pkl.dump(exp_by_cells, chunked_exp_mat_out)
    cell_data_dict = {}
    with open(cell_annotation_path, 'rt') as cell_annotations_in:
        cell_reader = csv.reader(cell_annotations_in)
        # Skips header
        next(cell_reader)
        row_num = 0
        for cell_row in cell_reader:
            # Iterate first because cells are 1-indexed
            row_num += 1
            exp_level = cell_exp_levels[row_num]
            tsne = [float(cell_row[idx]) if cell_row[idx] != "NA" else "NA" for idx in [13, 14]]
            # Weirdly, casting these strings to bools doesn't work. So we just do a manual check...lol
            if cell_row[19] == "TRUE":
                doublet = True
            else:
                doublet = False
            cluster = cell_row[22]
            traj = cell_row[23]
            umap = [float(cell_row[idx]) if cell_row[idx] != "NA" else "NA" for idx in [24, 25, 26]]
            cell_data_dict[row_num] = [exp_level, tsne, doublet, cluster, traj, umap]
    # Don't count first entry, as it is an empty initialization entry due to our cells being 1-indexed
    exp_mean = np.mean(cell_exp_levels[1:])
    exp_stdev = np.std(cell_exp_levels[1:])
    exp_mat_filtered = [[], [], []]
    discarded_cells, empty_cells, processed_cells = 0, 0, 0
    for cell_id, exp_data in enumerate(exp_by_cells):
        if cell_id not in cell_data_dict:
            empty_cells += 1
            continue
        cell_data = cell_data_dict[cell_id]
        cell_exp = cell_data[0]
        cell_doublet = cell_data[2]
        # Checks if cell aggregate expression level is two standard deviations away from the mean. If it is, the
        # cell's "doublet" bool is changed to True, so that we can know it was discarded in later analysis.
        if cell_doublet or (np.abs(cell_exp-exp_mean)/exp_stdev) > 2:
            cell_data_dict[cell_id][2] = True
            discarded_cells += 1
        elif not cell_exp:
            cell_data_dict[cell_id][2] = True
            empty_cells += 1
        else:
            # Reformats our by-cell sparse matrix to the more generalized sparse matrix list-of-lists form.
            for gene_id, exp in exp_data:
                exp_mat_filtered[0].append(gene_id)
                exp_mat_filtered[1].append(cell_id)
                exp_mat_filtered[2].append(exp)
            processed_cells += 1
    log(f"Discarded cell count: {str(discarded_cells)}")
    log(f"Empty cell count: {str(empty_cells)}")
    log(f"Processed cell count: {str(processed_cells)}")
    with open(filtered_mat_path, 'wb') as filtered_mat_out:
        pkl.dump(exp_mat_filtered, filtered_mat_out)
    with open(cell_data_path, 'wb') as cell_data_out:
        pkl.dump(cell_data_dict, cell_data_out)
    log(f"Cells processed in {str(perf_counter() - start)}")
    return cell_data_dict, exp_mat_filtered


def generate_gene_dict(gene_annotation_path: str, filt_exp_mat: list, gene_data_path: str, double_filt_mat_path: str):
    """Given a file of gene annotations and an expression matrix with bad cells filtered out, stores relevant gene data
    in a pickled dictionary and calculates the variance across all good cells of reads in each gene, and picks the top
    2000 genes by variance, returning a matrix of just these genes (across good cells)."""
    start = perf_counter()
    sort_exp_mat_path = "intermediates/exp_mat_sort_genes.pkl"
    gene_exp_levels, exp_by_genes = aggregate_expression_level("genes", filt_exp_mat, sort_exp_mat_path)
    gene_data_dict = {}
    with open(gene_annotation_path, 'rt') as gene_annotations_in:
        gene_reader = csv.reader(gene_annotations_in)
        # Skips header
        next(gene_reader)
        row_num = 0
        for gene_row in gene_reader:
            # Iterate first because genes are 1-indexed
            row_num += 1
            gene_name = gene_row[2]
            gene_exp = gene_exp_levels[row_num]
            # We need to append additional info to this later, so we'll initialize the entry as a list.
            gene_data_dict[row_num] = [gene_name, gene_exp]
    with open(gene_data_path, 'wb') as gene_data_out:
        pkl.dump(gene_data_dict, gene_data_out)
    variance_list = []
    for gene_id, gene_data in enumerate(exp_by_genes):
        # Skips genes with no data.
        if not gene_data:
            continue
        # Variance of expression levels across expressing cells for this gene.
        gene_variance = np.var([entry[1] for entry in gene_data])
        variance_list.append(gene_variance)
    # Anything with variance more than or equal to this value is analyzed as one of the top 2000 most variant genes.
    variance_cutoff = sorted(variance_list[-2000])
    exp_mat_filtered = [[], [], []]
    for gene_id, gene_data in enumerate(exp_by_genes):
        # Skips genes with no data.
        if not gene_data:
            continue
        # Variance of expression levels across expressing cells for this gene.
        gene_variance = np.var([entry[1] for entry in gene_data])
        if gene_variance < variance_cutoff:
            gene_data_dict[gene_id].append(False)
            continue
        # For the top 2000 most variant gene, we write their expression level into our familiar sparse matrix format.
        else:
            gene_data_dict[gene_id].append(True)
            for (cell_id, exp) in gene_data:
                exp_mat_filtered[0].append(gene_id)
                exp_mat_filtered[1].append(cell_id)
                exp_mat_filtered[2].append(exp)
    with open(double_filt_mat_path, 'wb') as filt_mat_out:
        pkl.dump(exp_mat_filtered, filt_mat_out)
    log(f"Genes processed in {str(perf_counter() - start)}")
    return gene_data_dict, exp_mat_filtered


def main():
    """Manager function. Checks if a variety of files exist, and if they don't, generates them. Details of each called
    function are in the corresponding docstrings."""
    # Creates any missing data directories for us.
    log("Log initialized...")
    data_directories = ["sourcefiles", "intermediates"]
    for directory in data_directories:
        if not os.path.isdir(directory):
            os.mkdir(directory)
    source_files = ["gene_count.txt", "cell_annotate.csv", "gene_annotate.csv"]
    source_paths = [f"sourcefiles/{file}" for file in source_files]
    # Checks that we have all of the necessary sourcefiles in the right relative paths.
    for source_path in source_paths:
        if not os.path.exists(source_path):
            exit(f"Sourcefile {source_path} missing. Please provide the appropriate file.")
    sparse_mat_path = "intermediates/sparse_mat.pkl"
    if os.path.exists(sparse_mat_path):
        with open(sparse_mat_path, 'rb') as sparse_mat_in:
            exp_matrix = pkl.load(sparse_mat_in)
    else:
        exp_matrix = generate_exp_matrix(source_paths[0], sparse_mat_path)
    cell_data_path = "intermediates/cell_data_dict.pkl"
    cell_filt_mat_path = "intermediates/cell_filt_exp_mat.pkl"
    if os.path.exists(cell_data_path) and os.path.exists(cell_filt_mat_path):
        with open(cell_data_path, 'rb') as cell_data_in:
            cell_data_dict = pkl.load(cell_data_in)
        with open(cell_filt_mat_path, 'rb') as filtered_mat_in:
            cell_filt_exp_mat = pkl.load(filtered_mat_in)
    else:
        cell_data_dict, cell_filt_exp_mat = process_cell_data(source_paths[1], exp_matrix, cell_data_path,
                                                              cell_filt_mat_path)
    gene_data_path = "intermediates/gene_data_dict.pkl"
    double_filt_mat_path = "intermediates/double_filt_exp_mat.pkl"
    if os.path.exists(gene_data_path):
        with open(gene_data_path, 'rb') as gene_data_in:
            gene_data_dict = pkl.load(gene_data_in)
        with open(double_filt_mat_path, 'rb') as double_filt_mat_in:
            double_filt_mat = pkl.load(double_filt_mat_in)
    else:
        gene_data_dict, double_filt_mat = generate_gene_dict(source_paths[2], cell_filt_exp_mat, gene_data_path,
                                                             double_filt_mat_path)
    # Checks to see if our script executed successfully
    if cell_data_dict and gene_data_dict and double_filt_mat:
        print("All Done!")


if __name__ == "__main__":
    main()
