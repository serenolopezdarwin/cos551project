"""
Given a filtered cell-by-gene expression sparse matrix, calculates the most enriched genes in each of our defined
cell types and plots them in percentage-by-log-fold-change dot plots.
"""
# Imports
from cos551 import *
import numpy as np
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.lines as lns
# noinspection PyUnresolvedReferences
from mpl_toolkits import mplot3d
import os
import pandas as pd
import pickle as pkl
import seaborn as sns
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
TRAJECTORIES = ["Endothelial trajectory", "Mesenchymal trajectory", "Neural tube and notochord trajectory",
                "Neural crest (melanocyte) trajectory 1", "Haematopoiesis trajectory", "Epithelial trajectory",
                "Epithelial trajectory", "Hepatocyte trajectory", "Neural crest (PNS glia) trajectory 2",
                "Neural crest (PNS neuron) trajectory 3", "Lens trajectory"]
# Marker genes originally identified in Cao et al 2019.
TEST_GENES = ["Col6a6", "Glis1", "Nr1h5", "Col9a1", "Ntng1", "Trp63", "Pth2r", "Fndc3c1", "Mybl1", "Tfap2d",
              "C130060K24Rik", "Dmbx1", "Mylk4", "Foxb1", "Npy", "Rab5a", "Il31ra", "Pax2", "Id4", "Emcn", "Lamc3",
              "Tspan8", "Mpz", "Ppp1r1c", "Cpa2", "Hbb-bh1", "Dlx6", "Eomes", "A1cf", "Metrnl", "Ms4a4a", "Gmnc",
              "Uts2b", "Myh6", "Gp1ba", "Tyr", "Cryba2", "Lcn2"]


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
    total_counts = {}
    for cell_type in CELL_TYPES:
        # This will store the expressed and unexpressed counts for each of our genes. First entry is total count,
        # second is expressed count. Needs to be a dict to be unique in memory and a list to be mutable.
        expression_counts[cell_type] = {}
        for gene_id in gene_ids:
            expression_counts[cell_type][gene_id] = 0
        total_counts[cell_type] = 0
    filt_mat_chunk_cell_path = "intermediates/filt_mat_chunk_cell.pkl"
    if os.path.exists(filt_mat_chunk_cell_path):
        with open(filt_mat_chunk_cell_path, 'rb') as filt_mat_by_cells_in:
            matrix_by_cells = pkl.load(filt_mat_by_cells_in)
    else:
        _, matrix_by_cells = aggregate_expression_level("cells", filt_mat)
        with open(filt_mat_chunk_cell_path, 'wb') as filt_mat_by_cells_out:
            pkl.dump(matrix_by_cells, filt_mat_by_cells_out)
    # Just for timing introspection of loop.
    log("Calculating gene expression percentages...")
    for cell_id, cell_sparse_mat in enumerate(matrix_by_cells):
        # Skips empty cells
        if cell_id not in cell_data_dict or not cell_sparse_mat:
            continue
        cell_type = cell_data_dict[cell_id][3]
        # Skips cells with no cell type membership
        if cell_type == "NA":
            continue
        # Tracks total number of cells in this cluster
        total_counts[cell_type] += 1
        # Builds a list of the genes expressed in this cell and translates their indexes to our 1-1999 standard.
        expressed_genes = []
        for gene_id, _ in cell_sparse_mat:
            expressed_genes.append(gene_id)
            expression_counts[cell_type][gene_id] += 1
        # expressed_gene_idxs = swap_gene_indexing(expressed_genes, gene_key)
        if cell_id % 20586 == 0:
            percent = cell_id // 20586
            log(f"{percent}% Done.")
    expression_percentages = {}
    for cell_type, expression_data in expression_counts.items():
        expression_percentages[cell_type] = {}
        total_cells = total_counts[cell_type]
        for gene_id, gene_exp in expression_data.items():
            expression_percentage = gene_exp / total_cells
            expression_percentages[cell_type][gene_id] = expression_percentage
    return expression_percentages


@log_time("Expression levels calculated")
def calculate_fold_changes(filt_mat: list, cell_data_dict: dict, gene_ids: list) -> dict:
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
    # Will hold the mean expression level of each cluster
    cluster_expression = {}
    cluster_mean_exp = {}
    # Initializes the clusters for our fold_changes dictionary.
    for cell_type in CELL_TYPES:
        fold_changes[cell_type] = {}
        cluster_expression[cell_type] = [0, 0]
        cluster_mean_exp[cell_type] = 0
    for cell_id in cell_ids:
        # Catches and skips erroneous cells
        if cell_id not in cell_data_dict:
            continue
        cell_data = cell_data_dict[cell_id]
        cell_doublet = cell_data[2]
        cell_cluster = cell_data[3]
        # Skips badly annotated cells.
        if cell_doublet or cell_cluster == "NA":
            continue
        good_cells[cell_id] = cell_cluster
        cell_exp = cell_data[0]
        cluster_expression[cell_cluster][0] += cell_exp
        cluster_expression[cell_cluster][1] += 1
    # Calculates mean expression level of cells in each cluster
    for cell_type, exp_data in cluster_expression.items():
        cluster_mean_exp[cell_type] = exp_data[0] / exp_data[1]
    # Normalizes all expression levels by the highest expression level.
    max_exp = max([exp for exp in cluster_mean_exp.values()])
    for cell_type, exp in cluster_mean_exp.items():
        cluster_mean_exp[cell_type] = exp / max_exp
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
            if cell not in exp_dict:
                exp = 0
            else:
                exp = exp_dict[cell]
            # Adds pseudocount of 0.01
            adjusted_exp = exp + 0.01
            cluster_exp[cluster].append(adjusted_exp)
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
            # Normalizes by average expression level of each cluster.
            cluster_exp = cluster_mean_exp[cluster]
            # We can get negatives so we want the absolute difference between these two.
            enrichment = (cluster_mean - bg_mean) / cluster_exp
            fold_changes[cluster][gene_id] = enrichment
        if processed_genes % 20 == 0:
            percent_processed = processed_genes // 20
            log(f"{percent_processed}% done.")
    # Normalizes all fold-changes to be positive.
    for gene_id in gene_ids:
        min_fc = min([fold_changes[cell_type][gene_id] for cell_type in CELL_TYPES])
        for cell_type in CELL_TYPES:
            fold_changes[cell_type][gene_id] -= min_fc
    return fold_changes


def generate_dot_plot(exp_perc: dict, fold_changes: dict, gene_ids: list, gene_data_dict: dict) -> list:
    """Based on our calculated fold changes and expression percentages, generates a dot plot of each cluster's most
    differentially expressed gene, sized by its expression percentage and colored by its adjusted expression level
    across all clusters. Returns a list of all marker genes we identify."""
    # Creates a path for our figures.
    if not os.path.exists("figures/"):
        os.mkdir("figures/")
    # Holds the genes with greatest fold-change for each cluster.
    max_fc = {}
    for cell_type in CELL_TYPES:
        max_fc[cell_type] = (0, 0, 0)
    for gene_id in gene_ids:
        fold_change_list = [fold_changes[cell_type][gene_id] for cell_type in CELL_TYPES]
        top_two = sorted(fold_change_list)[-2:]
        max_entry = top_two[1]
        second_entry = top_two[0]
        enrichment = max_entry / second_entry
        # Finds cluster with highest fold change and adds the gene id and two entries to that place in max_fc
        for cell_type in CELL_TYPES:
            if fold_changes[cell_type][gene_id] == max_entry:
                if enrichment > max_fc[cell_type][1] and max_entry > max_fc[cell_type][2]:
                    max_fc[cell_type] = (gene_id, enrichment, max_entry)
    # Makes a dataframe with columns of each representative gene and rows of each cell type.
    fc_dict = {}
    for cell_type in CELL_TYPES:
        rep_gene = max_fc[cell_type][0]
        fc_list = []
        ep_list = []
        for ct in CELL_TYPES:
            ep = exp_perc[ct][rep_gene]
            fc = fold_changes[ct][rep_gene]
            fc_list.append(fc * ep)
            ep_list.append(ep)
        min_fc = min([fc for fc in fc_list if fc])
        # Pseudocount of .1 to separate zero-values from very low values.
        adjusted_fc = [np.log10((fc / min_fc + .1)) + 1 for fc in fc_list]
        # Converts our two lists into a list of tuples.
        gene_data = [(fc, ep) for fc, ep in zip(adjusted_fc, ep_list)]
        fc_dict[cell_type] = gene_data
    # Gets the cluster name (and corresponding idx), gene name (and corresponding idx), fold change and exp perc.
    gene_idxs = []
    cluster_idxs = []
    full_fc = []
    full_ep = []
    gene_names = []
    for gene_idx, cell_type in enumerate(CELL_TYPES):
        gene_data = fc_dict[cell_type]
        rep_gene = max_fc[cell_type][0]
        gene_name = gene_data_dict[rep_gene][0]
        gene_names.append(gene_name)
        # We appended to these lists in the order of CELL_TYPES so they can be our indices.
        for idx, (fc, ep) in enumerate(gene_data):
            gene_idxs.append(gene_idx)
            # Flip the indices because we want to label clusters from top to bottom.
            cluster_idxs.append(37 - idx)
            full_fc.append(fc)
            full_ep.append(ep)
    point_frame = pd.DataFrame(data={"gene_idx": gene_idxs, "cluster_idx": cluster_idxs, "Scaled Expression": full_fc,
                                     "Expression Percentage": full_ep})
    # Sets plot size.
    mpl.rcParams['figure.figsize'] = 10, 10
    # Colors by fold change, sizes by expression percentage.
    dot_plot = sns.scatterplot(data=point_frame, x="gene_idx", y="cluster_idx", size="Expression Percentage",
                               hue="Scaled Expression", sizes=(1, 400), edgecolor="none")
    dot_plot.set(xlabel=None, ylabel=None)
    # Ensures we label every cluster and gene.
    seq_along = [n for n in range(0, len(CELL_TYPES))]
    plt.xticks(seq_along, gene_names, rotation=90)
    plt.yticks(seq_along, CELL_TYPES[::-1])
    # Puts the legend below the plot and adjusts margins to show it.
    plt.legend(bbox_to_anchor=(0.45, -.4), loc='lower center', borderaxespad=0., ncol=2)
    fig = dot_plot.get_figure()
    fig.savefig("figures/dot_plot.png", bbox_inches='tight')
    plt.close()
    return gene_names


def merge_colormaps(cm_names: list, num_colors: int) -> mpl.colors.ListedColormap:
    """Given any number of matplotlib colormaps, merges them and returns a colormap with num_colors entries from each of
    the input colormaps, from first to last. Iterates over the colormaps at once rather than sequentially."""
    colors = []
    cm_count = len(cm_names)
    # Makes a list of each colorant's color list.
    cms = [[c for c in plt.get_cmap(cm_name).colors] for cm_name in cm_names]
    i = 0
    while i < num_colors:
        cm_idx = i % cm_count
        try:
            color = cms[cm_idx][i // cm_count]
            colors.append(color)
            i += 1
        # If we run out of one colormap, simply move onto the next and extend the loop by one iteration.
        except IndexError:
            i += 1
            num_colors += 1
            pass
    colormap = col.ListedColormap(colors)
    return colormap


def generate_tsne_plot(cell_data_dict: dict):
    """Given a dictionary of cell cluster and tsne data, plots the cells in tsne space and colors by cluster."""
    # Makes a new qualitative colormap with colors for each of our cell types.
    tsne_colors = merge_colormaps(['Set1', 'Pastel1', 'Dark2', 'Set2', 'Pastel2'], len(CELL_TYPES)).colors
    tsne_1 = []
    tsne_2 = []
    cluster_labels = []
    cell_ages = []
    for cell_data in cell_data_dict.values():
        cell_doublet = cell_data[2]
        cell_cluster = cell_data[3]
        cell_age = cell_data[6]
        # Skips bad cells
        if cell_doublet or cell_cluster not in CELL_TYPES:
            continue
        cell_tsne = cell_data[1]
        tsne_1.append(cell_tsne[0])
        tsne_2.append(cell_tsne[1])
        cluster_labels.append(cell_cluster)
        cell_ages.append(cell_age)
    point_frame = pd.DataFrame(data={"tsne1": tsne_1, "tsne2": tsne_2, "Cluster": cluster_labels, "Age": cell_ages})
    # Sets plot size
    sns.set(rc={'figure.figsize': (20, 20)})
    # Colors by cluster membership and sets static small size. Takes forever to plot.
    tsne_plot = sns.scatterplot(data=point_frame, x="tsne1", y="tsne2", size=1, hue=cluster_labels, palette=tsne_colors,
                                edgecolor="none")
    fig = tsne_plot.get_figure()
    # Removes the "size" category from the legend.
    handles, labels = tsne_plot.get_legend_handles_labels()
    new_handles = [0] * len(CELL_TYPES)
    new_labels = [0] * len(CELL_TYPES)
    for handle, cluster in zip(handles[:-1], labels[:-1]):
        cluster_idx = CELL_TYPES.index(cluster)
        new_handles[cluster_idx] = handle
        # This is a PyCharm bug, there is no instance where cluster_idx will be a string here.
        # noinspection PyTypeChecker
        new_labels[cluster_idx] = f"{str(cluster_idx + 1)} - {cluster}"
    #  Puts the legend to the right of the plot and adjusts margins to show it. Also removes last entry (size).
    plt.legend(new_handles, new_labels, bbox_to_anchor=(1.4, 1), loc='upper right', borderaxespad=0, prop={'size': 20})
    # We don't need axes on a TSNE plot.
    plt.axis('off')
    fig.savefig("figures/tsne_plot.png", bbox_inches='tight')
    plt.close()
    # By-age tSNE plots. Take forever to plot.
    sns.set(font_scale=5)
    tsne_age_plots = sns.FacetGrid(data=point_frame, col="Age", hue="Cluster", palette=tsne_colors,
                                   col_order=['9.5', '10.5', '11.5', '12.5', '13.5'], height=15)
    # Removes the axis from every subplot
    for _, ax in tsne_age_plots.axes_dict.items():
        ax.axis('off')
    # Re-orders the plots in numerical age order.
    tsne_age_plots.map(sns.scatterplot, "tsne1", "tsne2", edgecolor="none", size=1)
    tsne_age_plots.figure.savefig("figures/tsne_age_plot.png")
    plt.close()


def export_legend(legend, filename="legend.png", exp1=-5, exp2=-5, exp3=5, exp4=5):
    """Given a matplotlib legend object, exports an image of the legend alone to the designated filename.
    NOT MY CODE: from https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib"""
    expand = [exp1, exp2, exp3, exp4]
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)
    plt.close()


def generate_umap_plot(cell_data_dict: dict):
    """"""

    trajectories = []
    umap = [[], [], []]
    for cell_data in cell_data_dict.values():
        traj = cell_data[4]
        # Optimized skipping of invalid trajectories.
        try:
            traj_idx = TRAJECTORIES.index(traj)
        except ValueError:
            continue
        # Adds umap coordinates to our umap lists.
        for idx in range(3):
            coord = cell_data[5][idx]
            umap[idx].append(coord)
        trajectories.append(traj_idx)

    for i in range(0, 180, 20):
        for j in range(0, 360, 45):
            fig = plt.figure(figsize=(100, 100))
            ax = plt.axes(projection='3d')
            ax.scatter(umap[0], umap[1], umap[2], c=trajectories, cmap='Set1')
            ax.view_init(i, j)
            # Coordinates once again don't really matter here.
            ax.set_facecolor('white')
            ax.axis('off')
            fig.savefig(f"figures/umap_test/umap_{str(i)}_{str(j)}.png", bbox_inches='tight')
            plt.close()
            print(f"Elev: {i} / 180, Azi: {j} / 360")
    # Makes and plots a legend that identifies our UMAP clusters.
    handle_colors = plt.get_cmap("Set1").colors[0:len(TRAJECTORIES)]
    handles = []
    labels = []
    for idx, color in enumerate(handle_colors):
        handle = lns.Line2D([0], [0], marker='o', color='w', label='Circle',
                            markerfacecolor=color, markersize=10)
        handles.append(handle)
        labels.append(f"{str(idx + 1)} - {TRAJECTORIES[idx]}")
    legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=True)
    export_legend(legend, "figures/umap_legend.png")


def compare_markers(markers: list, gene_ids: list, genes_by_name: dict):
    """Given a list of identified marker genes, compares them to the originally identified markers and tries to figure
    out why they escaped our classification."""
    re_genes = 0
    skipped_genes = 0
    filtered_genes = 0
    for idx, gene_name in enumerate(TEST_GENES):
        gene_id = genes_by_name[gene_name][0]
        cluster = CELL_TYPES[idx]
        marker = markers[idx]
        if gene_name in markers:
            log(f"{gene_name} re-identified for {cluster}.")
            re_genes += 1
        elif gene_id in gene_ids:
            log(f"{gene_name} not identified for {cluster}. Instead: {marker}")
            skipped_genes += 1
        else:
            log(f"{gene_name} filtered out for {cluster}. Instead: {marker}")
            filtered_genes += 1
    log(f"{str(re_genes)} genes re-identified.")
    log(f"{str(skipped_genes)} genes included but skipped.")
    log(f"{str(filtered_genes)} genes filtered out.")


def main() -> None:
    """Manager function.Details of each called function are in the corresponding docstrings."""
    # Using python two breaks our division, so this makes sure we're on the right version.
    if sys.version_info[0] < 3:
        raise VersionError("Must be using Python 3. Please update to Python 3 before running this script.")
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
            genes_by_name = invert_hash(gene_data_dict, array_vals=True, identifier=0).inverse
        log("Data Loaded.")
    gene_ids = list(set(filt_mat[0]))
    expression_percentage_path = "intermediates/expression_percentages.pkl"
    if os.path.exists(expression_percentage_path):
        with open(expression_percentage_path, 'rb') as exp_perc_in:
            expression_percentages = pkl.load(exp_perc_in)
    else:
        expression_percentages = count_gene_percentages(filt_mat, cell_data_dict, gene_ids)
        with open(expression_percentage_path, 'wb') as exp_perc_out:
            pkl.dump(expression_percentages, exp_perc_out)
    fold_change_path = "intermediates/fold_changes.pkl"
    if os.path.exists(fold_change_path):
        with open(fold_change_path, 'rb') as fc_in:
            fold_changes = pkl.load(fc_in)
    else:
        fold_changes = calculate_fold_changes(filt_mat, cell_data_dict, gene_ids)
        with open(fold_change_path, 'wb') as fc_out:
            pkl.dump(fold_changes, fc_out)
    # This script always makes plots.
    gene_names = generate_dot_plot(expression_percentages, fold_changes, gene_ids, gene_data_dict)
    generate_tsne_plot(cell_data_dict)
    compare_markers(gene_names, gene_ids, genes_by_name)


if __name__ == "__main__":
    main()
