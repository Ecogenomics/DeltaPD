################################################################################
#                                                                              #
#  This file is part of DeltaPD.                                               #
#                                                                              #
#  DeltaPD is free software: you can redistribute it and/or modify             #
#  it under the terms of the GNU Affero General Public License as published by #
#  the Free Software Foundation, either version 3 of the License, or           #
#  (at your option) any later version.                                         #
#                                                                              #
#  DeltaPD is distributed in the hope that it will be useful,                  #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of              #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#  GNU Affero General Public License for more details.                         #
#                                                                              #
#  You should have received a copy of the GNU Affero General Public License    #
#  along with DeltaPD.  If not, see <https://www.gnu.org/licenses/>.           #
#                                                                              #
################################################################################


import os
from collections import defaultdict

import dendropy
import ete3
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from phylodm.pdm import PDM
from tqdm import tqdm

from deltapd.config import desc_pbar
from deltapd.model import set_xy_labels


def split_qry_gid(qry_gid, qry_sep):
    """Splits the accession of the query genome (GID~GENE)"""
    return qry_gid.split(qry_sep)


def split_qry_label(qry_gid, qry_sep):
    """Splits the accession of the query genome (GID~GENE)"""
    return qry_gid.split(qry_sep)


def get_lca(tree_str, gid_to_tax):
    tree = dendropy.Tree.get_from_string(tree_str, schema='newick', preserve_underscores=True)

    # Index the tree with descendant taxonomies.
    for leaf_node in tree.leaf_node_iter():
        cur_tax = gid_to_tax[leaf_node.taxon.label]
        cur_parent = leaf_node
        while cur_parent is not None:
            if not hasattr(cur_parent, 'desc_tax'):
                cur_parent.desc_tax = [set() for _ in range(len(cur_tax))]
            if not hasattr(cur_parent, 'desc_taxa'):
                cur_parent.desc_taxa = set()
            for i in range(len(cur_tax)):
                cur_parent.desc_tax[i].add(cur_tax[i])
            cur_parent.desc_taxa.add(leaf_node.taxon.label)
            cur_parent = cur_parent.parent_node

    # Find the internal nodes that represent the ranks.
    d_ranks = defaultdict(list)
    for i in range(1, 7):
        done = set()
        for leaf_node in tree.leaf_node_iter():
            if leaf_node.taxon.label in done:
                continue

            # Traverse up the tree until the current rank is no longer monophyletic.
            cur_parent = leaf_node.parent_node
            last_node = leaf_node
            while cur_parent is not None and len(cur_parent.desc_tax[i]) == 1:
                last_node = cur_parent
                cur_parent = cur_parent.parent_node

            # The last node was monophyletic.
            d_ranks[list(last_node.desc_tax[i])[0]].append(last_node)
            done = done.union(last_node.desc_taxa)

    # Convert the nodes to the LCA
    out = defaultdict(set)
    for rank_name, node_list in d_ranks.items():
        for node in node_list:
            if node.is_leaf():
                out[rank_name].add(node.taxon.label)
            else:
                left = node.child_nodes()[0]
                right = node.child_nodes()[1]
                out[rank_name].add((left.leaf_nodes()[0].taxon.label,
                                    right.leaf_nodes()[0].taxon.label))
    return out


def generate_colour_palette(tax_list, metadata):
    sub_ranks = metadata.sub_ranks()

    # Unique ranks at level
    d_ranks_at_level = defaultdict(set)
    for cur_tax in tax_list:
        for i in range(7):
            d_ranks_at_level[i].add(cur_tax[i])

    # Determine the highest rank that differs.
    highest_rank = 0
    for rank_level, rank_set in d_ranks_at_level.items():
        if len(rank_set) > 1:
            highest_rank = rank_level
            break

    # Starting at the highest rank which differs, determine the immediate children of each rank.
    d_top_down_ranks = {x: [set() for _ in range(7)] for x in d_ranks_at_level[highest_rank]}
    [d_top_down_ranks[x][highest_rank].add(x) for x in d_top_down_ranks]
    for parent_rank, parent_list in d_top_down_ranks.items():
        for i in range(highest_rank, 6):
            previous_gids = parent_list[i]
            for previous_gid in previous_gids:
                for sub_rank in sub_ranks[previous_gid]:
                    parent_list[i + 1].add(sub_rank)

    # Find out how many sub ranks are under each.
    d_top_rank_cnt = defaultdict(lambda: 0)
    for parent_rank, parent_list in d_top_down_ranks.items():
        for i in range(highest_rank, 7):
            d_top_rank_cnt[parent_rank] += len(parent_list[i])

    # Generate a colour pallete for each of the top ranks, big enough to have a colour
    # for each of the descendants.
    d_top_rank_palette = dict()
    for i, (top_rank, top_rank_cnt) in enumerate(d_top_rank_cnt.items()):
        cur_start = (3 / len(d_top_rank_cnt)) * i
        d_top_rank_palette[top_rank] = sns.cubehelix_palette(top_rank_cnt, start=cur_start, rot=0, dark=0.6,
                                                             reverse=True)

    # Assign each of the desc ranks these colours
    d_out = dict()
    for parent_rank, parent_list in d_top_down_ranks.items():
        n_used = 0
        for i in range(highest_rank, 7):
            for cur_rank in parent_list[i]:
                d_out[cur_rank] = rgb_to_hex(d_top_rank_palette[parent_rank][n_used])
                n_used += 1

    return d_out


def rgb_to_hex(rgb):
    r, g, b = [int(x * 255) for x in rgb]
    return '#%02x%02x%02x' % (r, g, b)


def get_pdms(path_r_tree, path_q_tree, cpus):
    """Creates a Phylogenetic Distance Matrix (PDM) for both the reference,
    and query trees.

    Parameters
    ----------
    path_r_tree : str
        The path to the reference tree in Newick format.
    path_q_tree : str
        The path to the query tree in Newick format.
    cpus : int
        The number of CPUs to use when loading the PDM.

    Returns
    -------
    PDM
        The PDM object for the reference tree.
    np.ndarray
        The symmetric distance matrix for the reference tree.
    np.ndarray
        Each of the positional labels for the reference tree.
    PDM
        The PDM object for the query tree.
    np.ndarray
        The symmetric distance matrix for the query tree.
    np.ndarray
        Each of the positional labels for the query tree.
    """
    with tqdm(total=4, bar_format=desc_pbar) as p_bar:
        p_bar.set_description_str('Generating reference PDM from Newick file.')
        r_pdm = PDM.get_from_newick_file(path_r_tree, method='pd', cpus=cpus)
        p_bar.update()

        p_bar.set_description_str('Converting reference PDM to symmetric matrix.')
        r_labels, r_mat = r_pdm.as_matrix(normalised=True)
        r_labels = np.array(r_labels)
        p_bar.update()

        p_bar.set_description_str('Generating query PDM from Newick file.')
        q_pdm = PDM.get_from_newick_file(path_q_tree, method='pd', cpus=cpus)
        p_bar.update()

        p_bar.set_description_str('Converting query PDM to symmetric matrix.')
        q_labels, q_mat = q_pdm.as_matrix(normalised=True)
        q_labels = np.array(q_labels)
        p_bar.update()

        return r_pdm, r_labels, r_mat, q_pdm, q_labels, q_mat


def get_mapping(r_labels, q_labels, get_rep_label_dict, qry_sep):
    """Calculates the q -> r (1:1) and r -> q (1:n) taxon mapping.

    Parameters
    ----------
    r_labels : np.ndarray
        The array containing the reference labels.
    q_labels : np.ndarray
        The array containing the query labels.
    get_rep_label_dict : Dict[str, str]
        The genome accession to the representative accession mapping.
    qry_sep : str
        The separator between the query taxon accession and gene id.

    Returns
    -------
    qry_idx_to_rep_idx : np.ndarray
        For a given query taxon, return the index of the representative taxon.
    rep_idx_to_qry_idxs : Dict[int, frozenset]
        For a given rep taxon, return the set of query taxa.
    """

    # Setup output variables.
    qry_idx_to_rep_idx = np.zeros(q_labels.size, dtype=np.int64)
    rep_idx_to_qry_idxs = dict()

    # Pre-index the labels.
    r_label_to_idx = {v: i for i, v in enumerate(r_labels)}

    # Populate the values.
    for qry_idx, qry_label in enumerate(q_labels):
        rep_label = get_rep_label_dict[split_qry_label(qry_label, qry_sep)[0]]
        rep_idx = r_label_to_idx[rep_label]
        qry_idx_to_rep_idx[qry_idx] = rep_idx
        if rep_idx not in rep_idx_to_qry_idxs:
            rep_idx_to_qry_idxs[rep_idx] = set()
        rep_idx_to_qry_idxs[rep_idx].add(qry_idx)

    # Make the values immutable.
    rep_idx_to_qry_idxs = {k: frozenset(v) for k, v in rep_idx_to_qry_idxs.items()}
    return qry_idx_to_rep_idx, rep_idx_to_qry_idxs


def get_knn(r_idx_to_q_idxs, r_mat, q_mat, k):
    """Return the k nearest rep taxa for each representative taxon.

    Returns
    -------
    rep_idx_to_rep_idxs : Dict[int, frozenset]
        For a given rep taxon index, return the set of rep taxa indexes that are it's knns.
        e.g. `rep_idx_to_rep_idxs[0] = {0, 1, 55}`
    r_to_q : np.ndarray
        For a given rep index, return a mask of query indexes that are represented by it.
        e.g. `r_to_q[0, :] = [T, T, F, F, T]`
    """

    # Define the output variables.
    rep_idx_to_rep_idxs = dict()
    r_to_q = np.zeros((r_mat.shape[0], q_mat.shape[0]), dtype=np.bool)

    # Sort the reference matrix.
    r_mat_sort = np.argsort(r_mat, axis=1)

    # Iterate over each of the rep idxs that represent at least one query genome.
    for rep_idx in tqdm(r_idx_to_q_idxs, smoothing=0.1):

        # Iterate over each of the closest -> furthest rep idxs.
        rep_nn_set, qry_nn_set = set(), set()
        for rep_nn_idx in r_mat_sort[rep_idx]:

            # Stop processing once those k neighbours have been found.
            if len(rep_nn_set) >= k:
                break

            # Only interested in those that rep a query.
            if rep_nn_idx in r_idx_to_q_idxs:
                rep_nn_set.add(int(rep_nn_idx))

                # Also record all of the query idxs that this nn reps.
                qry_nn_set = qry_nn_set.union(r_idx_to_q_idxs[rep_nn_idx])

        # Store these indices.
        rep_idx_to_rep_idxs[rep_idx] = frozenset(rep_nn_set)
        for qry_idx in qry_nn_set:
            r_to_q[rep_idx, qry_idx] = True

    return rep_idx_to_rep_idxs, r_to_q


def get_unq_rep_idx_mapping(rep_idx_to_rep_idxs):
    """Given the dictionary that gives the rep indexes for each knn, return
    only the rep idxs that give a unique mapping to their knn rep idxs.

    Parameters
    ----------
    rep_idx_to_rep_idxs : Dict[int, frozenset]
        Given a representative index, return the rep idxs that are its knns.

    Returns
    -------
    np.ndarray
        An integer array containing only the rep idxs that give a unique mapping.
    """
    out = set()
    seen = set()
    for rep_idx, rep_idx_set in rep_idx_to_rep_idxs.items():
        if rep_idx_set in seen:
            continue
        out.add(rep_idx)
        seen.add(rep_idx_set)
    return np.array(list(out))


def get_largest_jk_size(r_idx_to_k_q_idxs_mask):
    """Determine the largest number of query taxa that are represented
    by any representative taxon.

    Parameters
    ----------
    r_idx_to_k_q_idxs_mask : np.ndarray
        Given a rep index, return a mask of qry indexes represented.

    Returns
    -------
    int
        The largest number of query genomes represented in all sets.
    """
    return int(np.max(np.sum(r_idx_to_k_q_idxs_mask, axis=1)))


def get_qry_taxon_groups(q_labels, qry_sep):
    """Group all of the query taxa based on the common accessions.

    Parameters
    ----------
    q_labels : np.ndarray
        All of the query labels.
    qry_sep : str
        The separator between the query taxon accession and gene id.

    Returns
    -------
    np.ndarray:
        An array of integers containing the groupings.
        e.g.: A~123, A~111, B~123, C~55 = [0, 0, 1, 2]
    """
    all_labels = dict()
    i = 0
    for val in q_labels:
        start = val.split(qry_sep)[0]
        if start not in all_labels:
            all_labels[start] = i
            i += 1
    return np.array([all_labels[x.split(qry_sep)[0]] for x in q_labels])


def generate_scatterplot(base_model, q_label_outliers, jk_models, out_dir, rep_label, jk_sets, base_xy):
    # Create a plot, containing all of the base x and y points.
    # For any X and Y points that go to an outlier, put those in a different legend.

    # For each of the outliers within this, plot the scatterplot without these points as subplot
    x_in, y_in = list(), list()
    x_out, y_out = list(), list()
    for cur_x, cur_y, cur_label in zip(base_xy['x'], base_xy['y'], base_xy['labels']):
        if any([x in q_label_outliers for x in cur_label]):
            x_out.append(cur_x)
            y_out.append(cur_y)
        else:
            x_in.append(cur_x)
            y_in.append(cur_y)

    # Generate a grid to display each of the jackknife groups.
    grid_shape = int(np.ceil(np.sqrt(len(jk_models) + 1)))

    alpha = 0.6

    # Create subplots
    fig, subplots = plt.subplots(grid_shape, grid_shape, figsize=(30, 30), sharex=True, sharey=True)
    ordered_jk_sets = sorted(jk_models.items(), key=lambda x: -abs(x[1]['influence']))
    # Add this to account for the base model being zero
    ordered_jk_sets.insert(0, None)
    next_plot = 0
    x_space = np.linspace(0, max(base_xy['x'].max(), base_xy['y'].max()), 10)
    title_size = 8
    for row_idx in range(grid_shape):
        for col_idx in range(grid_shape):
            ax = subplots[row_idx, col_idx]

            # Set the remaining plots invisible.
            if next_plot >= len(ordered_jk_sets):
                ax.set_visible(False)
                continue

            # First plot is always the base model.
            if next_plot == 0:
                ax.scatter(x_in, y_in, label='Inliers', marker='o', alpha=alpha, edgecolors='white')
                ax.scatter(x_out, y_out, label='Outliers', marker='x', alpha=alpha, edgecolors='white')

                # Plot the fit model.
                fit = np.poly1d([base_model['gradient'], 0])
                ax.plot(x_space, fit(x_space), 'm--')
                title_objs = ['Base model',
                              f'[y={base_model["gradient"]:.4f}x]',
                              f'[MSE={base_model["mse"]:.4f}]',
                              f'[R2={base_model["r2"]:.4f}]']


            else:
                cur_q_label, cur_jk_model = ordered_jk_sets[next_plot]

                ax.scatter(jk_sets[cur_q_label]['x'],
                           jk_sets[cur_q_label]['y'], marker='o',
                           alpha=alpha, edgecolors='white')

                title_objs = [cur_q_label,
                              f'[y={cur_jk_model["gradient"]:.4f}x]',
                              f'[MSE={cur_jk_model["mse"]:.4f}]',
                              f'[R2={cur_jk_model["r2"]:.4f}]',
                              f'[Inf={cur_jk_model["influence"]:.4f}]']

                # Plot the fit model.
                fit = np.poly1d([cur_jk_model["gradient"], 0])
                ax.plot(x_space, fit(x_space), 'm--')

            # Common
            ax.set_title('  '.join(title_objs), size=title_size)
            ax.plot(x_space, x_space, 'k--', label='y=x', alpha=0.2)
            ax.set_aspect('equal')
            ax.grid(True)

            next_plot += 1

    if True:
        for outlier_label in q_label_outliers:
            cur_dir = os.path.join(out_dir, outlier_label)
            if not os.path.isdir(cur_dir):
                os.makedirs(cur_dir)
            fig_path = f'{os.path.join(cur_dir, rep_label)}.png'
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

    plt.close(fig)

    return


def plt_tree(path_r_tree, path_q_tree, cur_q_labels, cur_r_labels, metadata,
             q_label_outliers, out_dir, q_labels, r_labels, base_model, jk_models,
             rep_idx, q_idx_to_r_idx, rep_label, ete3_scale, qry_sep, inf_thresh):
    # Determine what nodes can be used to annotate each rank.
    r_tree_str = retain_taxa_in_tree(path_r_tree, cur_r_labels)
    q_tree_str = retain_taxa_in_tree(path_q_tree, cur_q_labels)

    # Get the taxonomy of each taxon.
    r_gid_to_tax = {x: metadata.get_gid_taxonomy(x) for x in cur_r_labels}
    q_gid_to_tax = {x: metadata.get_gid_taxonomy(split_qry_gid(x, qry_sep)[0]) for x in cur_q_labels}

    # Generate a colour palette for the top-most rank that differs.
    unique_ranks = set([tuple(x) for x in sorted({**r_gid_to_tax, **q_gid_to_tax}.values())])
    palette = generate_colour_palette(unique_ranks, metadata)

    draw_ete_reference_tree(unique_ranks,
                            tree_str=q_tree_str,
                            d_gid_to_tax=q_gid_to_tax,
                            palette=palette,
                            metadata=metadata,
                            jk_models=jk_models,
                            out_dir=out_dir,
                            rep_label=rep_label,
                            ete3_scale=ete3_scale,
                            inf_thresh=inf_thresh,
                            q_label_outliers=q_label_outliers)

    draw_ete_reference_tree(unique_ranks,
                            tree_str=r_tree_str,
                            d_gid_to_tax=r_gid_to_tax,
                            palette=palette,
                            metadata=metadata,
                            jk_models=None,
                            out_dir=out_dir,
                            rep_label=rep_label,
                            ete3_scale=ete3_scale,
                            inf_thresh=inf_thresh,
                            q_label_outliers=q_label_outliers)
    return


def count_taxonomy_ranks(unique_ranks):
    """For a given collection of tuples, containing the full taxonomy.
    Report the number of unique ranks at each level, and the highest rank
    at which they differ."""
    d_ranks_at_level = defaultdict(set)
    for cur_tax in unique_ranks:
        for i in range(7):
            d_ranks_at_level[i].add(cur_tax[i])

    # Determine the highest rank that differs.
    highest_rank = 0
    for rank_level, rank_set in d_ranks_at_level.items():
        if len(rank_set) > 1:
            highest_rank = rank_level
            break
    return d_ranks_at_level, highest_rank


def ete3_draw_influence(tree, jk_models, metadata, at_col, n_inf_max):
    """Adds an influence column to each leaf node."""

    # Generate colour palettes for the influence/delta function.
    n_heatmap = 100
    heatmap = [rgb_to_hex(x) for x in sns.color_palette('coolwarm', n_heatmap + 1)]

    all_contig_lens = [metadata.get_contig_len(x) for x in jk_models]
    contig_min = np.min(all_contig_lens)
    contig_max = np.max(all_contig_lens)

    def contig_len_to_idx(contig_len):
        return abs(int(((contig_len - contig_min) / (contig_max - contig_min)) * 100) - 100)

    for leaf in tree.traverse():
        if leaf.is_leaf():
            cur_inf = jk_models[leaf.name]['influence']
            cur_delta = jk_models[leaf.name]['delta']

            append_str = f'i/d [{cur_inf:.4f}/{cur_delta:.4f}]'

            # If this is an outlier node, mark it as such.
            if jk_models[leaf.name]['is_outlier']:
                append_str += '  ***'

            txt_delta = ete3.TextFace(append_str, fsize=6)

            # influence = between 0 and thresh
            inf_heat_idx = int((min(abs(cur_inf), n_inf_max) / n_inf_max) * n_heatmap)
            txt_delta.background.color = heatmap[inf_heat_idx]

            # If this is an outlier node, mark it as such.
            if jk_models[leaf.name]['is_outlier']:
                txt_delta.border.width = 1

            leaf.add_face(txt_delta, column=at_col, position="aligned")

            # Include the contig lens
            contig_face = ete3.TextFace(f'{metadata.get_contig_len(leaf.name):,}', fsize=5)
            contig_face.margin_left = 2
            contig_face.background.color = heatmap[contig_len_to_idx(metadata.get_contig_len(leaf.name))]
            leaf.add_face(contig_face, column=at_col + 1, position='aligned')


def ete3_colour_clades(d_rank_lca, palette, tree):
    """Colours the clades fo the tree based on palette"""
    for i, (cur_rank, cur_lca_set) in enumerate(d_rank_lca.items()):
        if cur_rank not in palette:
            continue
        for cur_lca in cur_lca_set:
            ns = ete3.NodeStyle()
            ns['bgcolor'] = palette[cur_rank]
            if isinstance(cur_lca, tuple):
                n = tree.get_common_ancestor(cur_lca)
            else:
                n = tree.search_nodes(name=cur_lca)[0]
            n.set_style(ns)


def ete3_add_tax_to_leaf_nodes(tree, d_gid_to_tax, palette, highest_rank):
    """Add a new column for each leaf node in the tree that contains the
    taxonomy, it is coloured by the palette. Only the taxonomy that differs
    is shown.
    """
    for leaf in tree.traverse():
        if leaf.is_leaf():
            cur_tax_list = d_gid_to_tax[leaf.name]
            cols_used = 0
            for i in range(highest_rank, len(cur_tax_list)):
                cur_tax = cur_tax_list[i]
                new_face = ete3.TextFace(cur_tax, fsize=6)
                new_face.margin_left = 2
                if cur_tax in palette:
                    new_face.background.color = palette[cur_tax]
                leaf.add_face(new_face, column=cols_used, position='aligned')
                cols_used += 1


def draw_ete_reference_tree(unique_ranks, tree_str, d_gid_to_tax, palette,
                            metadata, out_dir, rep_label, q_label_outliers,
                            ete3_scale, inf_thresh, jk_models=None):
    """Draw the tree using ETE3.
    """
    # Load the tree.
    tree = ete3.Tree(tree_str, quoted_node_names=True, format=1)

    # Unique ranks at level
    d_ranks_at_level, highest_rank = count_taxonomy_ranks(unique_ranks)

    # Using the LCA, colour each monophyletic clade with the taxonomy.
    d_rank_lca = get_lca(tree_str, d_gid_to_tax)
    ete3_colour_clades(d_rank_lca, palette, tree)

    #  Add the taxonomy to each leaf node.
    ete3_add_tax_to_leaf_nodes(tree, d_gid_to_tax, palette, highest_rank)

    # Maybe add another column to each leaf node that contains the influence values.
    if jk_models:
        ete3_draw_influence(tree, jk_models, metadata, at_col=7 - highest_rank, n_inf_max=inf_thresh)

    ts = ete3.TreeStyle()
    ts.scale = ete3_scale

    # include the "common" ranks for each of the taxa as the title.
    common_ranks = ';'.join([list(d_ranks_at_level[x])[0] for x in range(highest_rank)])
    ts.title.add_face(ete3.TextFace(common_ranks, fsize=10), column=0)

    # Show the tree, or render to disk.
    if out_dir:
        for outlier_label in q_label_outliers:
            cur_dir = os.path.join(out_dir, outlier_label)
            if not os.path.isdir(cur_dir):
                os.makedirs(cur_dir)
            if jk_models:
                f_name = os.path.join(cur_dir, f'{rep_label}_qry.png')
            else:
                f_name = os.path.join(cur_dir, f'{rep_label}_rep.png')
            tree.render(f_name, w=2000, tree_style=ts)
    else:
        tree.show(tree_style=ts)


def retain_taxa_in_tree(path, taxa):
    """Keep only the specified taxa in the tree.

    Parameters
    ----------
    path : str
        The path to the tree in Newick format.
    taxa : List[str]
        A list of taxa to keep.

    Returns
    -------
    str
        The tree in Newick format.
    """
    t = dendropy.Tree.get_from_path(path, schema='newick', preserve_underscores=True)
    t.retain_taxa_with_labels(taxa)
    return t.as_string(schema='newick')


def write_model_to_disk(out_dir, unq_rep_idxs_2_filtered, r_labels, result_base, r_idx_to_k_q_idxs_mask, q_labels,
                        result_jk, inf_thresh, dif_thresh):
    """Write the jackknife influence value and model parameters to disk."""
    outliers = defaultdict(set)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(os.path.join(out_dir, 'model.tsv'), 'w') as fh:
        fh.write('reference_taxon\tquery_taxon\tr2\tmse\tgradient\tinfluence\tdelta\tis_outlier\n')
        for job_idx, rep_idx in enumerate(unq_rep_idxs_2_filtered):
            c = list()
            c.append(r_labels[rep_idx])
            c.append('NULL')
            c.append(result_base[job_idx, 0])
            c.append(result_base[job_idx, 1])
            c.append(result_base[job_idx, 2])
            c.append('NULL')
            c.append('NULL')
            c.append('NULL')
            fh.write('\t'.join(map(str, c)))
            fh.write('\n')

            num_jk = np.argwhere(r_idx_to_k_q_idxs_mask[rep_idx]).size
            for jk_idx in range(num_jk):
                is_cur_outlier = abs(result_jk[job_idx, 3, jk_idx]) >= inf_thresh and result_jk[job_idx, 4, jk_idx] > dif_thresh
                qry_idx = int(result_jk[job_idx, 5, jk_idx])

                d = list()
                d.append(r_labels[rep_idx])
                d.append(q_labels[qry_idx])
                d.append(result_jk[job_idx, 0, jk_idx])
                d.append(result_jk[job_idx, 1, jk_idx])
                d.append(result_jk[job_idx, 2, jk_idx])
                d.append(result_jk[job_idx, 3, jk_idx])
                d.append(result_jk[job_idx, 4, jk_idx])
                d.append(is_cur_outlier)
                fh.write('\t'.join(map(str, d)))
                fh.write('\n')

                if is_cur_outlier:
                    outliers[job_idx].add(jk_idx)
    return outliers


def plot_outliers(result_base, result_jk, d_outliers, unq_rep_idxs_2_filtered,
                  r_idx_to_k_q_idxs_mask, r_mat, q_mat, q_to_r_idx, qry_taxon_groups,
                  q_labels, r_labels, out_dir, path_r_tree, path_q_tree, metadata,
                  r_idx_to_k_r_idxs, ete3_scale, qry_sep, cpus, inf_thresh):
    """Create a scatterplot, and tree plot for all of the outliers."""

    for job_idx, jk_idxs in tqdm(d_outliers.items(), smoothing=0.1):

        # very hackish, but dont' need x and y for the tree so just force each fo the jk values into the dict
        all_jk_models = dict()
        for cur_jk_idx, jk_vals in enumerate(result_jk[job_idx].T):
            all_jk_models[q_labels[int(jk_vals[5])]] = {'gradient': jk_vals[2],
                                                        'r2': jk_vals[0],
                                                        'mse': jk_vals[1],
                                                        'influence': jk_vals[3],
                                                        'delta': jk_vals[4],
                                                        'is_outlier': cur_jk_idx in jk_idxs}

        # Generate the base model X and Y values (using all represented idxs)
        rep_idx = unq_rep_idxs_2_filtered[job_idx]
        cur_k_q_idxs = np.argwhere(r_idx_to_k_q_idxs_mask[rep_idx]).flatten()
        max_n_xy = cur_k_q_idxs.size ** 2

        base_x = np.empty(max_n_xy, dtype=np.float64)
        base_y = np.empty(max_n_xy, dtype=np.float64)
        base_qa = np.empty(max_n_xy, dtype=np.int64)
        base_qb = np.empty(max_n_xy, dtype=np.int64)
        base_ra = np.empty(max_n_xy, dtype=np.int64)
        base_rb = np.empty(max_n_xy, dtype=np.int64)

        n_points = set_xy_labels(base_x, base_y, base_qa, base_qb, base_ra, base_rb,
                                 cur_k_q_idxs, r_mat, q_mat, q_to_r_idx, qry_taxon_groups)

        base_x = base_x[0:n_points]
        base_y = base_y[0:n_points]
        base_qa = base_qa[0:n_points]
        base_qb = base_qb[0:n_points]
        base_ra = base_ra[0:n_points]
        base_rb = base_rb[0:n_points]

        rep_label = r_labels[rep_idx]

        base_model = {'gradient': result_base[job_idx, 2],
                      'r2': result_base[job_idx, 0],
                      'mse': result_base[job_idx, 1]}

        base_labels = np.empty((base_qa.size, 4), dtype=base_qa.dtype)
        base_labels[:, 0] = base_qa
        base_labels[:, 1] = base_qb
        base_labels[:, 2] = base_ra
        base_labels[:, 3] = base_rb

        base_q = q_labels[base_labels[:, [0, 1]]]
        base_r = r_labels[base_labels[:, [2, 3]]]
        base_xy = {'x': base_x, 'y': base_y, 'labels': np.concatenate((base_q, base_r), axis=1)}

        cur_k_r_idx = np.array(list(r_idx_to_k_r_idxs[rep_idx]))

        for jk_idx in jk_idxs:
            cur_jk_data = result_jk[job_idx, :, jk_idx]
            cur_idx = int(cur_jk_data[5])

            # Generate the comparison showing the base XY and the jk XY

            # Get all the idxs without this point
            cur_jk_idxs = cur_k_q_idxs[np.argwhere(cur_k_q_idxs != cur_idx).ravel()]
            assert (cur_k_q_idxs.size - 1 == cur_jk_idxs.size)

            jk_x = np.empty(max_n_xy, dtype=np.float64)
            jk_y = np.empty(max_n_xy, dtype=np.float64)
            jk_qa = np.empty(max_n_xy, dtype=np.int64)
            jk_qb = np.empty(max_n_xy, dtype=np.int64)
            jk_ra = np.empty(max_n_xy, dtype=np.int64)
            jk_rb = np.empty(max_n_xy, dtype=np.int64)

            n_jk_points = set_xy_labels(jk_x, jk_y, jk_qa, jk_qb, jk_ra, jk_rb,
                                        cur_jk_idxs, r_mat, q_mat, q_to_r_idx, qry_taxon_groups)

            jk_x = jk_x[0:n_jk_points]
            jk_y = jk_y[0:n_jk_points]
            jk_qa = jk_qa[0:n_jk_points]
            jk_qb = jk_qb[0:n_jk_points]
            jk_ra = jk_ra[0:n_jk_points]
            jk_rb = jk_rb[0:n_jk_points]

            # move to outer loop?
            q_label_outliers = set()
            q_label_outliers.add(q_labels[cur_idx])

            jk_sets = dict()
            jk_labels = np.empty((jk_qa.size, 4), dtype=jk_qa.dtype)
            jk_labels[:, 0] = jk_qa
            jk_labels[:, 1] = jk_qb
            jk_labels[:, 2] = jk_ra
            jk_labels[:, 3] = jk_rb
            jk_sets[q_labels[cur_idx]] = {'x': jk_x, 'y': jk_y, 'labels': jk_labels}

            jk_models = dict()
            jk_models[q_labels[cur_idx]] = {'gradient': result_jk[job_idx, 2, jk_idx],
                                            'r2': result_jk[job_idx, 0, jk_idx],
                                            'mse': result_jk[job_idx, 1, jk_idx],
                                            'influence': result_jk[job_idx, 3, jk_idx],
                                            'delta': result_jk[job_idx, 4, jk_idx]}

            generate_scatterplot(base_model, q_label_outliers, jk_models, out_dir, rep_label, jk_sets, base_xy)

            plt_tree(path_r_tree, path_q_tree, q_labels[cur_k_q_idxs], r_labels[cur_k_r_idx],
                     metadata, q_label_outliers,
                     out_dir, q_labels, r_labels, base_model, all_jk_models, rep_idx,
                     q_to_r_idx, rep_label, ete3_scale=ete3_scale, qry_sep=qry_sep, inf_thresh=inf_thresh)

            pass
