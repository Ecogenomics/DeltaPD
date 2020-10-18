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

import logging
import os
import shutil
from collections import defaultdict
from datetime import datetime

import numpy as np
from jinja2 import FileSystemLoader, Environment

from deltapd import __version__
from deltapd.common import sha256_file
from deltapd.file import MetadataFile
from deltapd.method import get_pdms, get_mapping, get_knn, get_unq_rep_idx_mapping, \
    get_qry_taxon_groups, write_model_to_disk, plot_outliers, get_largest_jk_size
from deltapd.model import fit_c


def run(path_r_tree, path_q_tree, path_metadata, inf_thresh, dif_thresh, out_dir, k,
        plot, ete3_scale, max_taxa, qry_sep, path_msa_file, cpus):

    # Start logging.
    log = logging.getLogger('default')
    start_time = datetime.now()

    # Read the metadata file.
    metadata = MetadataFile(path_metadata, qry_sep)
    metadata.add_contig_lens(path_msa_file)

    # Dictionary containing keys required for the HTML report.
    template_data = {'outliers': defaultdict(list)}

    # Prepare the data.
    log.info('Generating phylogenetic distance matricies (PDMs).')
    r_pdm, r_labels, r_mat, q_pdm, q_labels, q_mat = get_pdms(path_r_tree, path_q_tree, cpus)

    log.info('Mapping query taxa to their representative taxa.')
    q_idx_to_r_idx, r_idx_to_q_idxs = get_mapping(r_labels, q_labels, metadata.get_rep_label_dict(), qry_sep)

    log.info(f'Generating the k={k} nearest neighbours for each reference taxon.')
    r_idx_to_k_r_idxs, r_idx_to_k_q_idxs_mask = get_knn(r_idx_to_q_idxs, r_mat, q_mat, k)

    log.info(f'Calculating model parameters using {cpus} CPU{"s" if cpus > 1 else ""}.')
    qry_taxon_groups = get_qry_taxon_groups(q_labels, qry_sep)
    unq_rep_idxs = get_unq_rep_idx_mapping(r_idx_to_k_r_idxs)
    largest_jk_size = get_largest_jk_size(r_idx_to_k_q_idxs_mask)

    log.info(f'Progressively removing outliers until no outliers are detected.')
    changes = list()
    i = 0
    base_dir = out_dir
    while True:
        # Set the current iteration parameters.
        out_dir = os.path.join(base_dir, 'deltapd_data', str(i))

        # Storage array for the base model parameters (r2, mse, grad)
        result_base = np.zeros((unq_rep_idxs.size, 3), dtype=np.float64)

        # Storage array for the jackknife result parameters (r2, mse, grad, inf, delta, qry_idx)
        result_jk = np.zeros((unq_rep_idxs.size, 6, largest_jk_size), dtype=np.float64)

        # Determine which representative taxa indexes represent fewer than MAX_Q_TAXA query taxa.
        # TODO: Find a way to include these taxa, for now this is done as there are a few cases
        # TODO: in the bacterial tree that cause it to blow out to an unreasonable time.
        rep_idx_lt_max_q = np.argwhere(np.sum(r_idx_to_k_q_idxs_mask[unq_rep_idxs, :], axis=1) <= max_taxa).ravel()
        unq_rep_idxs_filtered = unq_rep_idxs[rep_idx_lt_max_q]
        n_filtered = unq_rep_idxs.size - rep_idx_lt_max_q.size
        if n_filtered > 0:
            log.info(f'Filtered {n_filtered} representatives as they represented over {max_taxa} taxa.')

        # If outliers were detected in the previous run, omit them from future runs.
        if len(changes) > 0:
            for rep_idx, qry_idx in changes[-1]:
                r_idx_to_k_q_idxs_mask[rep_idx, qry_idx] = False

        log.info('Fitting the model.')
        fit_c(unq_rep_idxs_filtered, r_idx_to_k_q_idxs_mask, r_mat, q_mat, q_idx_to_r_idx, qry_taxon_groups,
              result_base, result_jk, cpus)

        log.info(f'Writing model to disk: {os.path.join(out_dir, "model.tsv")}')
        d_outliers = write_model_to_disk(out_dir, unq_rep_idxs_filtered, r_labels,
                                         result_base, r_idx_to_k_q_idxs_mask, q_labels,
                                         result_jk, inf_thresh, dif_thresh)

        # Iterate over all of the outliers detected and record that they were removed.
        cur_set = set()
        for job_idx, jk_idx_set in d_outliers.items():
            cur_rep_idx = unq_rep_idxs_filtered[job_idx]
            for jk_idx in jk_idx_set:
                qry_idx = int(result_jk[job_idx, 5, jk_idx])
                cur_set.add((cur_rep_idx, qry_idx))

                # Report these data.
                precision = 4
                cur_lst = template_data['outliers'][q_labels[qry_idx]]
                cur_lst.append({'i': i,
                                'ref_label': r_labels[cur_rep_idx],
                                'r2': round(result_jk[job_idx, 0, jk_idx], precision),
                                'mse': round(result_jk[job_idx, 1, jk_idx], precision),
                                'grad': round(result_jk[job_idx, 2, jk_idx], precision),
                                'inf': round(result_jk[job_idx, 3, jk_idx], precision),
                                'delta': round(result_jk[job_idx, 4, jk_idx], precision),
                                'base_r2': round(result_base[job_idx, 0], precision),
                                'base_mse': round(result_base[job_idx, 1], precision),
                                'base_grad': round(result_base[job_idx, 2], precision)})

        changes.append(frozenset(cur_set))

        if plot:
            log.info('Creating plots.')
            plot_outliers(result_base, result_jk, d_outliers, unq_rep_idxs_filtered, r_idx_to_k_q_idxs_mask,
                          r_mat, q_mat, q_idx_to_r_idx, qry_taxon_groups,
                          q_labels, r_labels, out_dir, path_r_tree, path_q_tree, metadata, r_idx_to_k_r_idxs,
                          ete3_scale=ete3_scale, qry_sep=qry_sep, cpus=cpus, inf_thresh=inf_thresh)

        # If no changes were detected, stop processing.
        if len(changes[-1]) == 0:
            break

        # Increment iteration variables.
        i += 1

    # Done - record endtime for report.
    end_time = datetime.now()

    # Report input arguments.
    template_data['path_r_tree'] = path_r_tree
    template_data['path_r_tree_sha256'] = sha256_file(path_r_tree)
    template_data['path_q_tree'] = path_q_tree
    template_data['path_q_tree_sha256'] = sha256_file(path_q_tree)
    template_data['path_metadata'] = path_metadata
    template_data['path_metadata_sha256'] = sha256_file(path_metadata)
    template_data['inf_thresh'] = inf_thresh
    template_data['dif_thresh'] = dif_thresh
    template_data['out_dir'] = base_dir
    template_data['k'] = k
    template_data['plot'] = plot
    template_data['ete3_scale'] = ete3_scale
    template_data['max_taxa'] = max_taxa
    template_data['qry_sep'] = qry_sep
    template_data['path_msa_file'] = path_msa_file
    template_data['path_msa_file_sha256'] = sha256_file(path_msa_file)
    template_data['cpus'] = cpus

    # Report metdata.
    template_data['exec_time_start'] = start_time
    template_data['exec_time_end'] = end_time
    template_data['exec_time_duration'] = end_time - start_time
    template_data['deltapd_version'] = __version__

    # Copy static files for the HTML report
    raw_dir = os.path.join(base_dir, 'deltapd_data')
    report_dir = os.path.join(base_dir, 'report_data')
    template_dir = os.path.join('deltapd', 'templates')

    if not os.path.isdir(raw_dir):
        os.makedirs(raw_dir)
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)
    if os.path.isdir(os.path.join(report_dir, 'static')):
        shutil.rmtree(os.path.join(report_dir, 'static'))
    shutil.copytree(os.path.join(template_dir, 'static'), os.path.join(report_dir, 'static'))

    # Generate the HTML report.
    with open(os.path.join(base_dir, 'report.html'), 'w') as fh:
        template_env = Environment(loader=FileSystemLoader(template_dir))
        template = template_env.get_template('index.html')
        fh.write(template.render(**template_data))

    return
