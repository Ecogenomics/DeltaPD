#cython: language_level=3, boundscheck=False, wraparound=False

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

import cython
from cython.parallel import parallel, prange
from libc.stdlib cimport malloc, free, realloc
from libc.math cimport sqrt
from deltapd.util cimport *
from tqdm import tqdm
from deltapd.struct cimport *
from libcpp cimport bool


def fit_c(long[:] unq_rep_idxs,
            bool[:, :] r_idx_to_k_q_idxs_mask,
            double[:, :] r_mat,
            double[:, :] q_mat,
            long[:] q_idx_to_r_idx,
            long[:] qry_taxon_groups,
            double[:, :] result_base,
            double[:, :, :] result_jk,
            Py_ssize_t cpus):
    """The main method, schedules jobs with OpenMP."""
    cdef Py_ssize_t job_idx, rep_idx
    cdef Py_ssize_t n_reps = unq_rep_idxs.size
    with tqdm(total=n_reps, smoothing=0.1) as p_bar:
        for job_idx in prange(n_reps, nogil=True, schedule='dynamic',
                              num_threads=cpus, chunksize=1):
            rep_idx = unq_rep_idxs[job_idx]
            fit_job(r_idx_to_k_q_idxs_mask[rep_idx],
                          r_mat, q_mat, q_idx_to_r_idx,
                          qry_taxon_groups,
                          result_base[job_idx],
                          result_jk[job_idx])
            with gil:
                p_bar.update()


cdef void fit_job(bool[:] r_to_k_q_idxs_mask,
                  double[:, :] r_mat,
                  double[:, :] q_mat,
                  long[:] q_idx_to_r_idx,
                  long[:] qry_taxon_groups,
                  double[:] result_base,
                  double[:, :] result_jk) nogil:
    """For a given representative, and the query genomes it represents. Do work."""

    # From the given slice, determine what query indexes are present.
    cdef VecIdx cur_k_q_idx
    set_mask_to_array(r_to_k_q_idxs_mask, &cur_k_q_idx)

    # Base model parameters (constant).
    cdef ModelParam base_param      # Parameters
    cdef VecXY base_xy              # X/Y coordinates.
    cdef VecDbl base_x, base_y
    cdef VecIdx base_qa, base_qb, base_ra, base_rb
    base_xy.x = &base_x
    base_xy.y = &base_y
    base_xy.qa = &base_qa
    base_xy.qb = &base_qb
    base_xy.ra = &base_ra
    base_xy.rb = &base_rb

    # Jackknife model parameters (overwritten each iteration).
    cdef ModelParam jk_param        # Parameters
    cdef VecXY cur_jk_xy            # X/Y coordinates.
    cdef VecDbl cur_jk_x, cur_jk_y
    cdef VecIdx cur_jk_qa, cur_jk_qb, cur_jk_ra, cur_jk_rb
    cur_jk_xy.x = &cur_jk_x
    cur_jk_xy.y = &cur_jk_y
    cur_jk_xy.qa = &cur_jk_qa
    cur_jk_xy.qb = &cur_jk_qb
    cur_jk_xy.ra = &cur_jk_ra
    cur_jk_xy.rb = &cur_jk_rb

    # Will be used when determining each jackknife set in the loop.
    cdef VecIdx cur_jk_set
    cur_jk_set.size = cur_k_q_idx.size - 1
    cur_jk_set.data = <Py_ssize_t*> malloc(cur_jk_set.size * sizeof(Py_ssize_t))
    if not cur_jk_set.data:
        with gil:
            raise MemoryError()

    # Set the base model parameters.
    set_xy(&base_xy, &cur_k_q_idx, r_mat, q_mat, q_idx_to_r_idx, qry_taxon_groups)
    fit_model(&base_param, &base_xy)
    result_base[0] = base_param.r2    # r2
    result_base[1] = base_param.mse   # mse
    result_base[2] = base_param.grad  # gradient

    # Jackknife each taxon.
    cdef Py_ssize_t jk_idx
    for jk_idx in range(cur_k_q_idx.size):

        # Generate the jackknife set (to an already allocated array)
        set_vec_without_index(&cur_jk_set, &cur_k_q_idx, jk_idx)

        # Set the jackknife model parameters.
        set_xy(&cur_jk_xy, &cur_jk_set, r_mat, q_mat, q_idx_to_r_idx, qry_taxon_groups)
        fit_model(&jk_param, &cur_jk_xy)
        result_jk[0, jk_idx] = jk_param.r2                    # r2
        result_jk[1, jk_idx] = jk_param.mse                   # mse
        result_jk[2, jk_idx] = jk_param.grad                  # gradient
        result_jk[4, jk_idx] = base_param.mse - jk_param.mse  # delta
        result_jk[5, jk_idx] = cur_k_q_idx.data[jk_idx]       # q_idx omitted.

        # Free memory.
        free(cur_jk_xy.x.data)
        free(cur_jk_xy.y.data)
        free(cur_jk_xy.qa.data)
        free(cur_jk_xy.qb.data)
        free(cur_jk_xy.ra.data)
        free(cur_jk_xy.rb.data)

    # Set the relative influence information, using the existing data.
    set_relative_influence(result_jk[:, 0:cur_k_q_idx.size])  # result_jk 3

    # Free memory.
    free(cur_k_q_idx.data)
    free(cur_jk_set.data)
    free(base_xy.x.data)
    free(base_xy.y.data)
    free(base_xy.qa.data)
    free(base_xy.qb.data)
    free(base_xy.ra.data)
    free(base_xy.rb.data)


cdef void set_relative_influence(double[:, :] result_jk) nogil:
    """Determine the relative influence of each data point in the jackknife
    set. Efron (1992)."""
    cdef double[:] metric_v = result_jk[1, :]  # mse
    cdef double x_bar = get_mv_mean(metric_v)
    cdef Py_ssize_t n = result_jk.shape[1]
    cdef Py_ssize_t i
    cdef double numerator, denominator

    denominator = 0.0
    for i in range(n):
        denominator += ((n - 1) * (x_bar - metric_v[i])) ** 2
    denominator = sqrt(denominator / (n - 1))

    for i in range(n):
        numerator = (n - 1) * (x_bar - metric_v[i])
        result_jk[3, i] =  numerator / denominator


cdef void fit_model(ModelParam *model_param, VecXY *vec_xy) nogil:
    """Given X/Y coordinates, calculate the model parameters."""

    # Calculate model parameters.
    cdef double gradient = theil_sen_gradient(vec_xy)

    # Calculate the estimate using the model gradient.
    cdef VecDbl y_hat
    y_hat.size = vec_xy.x.size
    y_hat.data = <double*> malloc(y_hat.size * sizeof(double))
    if not y_hat.data:
        with gil:
            raise MemoryError()
    try:
        set_y_hat(&y_hat, vec_xy.x, gradient)

        # Set the output variables.
        model_param.r2 = calc_r2(vec_xy.y, &y_hat)
        model_param.mse = calc_mse_norm(vec_xy.y, &y_hat)
        model_param.grad = gradient
    finally:
        free(y_hat.data)


def set_xy_labels(double[:] mv_x, double[:] mv_y, long[:] mv_qa, long[:] mv_qb,
                long[:] mv_ra, long[:] mv_rb, long[:] mv_q_idx,
                double[:, :] r_mat, double[:, :] q_mat, long[:] q_to_r_idx,
                long[:] q_to_label_start):
    """Determine the X and Y coordinates, plus label for each query to ref."""

    # Set the storage vectors.
    cdef VecXY vec_xy
    cdef VecDbl vec_x, vec_y
    cdef VecIdx vec_qa, vec_qb, vec_ra, vec_rb
    vec_xy.x = &vec_x
    vec_xy.y = &vec_y
    vec_xy.qa = &vec_qa
    vec_xy.qb = &vec_qb
    vec_xy.ra = &vec_ra
    vec_xy.rb = &vec_rb

    # Translate the index.
    cdef VecIdx cur_k_q_idx
    mv_to_vec_idxx(mv_q_idx, &cur_k_q_idx)

    # Run the method.
    set_xy(&vec_xy, &cur_k_q_idx, r_mat, q_mat, q_to_r_idx,
           q_to_label_start)

    # Copy the data back to the memory view.
    set_vec_dbl_to_mv(mv_x[:], &vec_x)
    set_vec_dbl_to_mv(mv_y[:], &vec_y)
    set_vec_idxx_to_mv(mv_qa[:], &vec_qa)
    set_vec_idxx_to_mv(mv_qb[:], &vec_qb)
    set_vec_idxx_to_mv(mv_ra[:], &vec_ra)
    set_vec_idxx_to_mv(mv_rb[:], &vec_rb)

    cdef Py_ssize_t n_points = vec_xy.x.size

    # Free the data.
    free(vec_xy.x.data)
    free(vec_xy.y.data)
    free(vec_xy.qa.data)
    free(vec_xy.qb.data)
    free(vec_xy.ra.data)
    free(vec_xy.rb.data)
    free(cur_k_q_idx.data)

    return n_points


cdef void set_xy(VecXY *vec_xy, VecIdx *cur_k_q_idx, double[:, :] r_mat, double[:, :] q_mat,
             long[:] q_to_r_idx, long[:]  q_to_label_start) nogil:
    """Calculate the X/Y coordinates given the set of query indexes."""

    # Allocate memory for the output arrays.
    cdef Py_ssize_t arr_size = cur_k_q_idx.size
    cdef double* out_x = <double*> malloc(arr_size * sizeof(double))
    cdef double* out_y = <double*> malloc(arr_size * sizeof(double))
    cdef Py_ssize_t* out_qa = <Py_ssize_t*> malloc(arr_size * sizeof(Py_ssize_t))
    cdef Py_ssize_t* out_qb = <Py_ssize_t*> malloc(arr_size * sizeof(Py_ssize_t))
    cdef Py_ssize_t* out_ra = <Py_ssize_t*> malloc(arr_size * sizeof(Py_ssize_t))
    cdef Py_ssize_t* out_rb = <Py_ssize_t*> malloc(arr_size * sizeof(Py_ssize_t))

    if not (out_x and out_y and out_qa and out_qb and out_ra and out_rb):
        with gil:
            raise MemoryError()

    # Setup variables for looping.
    cdef Py_ssize_t n_comb, a_idx, b_idx, cur_qa, cur_ra, cur_qb, cur_rb, cur_idx
    cdef Py_ssize_t n_indices = cur_k_q_idx.size

    cdef Py_ssize_t cnt = 0
    for a_idx in range(n_indices):
        cur_qa = cur_k_q_idx.data[a_idx]  # query a index
        cur_ra = q_to_r_idx[cur_qa]       # query a rep index

        for b_idx in range(a_idx):

            cur_qb = cur_k_q_idx.data[b_idx]  # query b index
            cur_rb = q_to_r_idx[cur_qb]       # query b rep index

            # Keep so long as it's not comparing the same genome.
            if cur_ra != cur_rb and q_to_label_start[cur_qa] != q_to_label_start[cur_qb]:

                # Check if more memory is required.
                if arr_size - 5 <= cnt:
                    arr_size *= 2
                    out_x = <double*> realloc(out_x, arr_size * sizeof(double))
                    out_y = <double*> realloc(out_y, arr_size * sizeof(double))
                    out_qa = <Py_ssize_t*> realloc(out_qa, arr_size * sizeof(Py_ssize_t))
                    out_qb = <Py_ssize_t*> realloc(out_qb, arr_size * sizeof(Py_ssize_t))
                    out_ra = <Py_ssize_t*> realloc(out_ra, arr_size * sizeof(Py_ssize_t))
                    out_rb = <Py_ssize_t*> realloc(out_rb, arr_size * sizeof(Py_ssize_t))
                    if not (out_x and out_y and out_qa and out_qb and out_ra and out_rb):
                        with gil:
                            raise MemoryError()

                # Store the values.
                out_y[cnt] = q_mat[cur_qa, cur_qb]
                out_x[cnt] = r_mat[cur_ra, cur_rb]
                out_qa[cnt] = cur_qa
                out_qb[cnt] = cur_qb
                out_ra[cnt] = cur_ra
                out_rb[cnt] = cur_rb
                cnt += 1

    # Set the output data.
    vec_xy.x.size = cnt
    vec_xy.y.size = cnt
    vec_xy.qa.size = cnt
    vec_xy.qb.size = cnt
    vec_xy.ra.size = cnt
    vec_xy.rb.size = cnt

    vec_xy.x.data = out_x
    vec_xy.y.data = out_y
    vec_xy.qa.data = out_qa
    vec_xy.qb.data = out_qb
    vec_xy.ra.data = out_ra
    vec_xy.rb.data = out_rb
