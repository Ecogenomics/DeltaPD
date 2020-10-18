#!python
#cython: language_level=3

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
from libc.stdlib cimport free
from deltapd.util cimport *
from deltapd.struct cimport *
from deltapd.model cimport set_relative_influence, set_xy


def test_fit_job():
    pass


def test_set_relative_influence(double[:, :] result_jk):
    set_relative_influence(result_jk)


def test_fit_model():
    pass


def test_set_xy(double[:] mv_x, double[:] mv_y, int[:] mv_qa,
                int[:] mv_qb, int[:] mv_ra, int[:] mv_rb, int[:] mv_q_idx,
                double[:, :] r_mat, double[:, :] q_mat, long[:] q_to_r_idx,
                long[:] q_to_label_start):

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
    mv_to_vec_idx(mv_q_idx, &cur_k_q_idx)

    # Run the method.
    set_xy(&vec_xy, &cur_k_q_idx, r_mat, q_mat, q_to_r_idx,
           q_to_label_start)

    # Copy the data back to the memory view.
    set_vec_dbl_to_mv(mv_x, &vec_x)
    set_vec_dbl_to_mv(mv_y, &vec_y)
    set_vec_idx_to_mv(mv_qa, &vec_qa)
    set_vec_idx_to_mv(mv_qb, &vec_qb)
    set_vec_idx_to_mv(mv_ra, &vec_ra)
    set_vec_idx_to_mv(mv_rb, &vec_rb)

    # Free the data.
    free(vec_xy.x.data)
    free(vec_xy.y.data)
    free(vec_xy.qa.data)
    free(vec_xy.qb.data)
    free(vec_xy.ra.data)
    free(vec_xy.rb.data)
    free(cur_k_q_idx.data)

