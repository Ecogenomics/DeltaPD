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


from deltapd.util cimport *
from deltapd.struct cimport *
from libcpp cimport bool


cdef void fit_job(bool[:] r_to_k_q_idxs_mask,
                  double[:, :] r_mat,
                  double[:, :] q_mat,
                  long[:] q_idx_to_r_idx,
                  long[:] qry_taxon_groups,
                  double[:] result_base,
                  double[:, :] result_jk) nogil

cdef void set_relative_influence(double[:, :] result_jk) nogil

cdef void fit_model(ModelParam *model_param, VecXY *vec_xy) nogil

cdef void set_xy(VecXY *vec_xy, VecIdx *cur_k_q_idx, double[:, :] r_mat,
                 double[:, :] q_mat, long[:] q_to_r_idx,
                 long[:] q_to_label_start) nogil