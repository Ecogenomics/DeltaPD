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

from deltapd.struct cimport VecDbl, VecXY, VecIdx
from libcpp cimport bool


cdef double get_mv_mean(double[:] mv) nogil
cdef void set_vec_without_index(VecIdx *cur_jk_set, VecIdx *cur_k_q_idx, Py_ssize_t jk_idx) nogil
cdef double get_stddev(VecDbl * x) nogil
cdef void set_y_hat(VecDbl *y_hat, VecDbl *y, double grad) nogil
cdef double calc_r2(VecDbl *y, VecDbl *y_hat) nogil
cdef double calc_mse_norm(VecDbl *y, VecDbl *y_hat) nogil
cdef double get_mean(VecDbl *vec) nogil
cdef double theil_sen_gradient(VecXY *vec_xy) nogil
cdef double get_median(VecDbl *vec) nogil
cdef double calc_mse(VecDbl *y, VecDbl *y_hat) nogil
cdef void set_mask_to_array(bool[:] mask, VecIdx *out) nogil

# Helper functions for tests.
cdef void mv_to_vec_idx(int[:] mv, VecIdx *vec)
cdef void mv_to_vec_idxx(long[:] mv, VecIdx *vec)
cdef void mv_to_vec_dbl(double[:] mv, VecDbl *vec)
cdef void set_vec_idx_to_mv(int[:] mv, VecIdx *vec)
cdef void set_vec_idxx_to_mv(long[:] mv, VecIdx *vec)
cdef void set_vec_dbl_to_mv(double[:] mv, VecDbl *vec)