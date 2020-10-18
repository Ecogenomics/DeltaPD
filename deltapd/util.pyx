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

from libc.math cimport floor, sqrt
from libcpp.algorithm cimport sort as stdsort
from deltapd.struct cimport VecDbl, VecXY, VecIdx
from libc.stdlib cimport malloc, free, realloc


cdef double get_mv_mean(double[:] mv) nogil:
    """Calculate the mean of a memoryview."""
    cdef Py_ssize_t i, n
    cdef double cur = 0
    n = mv.shape[0]
    for i in range(n):
        cur += mv[i]
    return cur / n


cdef void set_vec_without_index(VecIdx *cur_jk_set, VecIdx *cur_k_q_idx, Py_ssize_t jk_idx) nogil:
    """Overrides all values in a vector contain all values that are not jk_idx, from cur_k_qidx."""
    cdef Py_ssize_t i, cur = 0
    for i in range(cur_k_q_idx.size):
        if i != jk_idx:
            cur_jk_set.data[cur] = cur_k_q_idx.data[i]
            cur += 1


cdef double get_stddev(VecDbl * x) nogil:
    """Calculate the standard deviation of a vector."""
    cdef double tot = 0.0
    cdef double mean = get_mean(x)
    cdef Py_ssize_t i
    for i in range(x.size):
        tot += (x.data[i] - mean) ** 2
    return sqrt((1/x.size) * tot)


cdef void set_y_hat(VecDbl *y_hat, VecDbl *x, double grad) nogil:
    """Override the values in y_hat to be the estimate given params"""
    cdef Py_ssize_t i
    for i in range(x.size):
        y_hat.data[i] = x.data[i] * grad


cdef double calc_r2(VecDbl *y, VecDbl *y_hat) nogil:
    """Calculate the correlation coefficient."""
    cdef double numerator = 0.0
    cdef double denominator = 0.0
    cdef double y_avg = get_mean(y)
    cdef Py_ssize_t i
    for i in range(y.size):
        numerator += (y.data[i] - y_hat.data[i]) ** 2
        denominator += (y.data[i] - y_avg) ** 2
    return 1 - (numerator / denominator)


cdef double calc_mse(VecDbl *y, VecDbl *y_hat) nogil:
    """Calculate the get_mean squared error."""
    cdef double total = 0.0
    cdef Py_ssize_t i
    for i in range(y.size):
        total += (y.data[i] - y_hat.data[i]) ** 2
    return total / y.size


cdef double calc_mse_norm(VecDbl *y, VecDbl *y_hat) nogil:
    """Calculate the normalised MSE"""
    return sqrt(calc_mse(y, y_hat)) / get_stddev(y)


cdef double get_mean(VecDbl *vec) nogil:
    """Calculate the get_mean of a vector."""
    cdef double tot = 0.0
    cdef Py_ssize_t i
    for i in range(vec.size):
        tot += vec.data[i]
    return tot / vec.size


cdef double theil_sen_gradient(VecXY *vec_xy) nogil:
    """Calculate the gradient of X/Y coordinates using the Theil-Sen method."""
    cdef Py_ssize_t i
    cdef VecDbl frac
    frac.size = vec_xy.x.size
    frac.data = <double*> malloc(frac.size * sizeof(double))
    if not frac.data:
        with gil:
            raise MemoryError()
    try:
        for i in range(frac.size):
            if vec_xy.x.data[i] == 0:
                frac.data[i] = 0.0
            else:
                frac.data[i] = vec_xy.y.data[i] / vec_xy.x.data[i]
        return get_median(&frac)
    finally:
        free(frac.data)


cdef double get_median(VecDbl *vec) nogil:
    """Calculate the get_median of a vector, note: it will be sorted in-place."""
    if vec.size == 1:
        return vec.data[0]
    stdsort(&vec.data[0], &vec.data[vec.size])
    cdef Py_ssize_t idx_half = <Py_ssize_t>floor(vec.size / 2) - 1
    if vec.size % 2 == 0:
        return (vec.data[idx_half] + vec.data[idx_half+1]) / 2
    else:
        return vec.data[idx_half + 1]


cdef void set_mask_to_array(bool[:] mask, VecIdx *out) nogil:
    """Convert a vector of boolean values to a Vec struct containing
    the indices that they were True."""
    cdef Py_ssize_t size = 50
    cdef Py_ssize_t* data = <Py_ssize_t*> malloc(size * sizeof(Py_ssize_t))
    if not data:
        with gil:
            raise MemoryError()

    # Iterate over the mask.
    cdef Py_ssize_t i, cnt = 0
    for i in range(mask.shape[0]):
        if mask[i] is True:
            # Allocate more memory if required.
            if size - 5 <= cnt:
                size *= 2
                data = <Py_ssize_t*> realloc(data, size * sizeof(Py_ssize_t))
                if not data:
                    with gil:
                        raise MemoryError()
            # Store the value.
            data[cnt] = i
            cnt += 1

    # Set the output variables.
    out.size = cnt
    out.data = data


cdef void mv_to_vec_idx(int[:] mv, VecIdx *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    vec.size = mv.size
    vec.data = <Py_ssize_t*> malloc(vec.size * sizeof(Py_ssize_t))
    if not vec.data:
        raise MemoryError()
    for i in range(vec.size):
        vec.data[i] = mv[i]


cdef void mv_to_vec_idxx(long[:] mv, VecIdx *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    vec.size = mv.size
    vec.data = <Py_ssize_t*> malloc(vec.size * sizeof(Py_ssize_t))
    if not vec.data:
        raise MemoryError()
    for i in range(vec.size):
        vec.data[i] = mv[i]


cdef void mv_to_vec_dbl(double[:] mv, VecDbl *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    vec.size = mv.size
    vec.data = <double *> malloc(vec.size * sizeof(double))
    if not vec.data:
            raise MemoryError()
    for i in range(vec.size):
        vec.data[i] = mv[i]


cdef void set_vec_idx_to_mv(int[:] mv, VecIdx *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    for i in range(vec.size):
        mv[i] = vec.data[i]


cdef void set_vec_idxx_to_mv(long[:] mv, VecIdx *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    for i in range(vec.size):
        mv[i] = vec.data[i]


cdef void set_vec_dbl_to_mv(double[:] mv, VecDbl *vec):
    """A test helper function"""
    cdef Py_ssize_t i
    for i in range(vec.size):
        mv[i] = vec.data[i]


cpdef void sort_mat_rows(double[:, ::1] mat_in, double[:, :] mat_out) nogil:
    cdef Py_ssize_t row_idx
    for row_idx in range(mat_in.shape[0]):
        stdsort(&mat_in[row_idx, 0], (&mat_in[row_idx, 0]) + mat_in.shape[0])
        pass
    return
