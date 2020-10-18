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
from libc.stdlib cimport free
from deltapd.util cimport *


def test_get_mv_mean(double[:] mv):
    return get_mv_mean(mv)

def test_set_vec_without_index(int[:] mv_out, int[:] mv_all, Py_ssize_t omit):
    cdef VecIdx vec_out
    mv_to_vec_idx(mv_out, &vec_out)

    cdef VecIdx vec_all
    mv_to_vec_idx(mv_all, &vec_all)

    set_vec_without_index(&vec_out, &vec_all, omit)
    set_vec_idx_to_mv(mv_out, &vec_out)

    free(vec_all.data)
    free(vec_out.data)
    return

def test_get_stddev(double[:] mv):
    cdef VecDbl vec
    mv_to_vec_dbl(mv, &vec)

    cdef double std = get_stddev(&vec)
    free(vec.data)
    return std


def test_set_y_hat(double[:] mv_y_hat, double[:] mv_x, double grad):
    cdef VecDbl vec_y_hat
    mv_to_vec_dbl(mv_y_hat, &vec_y_hat)

    cdef VecDbl vec_x
    mv_to_vec_dbl(mv_x, &vec_x)

    set_y_hat(&vec_y_hat, &vec_x, grad)
    set_vec_dbl_to_mv(mv_y_hat, &vec_y_hat)

    free(vec_y_hat.data)
    free(vec_x.data)

def test_calc_r2(double[:] y, double[:] y_hat):
    cdef VecDbl y_vec
    mv_to_vec_dbl(y, &y_vec)

    cdef VecDbl y_hat_vec
    mv_to_vec_dbl(y_hat, &y_hat_vec)

    cdef double r2 = calc_r2(&y_vec, &y_hat_vec)

    free(y_hat_vec.data)
    free(y_vec.data)
    return r2

def test_calc_mse(double[:] y, double[:] y_hat):
    cdef VecDbl y_vec
    mv_to_vec_dbl(y, &y_vec)

    cdef VecDbl y_hat_vec
    mv_to_vec_dbl(y_hat, &y_hat_vec)

    cdef double mse = calc_mse(&y_vec, &y_hat_vec)

    free(y_hat_vec.data)
    free(y_vec.data)
    return mse

def test_calc_mse_norm(double[:] y, double[:] y_hat):
    cdef VecDbl y_vec
    mv_to_vec_dbl(y, &y_vec)

    cdef VecDbl y_hat_vec
    mv_to_vec_dbl(y_hat, &y_hat_vec)

    cdef double mse = calc_mse_norm(&y_vec, &y_hat_vec)

    free(y_hat_vec.data)
    free(y_vec.data)
    return mse

def test_get_mean(double[:] mv):
    cdef VecDbl vec
    mv_to_vec_dbl(mv, &vec)

    cdef mean = get_mean(&vec)
    free(vec.data)
    return mean

def test_theil_sen_gradient(double[:] x, double[:] y):
    cdef VecXY vec_xy
    cdef VecDbl vec_x
    cdef VecDbl vec_y

    vec_xy.x = &vec_x
    vec_xy.y = &vec_y

    mv_to_vec_dbl(x, &vec_x)
    mv_to_vec_dbl(y, &vec_y)

    cdef double grad = theil_sen_gradient(&vec_xy)

    free(vec_x.data)
    free(vec_y.data)
    return grad


def test_get_median(double[:] mv):
    cdef VecDbl vec
    mv_to_vec_dbl(mv, &vec)

    cdef median = get_median(&vec)
    free(vec.data)
    return median

def test_set_mask_to_array(bool[:] mv_mask, int[:] mv_out):
    cdef VecIdx vec_out
    set_mask_to_array(mv_mask, &vec_out)

    set_vec_idx_to_mv(mv_out, &vec_out)
    free(vec_out.data)
