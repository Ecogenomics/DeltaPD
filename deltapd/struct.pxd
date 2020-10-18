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

cdef struct VecDbl:
    Py_ssize_t size
    double* data

cdef struct VecIdx:
    Py_ssize_t size
    Py_ssize_t* data

cdef struct VecXY:
    VecDbl* x
    VecDbl* y
    VecIdx* qa
    VecIdx* qb
    VecIdx* ra
    VecIdx* rb

cdef struct ModelParam:
    double grad
    double r2
    double mse