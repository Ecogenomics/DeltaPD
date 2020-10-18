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

import shutil
import tempfile
import unittest

from sklearn import metrics
from sklearn.linear_model import TheilSenRegressor
import numpy as np

from deltapd.util import sort_mat_rows
from deltapd.util_test import *

PRECISION = 8
N_TESTS = 500


class TestUtil(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix='deltapd_test_')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_get_mv_mean(self):
        for n in range(1, N_TESTS):
            x = np.random.rand(n)
            true = np.mean(x)
            test = test_get_mv_mean(x)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_set_vec_without_index(self):
        vals = np.random.permutation(np.random.randint(0, N_TESTS * 10, N_TESTS))
        vals = vals.astype(np.int32)
        for i in range(len(vals)):
            true = np.array([x for ii, x in enumerate(vals) if ii != i])
            test = np.empty(true.size, dtype=np.int32)
            test_set_vec_without_index(test, vals, i)
            self.assertTrue(np.all(true == test))

    def test_get_stddev(self):
        for n in range(1, N_TESTS):
            x = np.random.rand(n)
            true = np.std(x)
            test = test_get_stddev(x)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_set_y_hat(self):
        for n in range(1, N_TESTS):
            grad = np.random.randint(0, 100) / 50
            x = np.random.rand(n)
            true = x * grad

            test = np.empty(x.size)
            test_set_y_hat(test, x, grad)

            self.assertTrue(np.allclose(true, test))

    def test_calc_r2(self):
        for n in range(2, N_TESTS):
            x = np.random.rand(n)
            y = np.random.rand(n)

            grad = np.random.randint(0, 100) / 50
            y_hat = x.flatten() * grad
            true = metrics.r2_score(y, y_hat)
            test = test_calc_r2(y, y_hat)

            self.assertAlmostEqual(true, test, PRECISION)

    def test_get_mse(self):
        for n in range(2, N_TESTS):
            x = np.random.rand(n)
            y = np.random.rand(n)

            grad = np.random.randint(0, 100) / 50
            y_hat = x.flatten() * grad
            # true = np.sqrt(metrics.mean_squared_error(y, y_hat)) / np.std(y)
            true = metrics.mean_squared_error(y, y_hat)
            test = test_calc_mse(y, y_hat)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_get_mse_norm(self):
        for n in range(2, N_TESTS):
            x = np.random.rand(n)
            y = np.random.rand(n)

            grad = np.random.randint(0, 100) / 50
            y_hat = x.flatten() * grad
            true = np.sqrt(metrics.mean_squared_error(y, y_hat)) / np.std(y)
            test = test_calc_mse_norm(y, y_hat)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_get_mean(self):
        for n in range(1, N_TESTS):
            x = np.random.rand(n)
            true = np.mean(x)
            test = test_get_mean(x)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_theil_sen_gradient(self):
        for n in range(1, N_TESTS):
            x = np.random.rand(n)
            y = np.random.rand(n)

            true = TheilSenRegressor(fit_intercept=False, max_subpopulation=n,
                                     max_iter=1e9, tol=1e-9).fit(x.reshape(-1, 1), y)
            true_grad = true.coef_
            # test_alt = np.get_median(np.divide(y, x, out=np.zeros_like(y), where=x != 0))
            test_grad = test_theil_sen_gradient(x, y)
            self.assertAlmostEqual(true_grad, test_grad, PRECISION)

    def test_median(self):
        for n in range(1, N_TESTS):
            arr = np.random.rand(n)
            true = np.median(arr)
            test = test_get_median(arr)
            self.assertAlmostEqual(true, test, PRECISION)

    def test_set_mask_to_array(self):

        # Check the base case of all False.
        x = np.zeros(3, dtype=np.bool)
        out = np.full(3, -1, dtype=np.int32)
        test_set_mask_to_array(x, out)
        self.assertTrue(np.all(out == -1))

        # Check random amounts (40% of array is True)
        for n in range(1, N_TESTS):
            x = np.zeros(n, dtype=np.bool)
            idx = np.random.permutation(range(n))[0:int(np.ceil(n * 0.4))]
            x[idx] = True

            out = np.full(x.size, -1, dtype=np.int32)
            test_set_mask_to_array(x, out)

            test = out[0:idx.size]
            true = np.argwhere(x).ravel()

            self.assertTrue(np.all(true == test))

    def test_sort_mat_rows(self):

        n = 5
        mat = np.empty((n, n), dtype=np.float64)
        indices = np.triu_indices_from(mat, 1)
        mat[indices] = np.random.random(indices[0].size)
        mat[np.diag_indices_from(mat)] = 0
        mat = mat + mat.T

        true = np.argsort(mat, axis=1)

        test = np.empty((n, n), dtype=np.float64)
        sort_mat_rows(mat, test)


        return



if __name__ == '__main__':
    unittest.main()
