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

import numpy as np

from deltapd.model_test import *

PRECISION = 8
N_TESTS = 500


class TestUtil(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix='deltapd_test_')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_fit_job(self):
        pass

    def test_test_set_relative_influence(self):
        for n in range(2, N_TESTS):
            # Generate random doubles for the calculation.
            v = np.random.random(n)

            # Calculate the influence.
            x_bar = np.mean(v)
            true = ((n - 1) * (x_bar - v)) / np.sqrt(np.sum(np.power((n - 1) * (x_bar - v), 2)) / (n - 1))

            # r2, mse, grad, delta, q_idx
            test_arr = np.zeros((6, n), dtype=np.float64)
            test_arr[1, :] = v
            test_set_relative_influence(test_arr)

            # Check
            self.assertTrue(np.allclose(test_arr[3, :], true))

    def test_test_fit_model(self):
        pass

    def test_test_set_xy(self):

        """
        Rules: Will always get an X and Y value back, unless:
        1. The genome name is the same, e.g. q_to_label_start[x] == "..."[y]
        2. The representative genome indexes are different.
        """

        for _ in range(N_TESTS):

            sizes = np.random.randint(10, 500, 2)

            n_rep = sizes.min()
            n_qry = sizes.max()

            r_mat = get_test_matrix(n_rep)
            q_mat = get_test_matrix(n_qry)

            # Get the indexes of the query taxa that are being used.
            n_qry_to_use = np.random.randint(1, n_qry)
            cur_k_q_idx = np.sort(np.random.permutation(n_qry)[0:n_qry_to_use]).astype(np.int32)

            # Randomly select representative genome to represent these query taxa.
            q_to_r_idx = np.random.randint(0, n_rep, n_qry)

            # Randomly select a grouping where the same number equals the same genome source.
            q_to_label_start = np.random.randint(0, n_qry, n_qry)

            # Output arrays at the maximum possible size.
            mv_x = np.zeros(n_qry ** 2, dtype=np.float64)
            mv_y = np.zeros(n_qry ** 2, dtype=np.float64)
            mv_qa = np.zeros(n_qry ** 2, dtype=np.int32)
            mv_qb = np.zeros(n_qry ** 2, dtype=np.int32)
            mv_ra = np.zeros(n_qry ** 2, dtype=np.int32)
            mv_rb = np.zeros(n_qry ** 2, dtype=np.int32)

            true_x = np.zeros(n_qry ** 2, dtype=np.float64)
            true_y = np.zeros(n_qry ** 2, dtype=np.float64)
            true_qa = np.zeros(n_qry ** 2, dtype=np.int32)
            true_qb = np.zeros(n_qry ** 2, dtype=np.int32)
            true_ra = np.zeros(n_qry ** 2, dtype=np.int32)
            true_rb = np.zeros(n_qry ** 2, dtype=np.int32)
            cnt = 0
            for a_idx in range(cur_k_q_idx.size):
                cur_qa = cur_k_q_idx[a_idx]
                cur_ra = q_to_r_idx[cur_qa]
                for b_idx in range(a_idx):
                    cur_qb = cur_k_q_idx[b_idx]
                    cur_rb = q_to_r_idx[cur_qb]

                    if cur_ra != cur_rb and q_to_label_start[cur_qa] != q_to_label_start[cur_qb]:
                        true_y[cnt] = q_mat[cur_qa, cur_qb]
                        true_x[cnt] = r_mat[cur_ra, cur_rb]
                        true_qa[cnt] = cur_qa
                        true_qb[cnt] = cur_qb
                        true_ra[cnt] = cur_ra
                        true_rb[cnt] = cur_rb
                        cnt += 1

            test_set_xy(mv_x, mv_y, mv_qa, mv_qb, mv_ra, mv_rb, cur_k_q_idx,
                        r_mat, q_mat, q_to_r_idx, q_to_label_start)

            self.assertTrue(np.allclose(true_x, mv_x))
            self.assertTrue(np.allclose(true_y, mv_y))
            self.assertTrue(np.all(mv_qa == true_qa))
            self.assertTrue(np.all(mv_qb == true_qb))
            self.assertTrue(np.all(mv_ra == true_ra))
            self.assertTrue(np.all(mv_rb == true_rb))


def get_test_matrix(n):
    mat = np.zeros((n, n), dtype=np.float64)
    mat_triu = np.triu_indices_from(mat, 1)
    mat[mat_triu] = np.random.random(mat_triu[0].size)
    return mat + mat.T


if __name__ == '__main__':
    unittest.main()
