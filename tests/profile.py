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


from deltapd.file import MetadataFile
from deltapd.main import run

if __name__ == '__main__':
    path_r_tree = '/Users/aaron/tmp/ar122_r95.tree'
    path_q_tree = '/Users/aaron/tmp/4_arc_gtr_f_g_trim.treefile'
    metadata_path = '/Users/aaron/tmp/ar122_metadata_r95.tsv'
    out_dir = '/tmp/deltapd2'
    msa_file = '/Users/aaron/tmp/ar122_ssu_r95.fna'

    inf_thresh = 2
    dif_thresh = 0.1
    k = 40
    cpus = 1

    metadata = MetadataFile(metadata_path)
    metadata.add_contig_lens(msa_file)
    run(path_r_tree, path_q_tree, metadata, inf_thresh, dif_thresh, out_dir, k, cpus)
