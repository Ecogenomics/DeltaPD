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

import os

from deltapd import __title__, __version__
from deltapd.cli import deltapd_parser, main_parser
from deltapd.exceptions import DeltaPDExit
from deltapd.main import run


def main():
    with deltapd_parser(main_parser(), __title__, __version__) as args:

        # Validate the arguments.
        if not os.path.isfile(args.r_tree):
            raise DeltaPDExit(f'The reference tree does not exist: {args.r_tree}')
        if not os.path.isfile(args.q_tree):
            raise DeltaPDExit(f'The query tree does not exist: {args.q_tree}')
        if not os.path.isfile(args.metadata):
            raise DeltaPDExit(f'The metadata file does not exist: {args.q_tree}')
        if not os.path.isfile(args.msa_file):
            raise DeltaPDExit(f'The msa file does not exist: {args.q_tree}')
        if args.influence_thresh < 0:
            raise DeltaPDExit(f'The influence threshold value must be >= 0: {args.influence_thresh}')
        if args.k <= 0:
            raise DeltaPDExit(f'The value of k must be greater than 0: {args.influence_thresh}')

        # Run the program.
        run(args.r_tree, args.q_tree, args.metadata, args.influence_thresh, args.diff_thresh, args.out_dir, args.k,
            args.plot, args.ete3_scale, args.max_taxa, args.qry_sep, args.msa_file, cpus=max(args.cpus, 1))


if __name__ == "__main__":
    main()
