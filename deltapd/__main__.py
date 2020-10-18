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

from deltapd import __title__, __version__, __url__
from deltapd.argparse import CustomArgParser, deltapd_parser
from deltapd.exceptions import DeltaPDExit
from deltapd.main import run


def main_parser():
    parser = CustomArgParser(prog=__title__, ver=__version__, url=__url__)

    req_args = parser.add_argument_group('I/O arguments (required)')
    req_args.add_argument('--r_tree', help='path to the reference tree')
    req_args.add_argument('--q_tree', help='path to the query tree')
    req_args.add_argument('--metadata', help='path to the metadata file')
    req_args.add_argument('--msa_file', help='path to the msa file used to infer the query tree')
    req_args.add_argument('--out_dir', help='path to output directory')

    opt_qry = parser.add_argument_group('Query tree arguments (optional)')
    opt_qry.add_argument('--max_taxa', type=int, default=1000,
                         help='if a ref taxon represents more than this number of qry taxa, ignore it')
    opt_qry.add_argument('--qry_sep', type=str, default='___',
                         help='query taxon separator in query tree, e.g. taxon___geneid')

    opt_out = parser.add_argument_group('Outlier arguments (optional)')
    opt_out.add_argument('--influence_thresh', help='outlier influence threshold value [0,inf)',
                         type=float, default=2)
    opt_out.add_argument('--diff_thresh', help='minimum change to base model to be considered an outlier',
                         type=float, default=0.1)
    opt_out.add_argument('--k', help='consider the query taxa represented by '
                                     'the ``k`` nearest neighbours for each representative taxon',
                         type=int, default=50)

    opt_plt = parser.add_argument_group('Plotting arguments (optional)')
    opt_plt.add_argument('--plot', help='generate outlier plots (slow)', action='store_true', default=False)
    opt_plt.add_argument('--ete3_scale', help='pixels per branch length unit', type=int, default=200)

    opt_args = parser.add_argument_group('Program arguments (optional)')
    opt_args.add_argument('--cpus', help='number of CPUs to use', type=int, default=1)
    opt_args.add_argument('--debug', help='output debugging information', action='store_true', default=False)

    return parser


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
