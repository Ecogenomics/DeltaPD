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

import re
from typing import Set, Dict, List


class MetadataFile(object):

    def __init__(self, path, qry_sep):
        self.path = path
        self.data = self.read()
        self.contig_lens = dict()
        self.qry_sep = qry_sep

    def get_rep_label_dict(self):
        return {k: v['gtdb_genome_representative'] for k, v in self.data.items()}

    def add_contig_lens(self, path):
        with open(path) as fh:
            re_hits = re.findall(r'>(.+)' + self.qry_sep + r'(.+) d__.+contig_len=(\d+)]', fh.read())
        for gid, gene_id, contig_len in re_hits:
            self.contig_lens[f'{gid}{self.qry_sep}{gene_id.replace("#", "_")}'] = int(contig_len)

    def get_contig_len(self, qry_label):
        return self.contig_lens[qry_label]

    def read(self):
        out = dict()
        with open(self.path, encoding='utf-8') as fh:
            indices = fh.readline().strip().split('\t')
            for line in fh:
                cols = line.strip().split('\t')
                out[cols[0]] = dict()
                for cur_idx, cur_val in zip(indices[1:], cols[1:]):

                    if cur_idx == 'gtdb_taxonomy':
                        new_val = cur_val.split(';')
                    else:
                        new_val = cur_val

                    out[cols[0]][cur_idx] = new_val
        return out

    def unique_ranks(self) -> Set[str]:
        out = set()
        for cur_val in self.data.values():
            out = out.union({x for x in cur_val['gtdb_taxonomy']})
        return out

    def sub_ranks(self) -> Dict[str, Set[str]]:
        out = dict()
        for cur_val in self.data.values():
            for i in range(len(cur_val['gtdb_taxonomy']) - 1):
                top_rank = cur_val['gtdb_taxonomy'][i]
                sub_rank = cur_val['gtdb_taxonomy'][i + 1]
                if top_rank not in out:
                    out[top_rank] = set()
                out[top_rank].add(sub_rank)
        return out

    def top_ranks(self) -> Dict[str, Set[str]]:
        out = dict()
        for cur_val in self.data.values():
            for i in range(len(cur_val['gtdb_taxonomy']) - 1):
                top_rank = cur_val['gtdb_taxonomy'][i]
                sub_rank = cur_val['gtdb_taxonomy'][i + 1]
                if sub_rank not in out:
                    out[sub_rank] = set()
                out[sub_rank].add(top_rank)
        return out

    def get_rep(self, gid):
        return self.data[gid]['gtdb_genome_representative']

    def get_gid_taxonomy(self, gid: str) -> List[str]:
        return self.data[gid]['gtdb_taxonomy']
