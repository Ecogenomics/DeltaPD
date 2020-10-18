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

import hashlib
import os


def assert_path_exists(path):
    """Create the specified path and raise an exception if it can't be created.

    Parameters
    ----------
    path : str
        The path to the directory that should be created if not exists.

    Raises
    ------
    OSError
        If the directory wasn't created.
    """
    try:
        if not os.path.isdir(path):
            os.makedirs(path)
            if not os.path.isdir(path):
                raise OSError(f'Unable to create the path: {path}')
    except Exception:
        raise OSError


def sha256_file(path):
    """Compute the SHA256 of a file.

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    str
        The SHA256 sum of the file.
    """
    if not os.path.isfile(path):
        raise OSError(f'The file does not exist: {path}')
    sha256 = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            data = fh.read(65536)
            if not data:
                break
            sha256.update(data)
    return sha256.hexdigest()
