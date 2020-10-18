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

import argparse
import logging
import sys
from contextlib import contextmanager
from gettext import gettext

from deltapd.exceptions import DeltaPDExit
from deltapd.logging import colour, logger_setup


class CustomArgParser(argparse.ArgumentParser):
    """A custom argument parser for the `argparse` module."""

    def __init__(self, *args, **kwargs):
        self._version = kwargs.pop('ver') if 'ver' in kwargs else None
        self._url = kwargs.pop('url') if 'url' in kwargs else None
        super(CustomArgParser, self).__init__(*args, **kwargs)

    def print_usage(self, file=None):
        """Prints the top-level commands.

        Parameters
        ----------
        file : sys.stdout, sys.stderr
            The file to write the help message to.
        """
        if file is None:
            file = sys.stdout
        lines = ['  ' + colour(f'{self.prog} v{self._version}', fg='blue'),
                 f'  {self._url}',
                 '']
        file.write('\n'.join(lines))

    def print_help(self, file=None):
        """Prints the help message.

        Parameters
        ----------
        file : sys.stdout, sys.stderr
            The file to write the help message to.
        """
        if file is None:
            file = sys.stdout
        file.write(self.format_help()[0].upper() + self.format_help()[1:] + '\n')

    def exit(self, status=0, message=None):
        """Prints the result of a failed interaction with `argparse`.

        Parameters
        ----------
        status : int
            The code to exit with.
        message: str, optional
            The message to display on exit.
        """
        if message:
            sys.stderr.write(f"{colour(message, fg='red')}\n")
        sys.exit(status)

    def error(self, message):
        """Prints the result of an error with `argparse`.

        Parameters
        ----------
        message : str
            The error message to display.
        """
        self.print_usage(sys.stderr)
        args = {'prog': self.prog, 'message': message}
        self.exit(2, gettext('[%(prog)s] Error: %(message)s\n') % args)

    def print_version(self):
        """Display the current version of the program.

        Examples
        --------
        >>> x.print_version()
        title v1.2.3
        """
        sys.stdout.write(f"{colour(f'{self.prog} v{self._version}', fg='blue')}\n")


@contextmanager
def deltapd_parser(parser, title, version):
    """A context manager to safely handle user arguments and error reporting.

    Arguments
    ---------
    parser : argparse.CustomArgParser
        The arguments to parse.
    title : str
        The program title.
    version : str
        The program version.

    Examples
    --------
    >>> with deltapd_parser(parser, 'foo', '1.2.3') as args:
            pass

    """
    try:
        if len(sys.argv) == 1:
            parser.print_usage()
            sys.exit(1)
        elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
            parser.print_version()
            sys.exit(0)
        else:
            args = parser.parse_args()
            logger_setup(args.out_dir if hasattr(args, 'out_dir') else None,
                         f'{title}.log', title, version, False,
                         hasattr(args, 'debug') and args.debug)
            yield parser.parse_args()

    except KeyboardInterrupt:
        log = logging.getLogger('default')
        log.error('Controlled exit resulting from keyboard interrupt.')
        sys.exit(1)
    except DeltaPDExit as e:
        log = logging.getLogger('default')
        log.error(str(e))
        sys.exit(1)
