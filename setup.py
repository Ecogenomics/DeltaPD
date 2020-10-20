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
import platform
import re
from distutils.command.build import build as build_orig

from setuptools import setup, find_packages, Extension


def read_meta():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'deltapd/__init__.py')
    with open(path) as fh:
        hits = re.findall(r'__(\w+)__ ?= ?["\'](.+)["\']\n', fh.read())
    return {k: v for k, v in hits}


def readme():
    with open('README.md') as f:
        return f.read()


class build(build_orig):

    def finalize_options(self):
        super().finalize_options()
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        for extension in self.distribution.ext_modules:
            extension.include_dirs.append(numpy.get_include())
        from Cython.Build import cythonize
        self.distribution.ext_modules = cythonize(self.distribution.ext_modules,
                                                  language_level=3)


# Environment specific compiler setttings.
compile_extra_args = ['-O3', '-ffast-math', '-march=native']
link_extra_args = list()
if platform.system() == "Windows":
    pass
elif platform.system() == "Darwin":
    compile_extra_args.extend(['-std=c++11', "-mmacosx-version-min=10.9", '-Xpreprocessor', '-fopenmp'], )
    link_extra_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.9", '-Xpreprocessor', '-fopenmp'])
else:
    compile_extra_args.extend(['-fopenmp'])
    link_extra_args.extend(['-fopenmp'])

# Select Cython modules to compile.
ext_modules = [Extension('deltapd.model', ['deltapd/model.pyx'],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args
                         ),
               Extension('deltapd.model_test', ['deltapd/model_test.pyx'],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args
                         ),
               Extension('deltapd.util', ['deltapd/util.pyx'],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args
                         ),
               Extension('deltapd.util_test', ['deltapd/util_test.pyx'],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args
                         )
               ]

meta = read_meta()
setup(name=meta['title'],
      version=meta['version'],
      description=meta['description'],
      long_description=readme(),
      long_description_content_type='text/markdown',
      author=meta['author'],
      author_email=meta['author_email'],
      url=meta['url'],
      license=meta['license'],
      project_urls={
          'Bug Tracker': meta['bug_url'],
          'Documentation': meta['doc_url'],
          'Source Code': meta['src_url'],
      },
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords=meta['keywords'],
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'deltapd = deltapd.__main__:main'
          ]
      },
      install_requires=['numpy', 'phylodm', 'tqdm', 'ete3', 'dendropy', 'matplotlib', 'jinja2'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])],
      ext_modules=ext_modules,
      cmdclass={'build': build}
      )
