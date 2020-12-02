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
from glob import glob

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext

try:
    from Cython.Build import cythonize

    CYTHON = True
    EXT = 'pyx'
    LANG = 'c++'
except ImportError:
    CYTHON = False
    EXT = 'cpp'
    LANG = ''


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


# Environment-specific compiler settings.
compile_extra_args = ['-O3', '-ffast-math', '-march=native']
link_extra_args = list()
if platform.system() == "Windows":
    pass  # TODO: Test this case.
elif platform.system() == "Darwin":
    compile_extra_args.extend(['-std=c++11', "-mmacosx-version-min=10.9", '-Xpreprocessor', '-fopenmp'], )
    link_extra_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.9", '-Xpreprocessor', '-fopenmp'])
else:
    compile_extra_args.extend(['-fopenmp'])
    link_extra_args.extend(['-fopenmp'])

# Include the Cython modules (either to translate, or pre-translated)
ext_modules = [
    Extension(name='deltapd.model',
              sources=[f'deltapd/model.{EXT}'],
              language=LANG,
              extra_compile_args=compile_extra_args,
              extra_link_args=link_extra_args),
    Extension(name='deltapd.model_test',
              sources=[f'deltapd/model_test.{EXT}'],
              language=LANG,
              extra_compile_args=compile_extra_args,
              extra_link_args=link_extra_args),
    Extension(name='deltapd.util',
              sources=[f'deltapd/util.{EXT}'],
              language=LANG,
              extra_compile_args=compile_extra_args,
              extra_link_args=link_extra_args),
    Extension(name='deltapd.util_test',
              sources=[f'deltapd/util_test.{EXT}'],
              language=LANG,
              extra_compile_args=compile_extra_args,
              extra_link_args=link_extra_args)
]

if CYTHON:
    ext_modules = cythonize(ext_modules, compiler_directives={'language_level': '3'})


# Read the package information.
def read_meta():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'deltapd/__init__.py')
    with open(path) as fh:
        hits = re.findall(r'__(\w+)__ ?= ?["\'](.+)["\']\n', fh.read())
    return {k: v for k, v in hits}


# Package the README
def readme():
    with open('README.md') as f:
        return f.read()


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
      entry_points={
          'console_scripts': [
              'deltapd = deltapd.__main__:main'
          ]
      },
      install_requires=['numpy', 'phylodm', 'tqdm', 'ete3', 'dendropy', 'matplotlib', 'jinja2', 'seaborn'],
      python_requires='>=3.6',
      setup_requires=['numpy', 'cython'],
      include_package_data=True,
      package_data={
          'deltapd': ['*.pyx', '*.pxd']
      },
      data_files=[
          ("", ["LICENSE"]),
          ('deltapd/templates', glob('deltapd/templates/*.html')),
          ('deltapd/templates/assets/css', glob('deltapd/templates/assets/css/*.css')),
          ('deltapd/templates/assets/js', glob('deltapd/templates/assets/js/*.js')),
          ('deltapd/templates/static/img', glob('deltapd/templates/static/img/*.svg')),
          ('deltapd/templates/static/font/Roboto', ['deltapd/templates/static/font/Roboto/LICENSE.txt',
                                                    'deltapd/templates/static/font/Roboto/Roboto-Regular.ttf'])
      ],
      ext_modules=ext_modules,
      cmdclass={"build_ext": build_ext}
      )
