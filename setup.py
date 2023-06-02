# pylint: disable=invalid-name, undefined-all-variable, deprecated-module, consider-using-f-string, unspecified-encoding, exec-used, missing-module-docstring
import distutils.command.build as _build  # pylint: disable=deprecated-module
import os
import sys
from distutils import spawn
from distutils.sysconfig import get_python_lib
from setuptools import setup


def extend_build():
    """Define external build cmd for cmake."""
    class build(_build.build):
        """Self-define build."""

        def run(self):
            """Self-define run."""
            cwd = os.getcwd()
            if spawn.find_executable('cmake') is None:
                sys.stderr.write("CMake is required to build this package.\n")
                sys.exit(-1)
            _source_dir = os.path.split(__file__)[0]
            _build_dir = os.path.join(_source_dir, 'build_setup_py')
            _prefix = get_python_lib()
            try:
                cmake_configure_command = [
                    'cmake',
                    '-H{0}'.format(_source_dir),
                    '-B{0}'.format(_build_dir),
                    '-DCMAKE_INSTALL_PREFIX={0}'.format(_prefix),
                ]
                _generator = os.getenv('CMAKE_GENERATOR')
                if _generator is not None:
                    cmake_configure_command.append('-G{0}'.format(_generator))
                spawn.spawn(cmake_configure_command)
                spawn.spawn(
                    ['cmake', '--build', _build_dir, '--target', 'install'])
                os.chdir(cwd)
            except spawn.DistutilsExecError:
                sys.stderr.write("Error while building with CMake\n")
                sys.exit(-1)
            _build.build.run(self)

    return build


_here = os.path.abspath(os.path.dirname(__file__))

_this_package = 'src/pydisort'

version = {}
with open(os.path.join(_here, _this_package, 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='pydisort',
    version=version['__version__'],
    description='Testing pybind packaging.',
    author='Zoey Hu',
    author_email='zoey.zyhu@gmail.com',
    license='GPL',
    packages=['pydisort'],
    package_dir={'': 'src'},
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6'
    ],
    cmdclass={'build': extend_build()})
