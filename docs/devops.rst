#############
DevOps Manual
#############

.. note::

   This is a duplicate of the content found in the README with a manual
   conversion to ReST syntax.

************
Dependencies
************

Compilers
=========

* c++11 compiler

Executables
===========

* CMake >= 3.14
* Doxygen >= 1.8.5
* LaTeX >= 2017

Python Modules (for documentation)
==================================

* Sphinx >= 3.0.4
* Breathe >= 4.18.1
* sphinx\_rtd\_theme >= 0.4.3

For convenience, the minimal Python environment requirements for the
documentation build are included in ``environment.yaml`` and
``requirements.txt``. A minimal anaconda environment for building the
documentation can be created from an existing anaconda installation with the
following commands.

.. code:: bash

   $ conda env create --file environment.yaml

You can learn more about Anaconda Python environment creation and management in
the [Anaconda
Documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

C++ Libraries
=============

* eigen >= 3.3.7
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
* vector\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/vector_tools
* stress\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/stress_tools
* solver\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/solver_tools
* constitutive\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/consitutive_tools

Constitutive Tools
------------------

All of the ``{error,vector,stress,solver,constitutive}_tools`` libraries are
pulled from their git repos by branch name and built with their respective cmake
files as part of the cmake build for this project.

Eigen
-----

https://gitlab.com/libeigen/eigen

Eigen must be "installed" following the ``eigen/INSTALL`` instructions. The
Eigen dependence is easiest to resolve if eigen is installed in the default
install directory.  However, if you don't have admin privileges, you can also
insall Eigen to your home directory in ``$HOME/.local/include`` or
``$HOME/include``.

Non-admin Eigen install
-----------------------

[Reference](https://unix.stackexchange.com/questions/36871/where-should-a-local-executable-be-placed)

.. code:: bash

   # Create personal include file directory
   $ pwd
   /home/$USER
   $ mkdir -p .local/include
   # Move to repository directory
   $ cd /preferred/path/to/repos
   # Example
   $ pwd
   /projects/$USER/w13repos
   # Clone eigen
   $ git clone https://gitlab.com/libeigen/eigen.git
   $ cd eigen
   $ git checkout 3.3.7
   # Create build directory
   $ mkdir build
   $ cd build
   # OPTIONAL. Set c++ compiler separate from system default
   $ export CXX=$(command -v g++)
   # Build eigen
   $ cmake3 .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
   $ make install

---

---

**************
Build and Test
**************

This repository is built with CMake and uses Doxygen + Sphinx + Breathe to build
the documentation.

> **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
> API is much more useful.

Two build scripts have been created for convenience, ``new_build.sh`` and
``build_docs.sh``. The first will build everything including the library binary,
the test binary, and the documentation. This is the same build script used by
``jenkins_build.sh`` for CI builds and testing. The ``build_docs.sh`` script
only builds the documentation. Both build scripts clobber existing build
directories, reset any bash environment variables, and run the cmake
configuration from scratch.

Build on sstelmo
================

1) Activate the correct python environment

   .. code:: bash

      $ module load python/2019.10-python-3.7
      $ sv3r

2) Build everything and run tests

   .. code:: bash

      $ pwd
      /path/to/cpp_stub/

      # Just perform the build (pick one)
      $ ./new_build.sh <cxx compiler>
      $ ./new_build.sh c++
      $ ./new_build.sh g++
      $ ./new_build.sh icpc

      # Perform tests from PWD
      $ ./build/src/cpp/tests/test_cpp_stub

      # Build and perform tests
      $ ./jenkins_build.sh

3) View test results

   .. code:: bash

      # As built directly to PWD
      $ cat results.tex

      # As built by jenkins_build.sh
      $ cat build/src/cpp/tests/*_results.tex
      $ cat *results.tex

4) Display docs

   .. code:: bash

      # Sphinx
      $ firefox build/docs/sphinx/index.html &

      # Doxygen
      $ firefox build/docs/doxygen/html/index.html &

Building the documentation
==========================

The documentation can be built with ``build_docs.sh``. The steps used in that
shell script are repeated here.

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code:: bash

      $ pwd
      /path/to/cpp_stub/
      $ mkdir build/
      $ cd build/

3) Run cmake3 configuration

   .. code:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ cmake3 ..

4) Build the docs

   .. code:: bash

      $ cmake3 --build docs

5) Documentation builds to:

   .. code:: bash

      cpp_stub/build/docs/sphinx/index.html

6) Display docs

   .. code:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/sphinx/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/doxygen/html/index.html &

---

---

***********************
Contribution Guidelines
***********************

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code:: bash

   git commit -m "FEAT: short intent of new feature"
   git commit -m "BUG: fixes nasty bug"
   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* bugfix/\<description>
* feature/\<description>
* release/\<description>

reStructured Text
=================

Sphinx reads in docstrings and other special portions of the code as
reStructured text. Developers should follow styles in this [Sphinx style
guide](https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#).
