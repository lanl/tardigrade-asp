.. _user_manual:

###########
User Manual
###########

***********
Quick Start
***********

This stub repo contains hooks for writing Abaqus :cite:`ABAQUS2019` subroutines, like those found in the `Abaqus UMAT
documentation`_, and a template UMAT c++ interface. However, this template repository does not yet have a meaningful c++
constitutive model to be the subject of a user manual.

Until this template repository includes a `CMake`_ ``--install`` definition, user's are referred to the :ref:`build`
section of the :ref:`devops_manual` for build instructions.

The template UMAT can be used after build with the following Abaqus options

.. code:: bash

   $ abaqus -job <my_input_file> -user relative/path/to/cpp_stub/build/src/cpp/umat.o

Until the template repository and all upstream c++ libraries are built as shared library objects it is recommended that
the subroutines are left in the project build directory. However, it is possible to copy the shared library files to any
other directory provided the upstream projects ``{error,vector,stress,solver,constitutive}_tools`` are present in the
build directory, e.g. ``cpp_stub/build/_deps/{error,vector,stress,solver,constitutive}_tools-build/``.

.. code:: bash

   $ pwd
   /path/to/my/abaqus/job
   $ cp /path/to/cpp_stub/build/src/cpp/{umat.o,libcpp_stub.so} .
   $ abaqus -job <my_input_file> -user umat.o

******************************
Input File Material Definition
******************************

cpp_stub requires 2 material constants and 2 state variables. The c++ cpp_stub interface, material constants, and state
variables are described in the :ref:`sphinx_api`. The fixed expectations for the abaqus interface are defined in the
"Variables" section of the :ref:`sphinx_api` for :ref:`cpp_stub_source`. A complete discussion about the constants and their
meaning is not included here. Instead users are directed to calibrated material parameters found in cpp_stub entries in the
`Granta/MIMS`_ `Material Database`_ :cite:`MIMS`. Material parameter calibration sets should be availble for download
with the correct Abaqus input file formatting from MIMS.

The cpp_stub project contains abaqus integration tests for the cpp_stub abaqus interface. These tests perform actual abaqus
simulations using the same dummy parameters used for unit and integration testing of the cpp_stub c++ code. The cpp_stub
Abaqus input files used for integration testing can be found in the cpp_stub source code repository with the following bash
command

.. code:: bash

   $ pwd
   /path/to/my/cpp_stub
   $ find . -path ./build -prune -false -o -name "*.inp"
   ./src/abaqus/single_element_c3d8.inp

The material definition from an integration test input file is included below for reference

.. warning::

   The material constants used in this example material definition are *NOT* calibrated for any real material data.
   For calibrated material parameters, see the cpp_stub entry for materials found in the `Granta/MIMS`_ `Material
   Database`_ :cite:`MIMS`.

.. literalinclude:: ../src/abaqus/single_element_c3d8.inp
   :linenos:
   :lines: 42-50
