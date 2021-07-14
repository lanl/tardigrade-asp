###################
C++ Stub repository
###################

*******************
Project Description
*******************

A stub repository for C++ development projects in W-13.

This repository will contain the necessary setup files to integrate C++ doc
strings, `CMake <https://cmake.org/cmake/help/v3.14/>`_,
`Doxygen <https://www.doxygen.nl/manual/docblocks.html>`_,
`Sphinx <https://www.sphinx-doc.org/en/master/>`_, and
`Breathe <https://breathe.readthedocs.io/en/latest/>`_ for a complete build
system with integrated documentation. It will also include the necessary hooks
to commonly used C++ libraries for constitutive modeling. This stub repository
also includes template hooks for integrating C++ code as Abaqus subroutines.

    **NOTE**: you can use this repo as a stub for fortran projects as well! A
    step-by-step may never happen because c++ is the future of W-13 subroutines.
    For now, refer to the following references:
   
    * `CMake for Fortran example CMakeLists.txt <https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/languages/fortran/ForFortranExample>`_
    * `CMake documentation starting point <https://cmake.org/cmake/help/v3.14/module/CheckFortranSourceRuns.html>`_
    * `Stack Overflow thread <https://stackoverflow.com/questions/12705562/using-cmake-with-fortran>`_
    * `Doxygen comments for Fortran <https://www.doxygen.nl/manual/docblocks.html#fortranblocks>`_

********************************************
Setting up a new project from this stub repo
********************************************

This section has a corresponding section in the formatted and hyperlinked html
documentation. The html documentation also includes setup instructions for
Jenkins Continuous Integration (CI) setup. Build and view instructions are
included in following sections, separated from project setup by horizontal bars.

Prerequisites
=============

Set up your profile on Bitbucket with ssh keys. You can follow the instructions
on the W-13 Confluence page
`Gitting started with W-13 Bitbucket <https://xcp-confluence.lanl.gov/display/GIT/Gitting+Started+W-13%27s+Git+Server>`_

Clone cpp\_stub into a local repository
=======================================

1. Navigate to the `upstream cpp\_stub repo <https://xcp-stash.lanl.gov/projects/MM/repos/cpp_stub/browse>`_

2. Copy the ssh url from the Bitbucket "Clone" button on the Bitbucket
   repository web page. It should look like the following:

   .. code-block:: bash

      ssh://git@xcp-stash.lanl.gov:7999/mm/cpp_stub.git

3. Navigate to your preferred repository directory in a terminal or use the
   example commands below

   .. code-block:: bash

      $ ssh -X sstelmo.lanl.gov
      $ pwd
      /home/<moniker>
      $ mkdir -p /projects/$USER/w13repos
      $ cd /projects/$USER/w13repos
      $ pwd
      /projects/<moniker>/w13repos

4. Clone the stub repository using the url copied in step 2.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos
      $ git clone ssh://git@xcp-stash.lanl.gov:7999/mm/cpp_stub.git

5. Rename the local repository directory for your project.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos
      $ ls cpp_stub -d
      cpp_stub
      $ mv cpp_stub my_project
      $ ls cpp_stub -d
      ls: cannot access 'cpp_stub': No such file or directory
      $ mv cpp_stub my_project
      my_project

6. Change to your project's repository directory

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos
      $ cd my_project
      $ pwd
      /projects/<moniker>/w13repos/my_project

Create a new repository on Bitbucket
====================================

    Note: These notes are a text copy of a variation on the
    `New Bitbucket Repo Guide <https://simulia.lanl.gov/ECMF-D/devops_guide.html#new-bitbucket-repo-guide>`_
    which can also be found in the
    `W-13 DevOps Manual <https://xcp-confluence.lanl.gov/display/COM/W-13+DevOps>`_

1. Navigate to the W-13 `Material Models <https://xcp-stash.lanl.gov/projects/MM>`_ Bitbucket project.

2. Create a new repository by clicking on the "+" sign in the upper left corner.

3. Enter a name for your project and click "Create repository"

4. Follow the "My code is already tracked by Git" instructions.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git remote set-url origin ssh://git@xcp-stash.lanl.gov:7999/mm/my_project.git
      $ git push -u origin --all
      $ git push origin --tags

5. Refresh the Bitbucket webpage and verify that the repository code was pushed
   correctly. You should see a list of source files and this Bitbucket parsed
   ``README.md`` displayed. You can also select the drop down branch menu to
   view a "master" and "dev" branch.

Update settings for your repository
===================================

Bitbucket repositories in the
`Material Models <https://xcp-stash.lanl.gov/projects/MM>`_ project inherit permissions and
settings from that project. This included read permission for the
``w13bitbucket``
`UNIX group <https://xcp-confluence.lanl.gov/pages/viewpage.action?pageId=150929410&searchId=Y1EVB37UN>`_.
For most developers, these inherited repository settings are appropriate and
only a small number of settings must be updated.

1. Click on the gear icon in the lower left sidebar.

2. From the "Repository details" landing page, update the default branch from
   "master" to "dev".

3. From the "Repository permissions" tab you can add additional permissions by
   user and UNIX group.

4. From the "Default reviewers" tab you can add yourself and any project
   co-owners as default reviewers.

Fork the upstream repository
============================

In the
`Forking Workflow <https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow>`_
the repository you just created in the
`Material Models <https://xcp-stash.lanl.gov/projects/MM>`_ project is called
the "upstream" repository. Throughout older W-13 documentation this may also be
called the "official" repository.

Bitbucket repositories that inherit permissions from W-13 projects use the
`Forking Workflow <https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow>`_
and limit permissions for pushing changes to the upstream repository. Now that
branches exist on this repository, no one will be able to push directly to
*existing branches* of the upstream respository.

1. Click the fork button in the left hand sidebar just above the gear icon.

2. Click "fork this repository" button.

3. Verify

   a. The "Project" points to your personal project space in Bitbucket. It will
      probably show your full name as it appears in the phonebook.

   b. The "Name" matches your project name

   c. The "Enable fork syncing" checkbox is checked

4. Click the "Fork Repository" button. You should land on a familiar looking
   source code view of your repository, but now located in your personal project
   space.

This repository is referred to as the "fork" or "remote" repository throughout
W-13 DevOps documentation.

Update the remote url in your local repository
==============================================

The final repo setup step is to update the remote url of the local clone of
``my_project``.  We will return to the terminal session.

1. Copy the url of your "remote" repository from the Bitbucket webpage. It
   should look like:

   .. code-block:: bash

      ssh://git@xcp-stash.lanl.gov:7999/~<moniker>/my_project.git

2. Return to your terminal session and update the remote repository for the
   final time.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git remote set-url origin ssh://git@xcp-stash.lanl.gov:7999/~<moniker>/my_project.git
      $ git push -u origin --all
      $ git push origin --tags

Update project name throughout repository
=========================================

    Note: the remaining steps are a truncated version of the W-13 Git project
    `contribution guide <https://simulia.lanl.gov/ECMF-D/contribution_guide.html>`_
    which can also be found in the
    `W-13 DevOps Manual <https://xcp-confluence.lanl.gov/display/COM/W-13+DevOps>`_.
    Critically, these steps will omit the Jira task creation and Bitbucket
    Pull-Request (PR) steps.  The Bitbucket PR steps may be reproduced using the
    contribution guide, but your project will have to create a Jira project prior to
    integrating the Jira workflow. Contact the xcp devops team
    <devops-help@lanl.gov> to create a Jira project. You can email the W-13 DevOps
    team <w13devops@lanl.gov> for notes about setup.

1. Create a feature branch for your project name updates

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git checkout -b feature/project-name-updates
      $ git branch
        dev
      * feature/project-name-updates
        master

2. Search for all instances of ``cpp_stub``. The list of occurrences will look
   quite long, but we can search and replace with ``sed`` to avoid manual file
   edits.

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project

      # Recursive, case-insensitive search and count occurrences
      $ grep -ri cpp_stub . --exclude-dir={build,.git} | wc -l
      57

      # Recursive, case-insensitive search and display
      $ grep -ri cpp_stub . --exclude-dir={build,.git}
      ...

      # Clean list of files with project name
      $ grep -ri cpp_stub . --exclude-dir={build,.git} -l
      ./CMakeLists.txt
      ./docs/api.rst
      ./docs/devops.rst
      ./README.md
      ./set_vars.sh
      ./src/cpp/cpp_stub.cpp
      ./src/cpp/cpp_stub.h
      ./src/cpp/tests/test_cpp_stub.cpp

3. Search and replace from command line

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project

      # Replace lower case occurrences in place
      $ sed -i 's/cpp_stub/my_project/g' $(grep -ri cpp_stub . --exclude-dir={build,.git} -l)
      $ grep -ri cpp_stub . --exclude-dir={build,.git} -l
      ./src/cpp/cpp_stub.h

      # Replace upper case occurrences in place
      $ sed -i 's/CPP_STUB/MY_PROJECT/g' $(grep -ri cpp_stub . --exclude-dir={build,.git} -l)

4. Verify no more occurrences of project name ``cpp_stub``

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ grep -ri cpp_stub . --exclude-dir={build,.git} | wc -l
      0
      $ grep -ri cpp_stub . --exclude-dir={build,.git}
      # no stdout to terminal because no files found
      $ grep -ri cpp_stub . --exclude-dir={build,.git} -l
      # no stdout to terminal because no files found

5. Search and replace camelcase project name occurrences, e.g. ``cppStub``.

   .. code-block:: bash

      $ grep -r cppStub . --exclude-dir={build,.git}
      ...
      $ sed -i 's/cppStub/myProject/g' $(grep -r cppStub . --exclude-dir={build,.git} -l)
      $ grep -r cppStub . --exclude-dir={build,.git} -l
      # no stdout to terminal because no files found

6. Find files containing the project in their file name

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ find . -type d \( -name .git -o -name build \) -prune -false -o -name "*cpp_stub*"
      ./src/cpp/cpp_stub.cpp
      ./src/cpp/cpp_stub.h
      ./src/cpp/tests/test_cpp_stub.cpp

7. Rename files after current project

   .. code-block:: bash

      $ rename 's/cpp_stub/myproject/' $(find . -type d \( -name .git -o -name build \) -prune -false -o -name "*cpp_stub*")

8. Commit and push your changes to your "remote" or "fork" repository

   .. code-block:: bash

      $ pwd
      /projects/<moniker>/w13repos/my_project
      # Add tracked files and message
      $ git commit -a -m "FEAT: replace cpp_stub with my_project through repository"
      $ git push origin feature/project-name-updates

You can also perform some cleanup in ``README.md`` to remove this walk-through.

From here, the W-13 best practice workflow would return to the Bitbucket webpage
and submit a Pull-Request from the ``feature/project-name-updates`` branch of
``\<moniker\>/my_project`` repository (a.k.a. fork or remote) to the ``dev`` branch
of ``Material Models/my_project`` repository (a.k.a. upstream or official).

After updating your project by merging to the upstream repository, the fork
syncing feature of Bitbucket will automatically update any identically named
branches in your fork repository. Best practices suggest you should limit the
upstream repository branches to clean ``dev`` and ``master`` branches and
*NEVER* develop directly on the ``dev`` and ``master`` branches of your fork
repository. Limit development work to ``feature/thing`` type branches on your
fork/remote repo and frequently commit changes and push from the local feature
branch back to the fork/remote repo.

Happy hacking!

************
Gitlab CI/CD
************

    Pending...

************
Dependencies
************

Compilers
=========

* c++11 compiler (listed version number has been tested at some point)

  * g++ >= GNU 4.8.5

Executables
===========

* `CMake <https://cmake.org/cmake/help/v3.14/>`_ >= 3.14
* `Doxygen <https://www.doxygen.nl/manual/docblocks.html>`_ >= 1.8.5
* `LaTeX <https://www.latex-project.org/help/documentation/>`_ >= 2017

Conda Environment
=================

For convenience, the minimal Python environment requirements for the
documentation build are included in ``configuration_files/environment.yaml``.
This file was created from the `pipreqs <https://github.com/bndr/pipreqs>`_
command line tool and Sphinx configuration inspection, e.g. the extension
packages.

.. code-block:: bash

   $ pwd
   path/to/cpp_stub/
   $ pipreqs --use-local --print --no-pin .

A minimal anaconda environment for building the documentation can be created
from an existing anaconda installation with the following commands.

.. code-block:: bash

   $ conda env create --file configuration_files/environment.yaml

You can learn more about Anaconda Python environment creation and management in
the
`Anaconda Documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_

C++ Libraries
=============

    **NOTE: Non-admin installations for Eigen and Boost are no longer required.** This project is built and deployed
    against C++ libraries managed in Conda. See the Conda environment file and README discussion for non-admin environment
    management.

* `Eigen <https://eigen.tuxfamily.org/dox/>`_ >= 3.3.7
* `BOOST <https://www.boost.org/doc/libs/1_53_0/>`_ >= 1.59.0
* error\_tools: https://re-git.lanl.gov/aea/material-models/error_tools
* vector\_tools: https://re-git.lanl.gov/aea/material-models/vector_tools
* abaqus\_tools: https://re-git.lanl.gov/aea/material-models/abaqus_tools
* constitutive\_tools: https://re-git.lanl.gov/aea/material-models/constitutive_tools
* stress\_tools: https://re-git.lanl.gov/aea/material-models/stress_tools
* solver\_tools: https://re-git.lanl.gov/aea/material-models/solver_tools

If not found on the current system or active Conda environment, all of the
``*_tools`` libraries are pulled from their git repos by branch name and built
with their respective cmake files as part of the cmake build for this project.

**************
Build and Test
**************

This project is built with `CMake <https://cmake.org/cmake/help/v3.14/>`_ and uses
`Sphinx <https://www.sphinx-doc.org/en/master/>`_ to build the documentation with
`Doxygen <https://www.doxygen.nl/manual/docblocks.html>`_ +
`Breathe <https://breathe.readthedocs.io/en/latest/>`_ for the c++ API.

Build on sstelmo
================

1) Activate the correct python environment

   .. code-block:: bash

      $ module load python/2020.07-python-3.8
      $ sv3r

2) Create a build directory

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/

      $ mkdir build
      $ cd build

3) Configure ``cmake3``

       This step only needs to be performed once unless you need to specify a
       new CMake configuration for a re-build. Most command line arguments and
       environment variables are stored in the CMake cache. Anything found in cache
       will not be re-configured unless you remove the cache file or clobber the build
       directory.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake3 ..

4) Build various portions of the project

       Most of the project will re-build only as necessary after source updates. Some portions of the documentation
       require a ``make clean`` after documentation source file updates to force a re-build.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Build everything
      $ cmake3 --build .

      # Build only the c++ primary libraries
      $ cmake3 --build src/cpp

5) Locate build files

       The build directory structure may change between version releases. Developers and users are encouraged to become
       familiar with the bash ``find``, ``grep``, and ``tree`` commands to locate build files.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # find c++ libraries and ignore intermediate files with similar extensions
      $ find . \( -name "*.o" -o -name "*.so" -o -name "*.a" \) | grep -vE "\.cpp\."

6) Clean build directory to force a re-build

       **HEALTH WARNING**
      
       The abaqus input files and bash scripts used for integration testing are
       built with the
       `CMake add_custom_target <https://cmake.org/cmake/help/latest/command/add_custom_target.html>`_
       feature. Consequently, the integration test target is *always considered
       out of date*. The integration test target copies all registered input files
       and the integration test bash script from source to build directory. This
       means the file copy operation is always performed when the integration test
       target is requested in the cmake build command, e.g. ``cmake --build .`` or
       ``cmake --build src/abaqus/tests``. This operation is computationally
       inexpensive with respect to building the ``cpp_stub`` source code.
      
       Input files are registered in the ``src/abaqus/tests/CMakeLists.txt`` file
       under the ``ABAQUS_INPUT_FILES`` CMake variable.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      $ make clean

Test on sstelmo
===============

4) Build tests of the project

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Build c++ tests
      $ cmake3 --build src/cpp/tests

      # Build Abaqus integration tests
      $ cmake3 --build src/abaqus/tests

5) Run the tests

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Run ctest
      $ ctest

      # Results print to screen
      # View details of most recent test execution including failure messages
      $ less Testing/Temporary/LastTest.log

Convenience build wrappers
==========================

Two build scripts have been created for convenience, ``new_build.sh`` and
``build_docs.sh``. The first will build everything including the library binary,
the test binary, and the documentation. This is the same build script used by
``jenkins_build.sh`` for CI builds and testing. The ``build_docs.sh`` script
only builds the documentation. Both build scripts clobber existing build
directories, reset any bash environment variables, and run the cmake
configuration from scratch.

2) Build everything and run tests

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/

      # Just perform the build (pick one)
      $ ./new_build.sh <cmake build type>
      $ ./new_build.sh None
      $ ./new_build.sh Release

      # Perform tests from PWD
      $ ./build/src/cpp/tests/test_cpp_stub

      # Build and perform tests
      $ ./jenkins_build.sh

3) View test results

   .. code-block:: bash

      # As built directly to PWD
      $ cat results.tex

      # As built by jenkins_build.sh
      $ cat build/src/cpp/tests/*_results.tex
      $ cat *results.tex

4) Display docs

   .. code-block:: bash

      # Sphinx
      $ firefox build/docs/sphinx/html/index.html &

      # Doxygen
      $ firefox build/docs/doxygen/html/index.html &

Building the documentation
==========================

    **HEALTH WARNING**
   
    **API Health Note**: The sphinx API docs are a work-in-progress. The doxygen
    API is much more useful.

The documentation can be built with ``build_docs.sh``. The steps used in that
shell script are repeated here.

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/
      $ mkdir build/
      $ cd build/

3) Run cmake3 configuration

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ cmake3 ..

4) Build the docs

   .. code-block:: bash

      $ cmake3 --build docs

5) Documentation builds to:

   .. code-block:: bash

      cpp_stub/build/docs/sphinx/html/index.html

6) Display docs

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/sphinx/html/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/doxygen/html/index.html &

*******************
Install the library
*******************

Build the entire before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake3 --build .

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake --install . --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install . --prefix /home/$USER/.local

      # Example install to conda environment
      $ conda active my_env
      $ cmake --install . --prefix ${CONDA_DEFAULT_ENV}

      # Example install to W-13 CI/CD conda environment performed by CI/CD institutional account
      $ cmake --install . --prefix /projects/python/release

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

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

`Sphinx <https://www.sphinx-doc.org/en/master/>`_ reads in docstrings and other
special portions of the code as reStructured text. Developers should follow
styles in this `Sphinx style guide
<https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#>`_.

Style Guide
===========

This project does not yet have a full style guide. Generally, wherever a style
can't be inferred from surrounding code this project falls back to
`PEP-8 <https://www.python.org/dev/peps/pep-0008/>`_-like styles. There are two
notable exceptions to the notional PEP-8 fall back:

1. `Doxygen <https://www.doxygen.nl/manual/docblocks.html>`_ style docstrings are
   required for automated, API from source documentation.
2. This project prefers expansive whitespace surrounding parentheses, braces, and
   brackets.

   * No leading space between a function and the argument list.
   * One space following an open paranthesis ``(``, brace ``{``, or bracket
     ``[``
   * One space leading a close paranthesis ``)``, brace ``}``, or bracket ``]``

An example of the whitespace style:

.. code-block:: bash

   my_function( arg1, { arg2, arg3 }, arg4 );

The following ``sed`` commands may be useful for updating white space, but must
be used with care. The developer is recommended to use a unique git commit
between each command with a corresponding review of the changes and a unit test
run.

* Trailing space for open paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([({[]\)\([^ ]\)/\1 \2/g' <list of files to update>

* Leading space for close paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([^ ]\)\([)}\]]\)/\1 \2/g' <list of files to update>

* White space between adjacent paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([)}\]]\)\([)}\]]\)/\1 \2/g' <list of files to update>
