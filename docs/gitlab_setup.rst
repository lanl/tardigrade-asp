########################
Bitbucket and repo setup
########################

.. include:: md2rst.txt

*************
Prerequisites
*************

1. Request access to `ASC RE Gitlab`_ through the `ASC RE Gitlab HPC Accounts`_
   page.
2. Set up your profile on `ASC RE Gitlab`_ with ssh keys. You can follow the
   instructions on the `ASC RE Gitlab User Documentation`_:

   * https://re-git.lanl.gov/help/ssh/README.md#generate-an-ssh-key-pair
   * https://re-git.lanl.gov/help/ssh/README.md#add-an-ssh-key-to-your-gitlab-account

***************************************
Clone cpp\_stub into a local repository
***************************************

1. Navigate to the `upstream repository`_.

2. Copy the ssh URL from the blue Gitlab "Clone" button on the
   `upstream repository`_ web page. The URL should look like the following:

   .. code:: bash
   
      ssh://git@re-git.lanl.gov:10022/aea/material-models/cpp_stub.git

3. Navigate to your preferred repository directory on your local computer. In a
   terminal, you can follow the example ``sstelmo`` session below

   .. code:: bash
   
      # Start an ssh session to sstelmo.lanl.gov
      $ ssh -X sstelmo.lanl.gov
   
      # Note the present working directory (pwd) is your home directory
      $ pwd
      /home/<moniker>
   
      # OPTIONAL: Create a project space repository directory
      $ mkdir -p /projects/$USER/w13repos
      # Change pwd to repository directory
      $ cd /projects/$USER/w13repos
      # Double check pwd is repository parent directory
      $ pwd
      /projects/<moniker>/w13repos

4. Clone the stub repository using the URL copied in step 2.

   .. code:: bash
   
      # Double check pwd is repository parent directory
      $ pwd
      /projects/<moniker>/w13repos
   
      # Clone the stub repository
      $ git clone ssh://git@re-git.lanl.gov:10022/aea/material-models/cpp_stub.git

5. Rename the local repository directory for your project.

   .. code:: bash
   
      # Double check pwd is repository directory
      $ pwd
      /projects/<moniker>/w13repos
   
      # Observe the stub repo directory name
      $ ls cpp_stub -d
      cpp_stub
   
      # Rename the stub repo directory after your project
      $ mv cpp_stub my_project
   
      # Observe that the stub repo directory no longer exists
      $ ls cpp_stub -d
      ls: cannot access 'cpp_stub': No such file or directory
   
      # Observe that your project directory exists
      $ ls my_project -d
      my_project

6. Change to your project's local repository directory

   .. code:: bash
   
      # Double check pwd is repository directory
      $ pwd
      /projects/<moniker>/w13repos
   
      # Change to your project directory
      $ cd my_project
   
      # Double check pwd is your project directory
      $ pwd
      /projects/<moniker>/w13repos/my_project

******************************************
Create a new upstream repository on Gitlab
******************************************

1. Navigate to the W-13 `Material Models`_ Gitlab sub-group.

2. Create a new repository by clicking on the blue "New project" button in the
   upper right corner of the sub-group main page.

   .. note::

      If you do not have the "Developer" or "Maintainer" role assigned to you in
      this sub-group, you will not be able to create a new project directly. You can
      request a role change from the `Material Models`_ sub-group owners. Sub-group
      owners may prefer to create a project for you and make you the owner of that
      project. You can check the `Material Models members`_ list for contact
      information.

3. On the "Create new project" page, follow the link for "Create blank project".

   .. note::

      Gitlab offers a feature to create template projects that may make this
      guide much simpler in the future. Contact the ``cpp_stub`` developers and `AEA
      Gitlab group`_ owners to discuss progress on simplified repository setup using
      templates.

3. Enter a name for your project in the "Project name" field. Optionally add a
   "project description" and click the blue "Create project" button.

4. Follow the "Push an existing Git repository" instructions at the bottom of
   the new project webpage.

   .. code:: bash
   
      $ pwd
      /projects/<moniker>/w13repos/my_project
      $ git remote rename origin old-origin
      $ git remote add origin ssh://git@re-git.lanl.gov:10022/aea/material-models/dummy.git
      $ git push -u origin --all
      $ git push -u origin --tags

5. Refresh the Gitlab project webpage and verify that the repository code was pushed
   correctly. You should see a list of source files and this Bitbucket parsed
   ``README.rst`` displayed. You can also review the "master" and "dev" branch from
   the left hand side bar "Repository" > "Branches" menu and the Git tags from the
   "Repository" > "Tags" menu.

6. Remove any issue branches from the ``cpp_stub`` project on the "Repository" >
   "Branches" menu. You should keep only the "master" and "dev" branches.

7. If everything looks correct on Gitlab project, you can clean up your local
   repository.

   .. warning::

      WARNING: the ``-D`` option FORCE deletes branches. Triple check the
      command and use with caution. If you're uncertain about this step, contact the
      cpp_stub developers for help.

   .. code:: bash

      # Remove the cpp_stub remote
      $ git remote remove old-origin

      # Ensure that you're on the master branch
      $ git checkout master

      # Remove all the cpp_stub issue branches
      $ git branch | grep -v "master\|dev" | xargs git branch -D 

***********************************
Update settings for your repository
***********************************

Gitlab repositories (a.k.a. 'projects') in the `Material Models`_ Gitlab
sub-group inherit permissions and settings from that sub-group.  This included
read permission for the ``w13bitbucket`` UNIX group (`W-13 Managed UNIX
Groups`_). For most developers, these inherited repository settings are
appropriate and only a small number of settings must be updated.

1. Click on the gear icon labelled "Settings" in the lower left sidebar.

2. Click on the "Repository" menu item that appears in the left sidebar

3. From the "Default branch" > "Expand" page, update the default branch from
   "master" to "dev".

4. From the "Protected branches" > "Expand" page, protect the "master" and "dev"
   branches according to the needs of your project.

5. From the "Project Information" > "Members" item at the top of the left side
   bar you can add additional permissions by user and UNIX group.

**********************************************
Update the remote URL in your local repository
**********************************************

The final repo setup step is to update the remote URL of the local clone of
``my_project``.  We will return to the terminal session.

1. Copy the URL of your "remote" repository from the Bitbucket webpage. It
should look like:

.. code:: bash

   ssh://git@re-git.lanl.gov:10022/aea/material-models/my_project.git

2. Return to your terminal session and update the remote repository for the
   final time.

.. code:: bash

   $ pwd
   /projects/<moniker>/w13repos/my_project
   $ git remote set-url origin ssh://git@re-git.lanl.gov:10022/aea/material-models/my_project.git
   $ git push -u origin --all
   $ git push origin --tags

*****************************************
Update project name throughout repository
*****************************************

.. note::

   Note: the remaining steps are a truncated version of the W-13 Git project
   `ECMF contribution guide`_.  which can also be found in the `W-13 DevOps
   Manual`_.  Critically, these steps will omit the Jira task creation and
   Bitbucket Pull-Request (PR) steps. The Bitbucket PR steps may be reproduced
   using the contribution guide, but your project will have to create a Jira
   project prior to integrating the Jira workflow. Contact the `AEA Gitlab group`_
   owners to create new sub-groups or projects.  You can email the W-13 DevOps team
   w13devops@lanl.gov for notes about setup.

1. Create a feature branch for your project name updates

.. code:: bash

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

.. code:: bash

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

.. code:: bash

   $ pwd
   /projects/<moniker>/w13repos/my_project

   # Replace lower case occurrences in place
   $ sed -i 's/cpp_stub/my_project/g' $(grep -ri cpp_stub . --exclude-dir={build,.git} -l)
   $ grep -ri cpp_stub . --exclude-dir={build,.git} -l
   ./src/cpp/cpp_stub.h

   # Replace upper case occurrences in place
   $ sed -i 's/CPP_STUB/MY_PROJECT/g' $(grep -ri cpp_stub . --exclude-dir={build,.git} -l)

4. Verify no more occurrences of project name ``cpp_stub``

.. code:: bash

   $ pwd
   /projects/<moniker>/w13repos/my_project
   $ grep -ri cpp_stub . --exclude-dir={build,.git} | wc -l
   0
   $ grep -ri cpp_stub . --exclude-dir={build,.git}
   # no stdout to terminal because no files found
   $ grep -ri cpp_stub . --exclude-dir={build,.git} -l
   # no stdout to terminal because no files found

5. Search and replace camelcase project name occurrences, e.g. ``cppStub``.

.. code:: bash

   $ grep -r cppStub . --exclude-dir={build,.git}
   ...
   $ sed -i 's/cppStub/myProject/g' $(grep -r cppStub . --exclude-dir={build,.git} -l)
   $ grep -r cppStub . --exclude-dir={build,.git} -l
   # no stdout to terminal because no files found

6. Find files containing the project in their file name

.. code:: bash

   $ pwd
   /projects/<moniker>/w13repos/my_project
   $ find . -type d \( -name .git -o -name build \) -prune -false -o -name "*cpp_stub*"
   ./src/cpp/cpp_stub.cpp
   ./src/cpp/cpp_stub.h
   ./src/cpp/tests/test_cpp_stub.cpp

7. Rename files after current project

.. code:: bash

   $ rename 's/cpp_stub/myproject/' $(find . -type d \( -name .git -o -name build \) -prune -false -o -name "*cpp_stub*")

8. Commit and push your changes to your "remote" or "fork" repository

.. code:: bash

   $ pwd
   /projects/<moniker>/w13repos/my_project
   # Add tracked files and message
   $ git commit -a -m "FEAT: replace cpp_stub with my_project through repository"
   $ git push origin feature/project-name-updates

You can also perform some cleanup in ``README.md`` to remove this walk-through.

From here, the W-13 best practice workflow would return to the Bitbucket webpage
and submit a Pull-Request from the ``feature/project-name-updates`` branch of
"\<moniker\>/my_project" repository (a.k.a. fork or remote) to the ``dev`` branch
of "Material Models/my_project" repository (a.k.a. upstream or official).

After updating your project by merging to the upstream repository, the fork
syncing feature of Bitbucket will automatically update any identically named
branches in your fork repository. Best practices suggest you should limit the
upstream repository branches to clean ``dev`` and ``master`` branches and
*NEVER* develop directly on the ``dev`` and ``master`` branches of your fork
repository. Limit development work to ``feature/thing`` type branches on your
fork/remote repo and frequently commit changes and push from the local feature
branch back to the fork/remote repo.

Happy hacking!
