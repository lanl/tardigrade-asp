.. _releasephilosophy:

Release Philosophy
==================

This section discusses topics related to |project| releases and version numbering.

Release and Deployment
++++++++++++++++++++++

.. warning::

   The |project| does not yet have a template ``cmake3 --install`` definition and is not yet deployable to the W-13
   Python Environments. Instead, users are directed to the :ref:`user_manual` and :ref:`build` section of the
   :ref:`devops_manual`. 

The |project| is built and installed as a c++ library in the `W-13 Python Environments`_ available hamming, sstelmo, and
any local linux machines with home and project drives mapped from the W-13 NFS server. These are Anaconda Python 3
environments with installed packages required for W-13 software development and engineering analysis.  There are two
versions of the W-13 Python Environments:

1) release
2) beta

The release environment contains the deployed final release versions of W-13 software projects, e.g. Toolbox and ECMF,
as tested against the installed Python packages. The beta environment contains the deployed developer release of W-13
software projects. While the deployed projects in beta have been unit and integration tested, the beta environment may
include updated or new Python modules that result in less stable behavior. The beta environment is used to test W-13
software projects against an updated Python environment before releasing the new environment.

Version Numbers
+++++++++++++++

The ECMF project follows the `PEP-440`_ standard for version numbering. The
final release version number uses the three component ("major.minor.micro")
scheme. The developer (a.k.a. dev or beta) version number follows the final
release number with an appended "+dev" local version number. The version numbers
correspond to git tags in the `Official ECMF Repo`_ which point to a static
release of the ECMF project.

Because the deployed release of the developer version is constantly updated
against development work, the version number found in the deployed release
contains additional information. During deployment, the developer version number
is appended with the git information from the most recent build. This
information contains the most recent git tag ("major.minor.micro+dev") followed
by the number of commits since the last final release and a short hash.

Major Number
------------

The major number is expected to increment infrequently. Incrementing the major
number requires a decision from the ECMF developement council.

Minor Number
------------

The minor number is updated for the following reasons:

* Third party software support
* New workflow modules
* Features spanning multiple module dependencies
* UX modifications

Incrementing the minor version requires a manual update to the release number on
a dedicated release commit. Minor version number increments should be made after
a decision from the ECMF development council. Review for content of the next
minor release takes place at the ECMF development council every other month and
requires approval by relevant stakeholders.

Micro Number
------------

The micro number is automatically incremented after any merge from the
development (dev) branch into the release (master) branch. The micro version
number indicates the following changes:

* Bug fixes
* Single module modifications
* UI modifications
* Quorum of completed Jira tasks with ECMF development concil approval by 2/3rd
  majority

Review for content of upcoming micro releases occurs at every ECMF development
council meeting: 24 per annum.

.. _releasebranchreq:

Release Branch Requirements
---------------------------

Major and minor version number updates require a release branch.
Releases correspond to a variety of bug fixes and features that characterize
the release, as documented in :ref:`whatsnew`.

Steps needed for a release include:

1. Update version number in ``ecmf/settings.py``, e.g., ``VERSION = 'X.Y.Z'``.
   Version bumps should be accompanied by resetting numbers to the right of the
   bump to zero, e.g., ``'0.2.30'`` to ``'0.3.0'`` and ``'1.2.30'`` to
   ``'2.0.0'``.

2. Modify ``docs/whats_new.rst`` to move version number for release PR commit and
   add description as relevant.
3. Commit changes and submit a pull request at `ECMF Bitbucket`_.
