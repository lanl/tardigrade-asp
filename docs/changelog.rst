.. _changelog:

Changelog
=========

0.0.2 (unreleased)
------------------

Breaking changes
~~~~~~~~~~~~~~~~
- Remove testing and support for intel ``icpc`` compiler (:jira:`516`, :pull:`9`). By `Kyle Brindley`_.

New Features
~~~~~~~~~~~~
- Add do-nothing template c++ Abaqus UMAT interface and sample Abaqus input file (:jira:`502`, :pull:`6`). By `Kyle Brindley`_.
- Use example c++ library in Abaqus UMAT template (:jira:`505`, :pull:`8`). By `Kyle Brindley`_.

Bug fixes
~~~~~~~~~

Documentation
~~~~~~~~~~~~~
- Add changelog to documentation (:jira:`450`, :pull:`11`). By `Kyle Brindley`_.
- Add direct CMake build instructions and minimal user manual (:jira:`519`, :pull:`12`. By `Kyle Brindley`_.

Internal Changes
~~~~~~~~~~~~~~~~
- Use BOOST and ctest for unit testing (:jira:`357`, :pull:`4`). By `Kyle Brindley`_.
- Update Jenkins CI configuration and store with version controlled repository (:jira:`442`, :pull:`5`). By `Kyle Brindley`_.
- Demonstrate c++ ``vector_tools`` library for unit testing (:jira:`506`, :pull:`7`). By `Kyle Brindley`_.
- Add integration tests for Abaqus UMAT interface (:jira:`504`, :pull:`10`). By `Kyle Brindley`_.

Enhancements
~~~~~~~~~~~~

0.0.1 (2020-10-26)
------------------

Breaking changes
~~~~~~~~~~~~~~~~

New Features
~~~~~~~~~~~~
- Create c++ stub repository targeting constitutive modeling (:jira:`332`, :pull:`1`). By `Kyle Brindley`_.

Bug fixes
~~~~~~~~~

Documentation
~~~~~~~~~~~~~

Internal Changes
~~~~~~~~~~~~~~~~
- Add continuous integration scripts (:jira:`333`, :pull:`2`). By `Kyle Brindley`_.

Enhancements
~~~~~~~~~~~~
