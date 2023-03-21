.. _changelog:


#########
Changelog
#########

******************
0.1.0 (unreleased)
******************

Breaking changes
================
- Remove linear elasticity, which now lives in stress tools (:issue:`1`, :merge:`9`). By `Kyle Brindley`_.

Internal Changes
================
- Change the build system to more closely mimic the current state of cpp_stub (:merge:`15`). By `Nathan Miller`_.

New Features
============
- Added linear elasticity for the internal particle (:merge:`3`). By `Nathan Miller`_.
- Added the calculation of the current distance between two points in the particles (:merge:`4`). By `Nathan Miller`_.
- Added the calculation of the decomposition of a vector into normal and tangential components w.r.t. another vector (:merge:`4`). By `Nathan Miller`_.
- Added the calculation of the linear traction separation energy (:merge:`4`). By `Nathan Miller`_.
- Added the calculation of the gradients of the traction separation energy w.r.t. the parameters (:merge:`5`). By `Nathan Miller`_.
- Added the calculation of the traction separation traction - Cauchy stress traction constraint equation (:merge:`6`). By `Nathan Miller`_.
- Added the calculation of Nanson's relation for mapping reference areas to the current configuration (:merge:`7`). By `Nathan Miller`_.
- Added the calculation of the Lagrangian for the overlap particle (:merge:`8`). By `Nathan Miller`_.
- Extended the Lagrangian for the overlap particle to have a radius other than 1 (:merge:`8`). By `Nathan Miller`_.
- Extended the number of gradients computed in the Lagrangian for the overlap particle (:merge:`8`). By `Nathan Miller`_.
- Added the calculation of the amount of overlap of a non-local and local particles (:merge:`8`). By `Nathan Miller`_.
- Added the calculation of the gradients for the local and non-local particles (:merge:`8`). By `Nathan Miller`_.
- Added the calculation of the fourth-order gradients for the Lagrangian for the overlap particle. This is required for the constraint equation. (:merge:`8`). By `Nathan Miller`_.
- Added the calculation of the third-order gradients for the solution of the overlap distance. This is required for the constraint equation. (:merge:`8`). By `Nathan Miller`_.
- Added the decomposition of a sphere for the purposes of integration and contact detection. (:merge:`10`). By `Nathan Miller`_.
- Added the capability to integrate 2D quadratic elements and surface meshes composed of quadratic elements. (:merge:`11`). By `Nathan Miller`_.
- Added a more general form of the distance calculation to support random deformations. (:merge:`12`). By `Nathan Miller`_.
- Initial commit of the asp base class with a linear elastic local particle energy definition. (:merge:`14`). By `Nathan Miller`_.
- Added the computation of the traction separation energy at a single surface point. (:merge:`15`). By `Nathan Miller`_.
- Changed the interface to the surface energy such that it sets an internal variable rather than returning a value (:merge:`16`). By `Nathan Miller`_.
