=============
Theory Manual
=============

We utilize a multiplicative decomposition of the deformation gradient of the form

.. math::

    \mathbf{F} = \mathbf{F}^e \mathbf{F}^d \mathbf{F}^i \mathbf{F}^{\theta}

where :math:`\mathbf{F}^e` is the elastic part, :math:`\mathbf{F}^d` is the damage part, :math:`\mathbf{F}^i` is the irrecoverable
strain ( plasticity) part, and :math:`\mathbf{F}^{\theta}` is the thermal part ( both irrecoverable and recoverable ) of the
deformation gradient.  Following Simo and Hughes :cite:`bib:simo2000`, the Helmholtz free energy is broken into
volumetric :math:`U^e` and isochoric :math:`\bar{\Psi}^e` parts where the viscoelastic response is assumed to be a
isochoric response alone. We define the Helmholtz free energy as

.. math::

    \Psi^e = U^e\left( J^e \right) + \sum_{i=1}^{N^{ve}} \bar{\Psi}_i \left( \mathbf{\bar{E}}^e,
    \mathbf{\bar{\Xi}}_i \right) + \bar{\Psi}_{\infty}^e\left( \bar{E}^e \right)

where :math:`J^e` is the Jacobian of the elastic part of the deformation gradient, :math:`\mathbf{\bar{E}}^e` is the
Isochoric elastic Green-Lagrange strain, :math:`\mathbf{\Xi}_i` is the :math:`i` th viscous backstrain, and
:math:`\bar{\Psi}_{\infty}^e` is the Helmholtz free energy associated with the fully relaxed configuration.

In order to enforce positivity in the Clausius-Duhem inequality we require that

.. math::

    -\frac{\partial ( \rho^0 \Psi ) }{\partial \mathbf{\xi} } \dot{\mathbf{\xi}} \geq 0

where :math:`\mathbf{\xi}` is the vector of strain-like internal state variables. We can define the potential function of
the irreversible processes as

.. math::

    \begin{align*}
    \frac{1}{2} \mathbf{\xi} \mathbf{H} \mathbf{\xi} &= -( \rho_0 \Psi )^{irr}
    \rightarrow \mathbf{q} = \mathbf{H} \mathbf{\xi}
    \end{align*}

where :math:`\mathbf{q}` is the vector of stress-like internal state variables and :math:`\mathbf{H}` is a collection of moduli
which can allow mixing between irrecoverable processes if desired.

The elastic damage evolves via

.. math::

    \begin{align}
    \dot{D} &= \dot{D}_0 \left( \theta \right) \left \langle \frac{f^d}{q^e} \right \rangle^{n^d}\\
    f^d &= \sigma^{vm} - A^d \bar{\sigma}\\
    \dot{\xi}^d &= \dot{D} h^d
    \end{align}

where :math:`\dot{D}_0\left( \theta \right)` is the temperature effect of the evolution of elastic damage, :math:`n^d`
controls the rate effect of the elastic damage evolution, :math:`\sigma^{vm}` is the Von-Mises stress,
:math:`\bar{\sigma}` is the mean stress, :math:`A^d` is the damage pressure term, and :math:`h^d` is the elastic damage
hardening curve. The elastic damage is applied via

.. math::

    \mathbf{E}^d = \frac{D}{1 - D}\mathbf{E}^e

The irrecoverable strain evolves via the irrecoverable velocity gradient which is defined as

.. math::

    \mathbf{L}^i = \gamma \mathbf{n}^i

where the plastic multiplier, :math:`\gamma` and flow direction :math:`\mathbf{n}^i` are defined as

.. math::

    \begin{align}
    \gamma &= \gamma_0 \left( \theta \right) \left\langle \frac{f^i}{q^i} \right\rangle^{n^i}\\
    f^i &= \sigma^{vm} - A^i \bar{\sigma}\\
    \mathbf{n}^i &= \frac{\partial}{\partial \mathbf{\sigma} } \left( \sigma^{vm} - B \bar{\sigma} \right)\\
    \dot{\xi}^i &= \gamma h^i
    \end{align}

where :math:`\gamma( \theta )` is the temperature effect on the evolution of :math:`\gamma`, :math:`n^i` is the term
which controls the rate effect of the irrecoverable strain evolution, :math:`A^i` is the irrecoverable deformation
pressure term, and :math:`B` is the flow direction parameter.

The irrecoverable thermal expansion is assumed to evolve via a driving force, :math:`q^{\theta d}`, which is assumed to
be proportional to the volumetric part of the reversible thermal strain :math:`\mathbf{E}^{\theta r}`. We define the
evolution of the irrecoverable thermal strain through an isotropic, with respect to temperature, hardening state
variable :math:`\xi^{\theta i}` and a kinematic hardening state variable :math:`\xi^{\theta k}` as

.. math::

    \begin{align}
    \mathbf{L}^{\theta i} &= \gamma^{\theta} \mathbf{n}\\
    \mathbf{n} &= \frac{\mathbf{e}^{\theta r} }{\left|\left| \mathbf{e}^{\theta r} \right|\right|}\\
    \gamma^{\theta} &= \frac{1}{\eta( \theta )} \left \langle f^{\theta} \right \rangle\\
    \dot{\xi}^{\theta i} &= \gamma^{\theta} h^{\theta i}\\
    \dot{\xi}^{\theta k} &= \gamma^{\theta} h^{\theta k} \text{sign}\left( q^{\theta d} - q^{\theta k} \right)\\
    f^{\theta} &= \left( \left| q^{\theta d} - q^{\theta k} \right| - q^{\theta i} \right ) f\left( \bar{\sigma} \right)
    \end{align}

where :math:`\eta` is a viscous damping term and :math:`f\left( \bar{\sigma} \right)` is a pressure suppression term.
