
Jaynes_Cummings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In cavity quantum electrodynamics (cQED), strong coupling between a photon mode and the electronic states of a
molecule inside a cavity occurs, provided that the coherent energy exchange between them
is faster than the decay of either component. In this regime, the parent molecular electronic
states hybridize with the cavity photon, and new quantized states, polaritons, emerge.
At the theoretical level, the simplest case of hybrid light-matter states occurring from
coupling between a two-level system and a photon is described by the Jaynes–Cummings (JC) model. :cite:`Jaynes1963`
The JC model considers the interaction between a two-level system with a lossless cavity mode with
frequency :math:`\omega_c` under rotating-wave approximation as

.. math::

   \hat{H}_\text{JC} = \hat{H}_\text{mol} + \hbar \omega_c \left(\hat{a}^\dagger_c \hat{a}_c + 1/2\right)
   + g \vec{\lambda} \cdot \vec{\mu} \left(\hat{a}^\dagger_c + \hat{a}_c\right)

.. note:: For accurate propagation of the nuclei, the gradients of transition dipole moments are required.
   Detailed description of forces and nonadiabatic coupling vectors with polaritonic states is in :cite:`Lee2024`

+------------------------+--------------------------------------------------+---------------------+
| Parameters             | Work                                             | Default             |
+========================+==================================================+=====================+
| **polariton**          | Polariton object                                 |                     |  
| (:class:`Polariton`)   |                                                  |                     |
+------------------------+--------------------------------------------------+---------------------+
| **coupling_strength**  | Coupling strength for cavity-matter interaction  | *0.001*             |
| *(double)*             |                                                  |                     |
+------------------------+--------------------------------------------------+---------------------+
| **force_level**        | Level for calculation of force of target         | *'hf'*              |
| *(string)*             | polaritonic state                                |                     |
+------------------------+--------------------------------------------------+---------------------+
| **l_check_crossing**   | Logical to check diabatic character of           | *False*             |
| *(boolean)*            | polaritonic states                               |                     |
+------------------------+--------------------------------------------------+---------------------+
| **l_crt**              | Logical to include the counter-rotating terms    | *False*             |
| *(boolean)*            |                                                  |                     |
+------------------------+--------------------------------------------------+---------------------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **coupling_strength** *(double)* - Default: *0.001*

  This parameter specifies the strength for the light-matter interaction.
  Users carefully controls this parameter to observe the strong coupling regime in cQED.

\

- **force_level** *(string)* - Default: *'hf'*

  This parameter specifies an option for calculation of the force of the target polaritonic state.
  Current version of PyUNIxMD supports three options which can require more computational efforts.
  The detailed expressions for these options are given in the literature. :cite:`Lee2024`

  + *'diag_only'*: The total force is approximated by diagonal term only in the derivatives of
    the Hamiltonian matrix elements from the Eq. (5) in the reference :cite:`Lee2024`.
    Because the accurate evaluation of the force needs the transition dipole gradients,
    the force is simply obtained by neglecting the off-diagonal elements in the expression.
  + *'hf'*: The total force is approximated by the Hellmann-Feynman contribution alone.
    In this case, only nonadiabatic coupling vectors between the BO states are needed
    instead of the calculation of transition dipole gradients. The detailed equation is given
    in the Eq. (6) of the reference :cite:`Lee2024`.
  + *'full'*: The total force is obtained by evaluating all contributions
    from the Eq. (5) in the reference :cite:`Lee2024`. In order to use this option,
    the gradients of the transition dipole moments must be provided from the QM methods.

\

- **l_check_crossing** *(boolean)* - Default: *False*

  When **l_check_crossing** is set to *True*, the diabatic character of adjacent polaritonic
  states will be checked. In general, the trivial crossings can occur due to the light-matter interaction
  even for the molecular system. To describe correct nuclear dynamics, the trivial crossing should
  be treated in polariton dynamics.

\

- **l_crt** *(boolean)* - Default: *False*

  When **l_crt** is set to *True*, the counter-rotating terms are included in the Hamiltonian.
  The JC model uses the rotating-wave approximation, which treats low-frequency interaction
  terms in the light-matter interaction. Due to this approximation, only lower polaritonic and
  upper polaritonic states are mixed. However, if the counter-rotating terms (high-frequency)
  are included (i.e., quantum Rabi model), the ground and highest polaritonic states can be mixed.

