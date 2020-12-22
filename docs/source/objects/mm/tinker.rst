
Tinker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tinker :cite:`Rackers2018` is a complete and general package for molecular mechanics and dynamics, with some special
features for biopolymers. Tinker has the ability to use any of several common parameter sets, such
as Amber (ff94, ff96, ff98, ff99, ff99SB), CHARMM (19, 22, 22/CMAP), Allinger MM (MM2-1991 and
MM3-2000), OPLS (OPLS-UA, OPLS-AA), Merck Molecular Force (MMFF), Liam Dang's polarizable model,
AMOEBA (2004, 2009, 2013, 2017, 2018) polarizable atomic multipole force fields, AMOEBA+ that adds
charge penetration effects, and HIPPO (Hydrogen-like Interatomic Polarizable POtential) force field.

.. note:: Our interface script is generated with 8.7 version of Tinker program. Tinker can be
   used for QM/MM calculation. UNI-xMD supports two types of QM/MM scheme. One is ``additive``
   scheme and the other is ``subtractive`` scheme, which are controlled by ``scheme`` variable.

.. note:: Currently, ``embedding`` and ``vdw`` variables control the charge-charge interactions
   and the van der Walls interactions, respectively. For ``embedding`` option, ``electrostatic``
   and ``mechanical`` are possible choices, while ``lennardjones`` becomes the only option for
   ``vdw`` option.

.. note:: For ``cell_par`` variable, it reads a list variable consisted of 6 elements,
   which correspond to lengths and angles.

+-------------------+------------------------------------------------+---------------------+
| Keywords          | Work                                           | Default             |
+===================+================================================+=====================+
| ``scheme``        | type of QM/MM scheme                           | ``None``            |
+-------------------+------------------------------------------------+---------------------+
| ``embedding``     | charge embedding options                       | ``None``            |
+-------------------+------------------------------------------------+---------------------+
| ``vdw``           | van der Walls interactions                     | ``None``            |
+-------------------+------------------------------------------------+---------------------+
| ``periodic``      | use periodicity in the calculations            | ``False``           |
+-------------------+------------------------------------------------+---------------------+
| ``cell_par``      | cell lattice parameters (lengths and angles)   | ``6 * [ 0.0 ]``     |
+-------------------+------------------------------------------------+---------------------+
| ``xyz_file``      | initial tinker.xyz file                        | ``./tinker.xyz``    |
+-------------------+------------------------------------------------+---------------------+
| ``key_file``      | initial tinker.key file                        | ``./tinker.key``    |
+-------------------+------------------------------------------------+---------------------+
| ``mm_path``       | path for MM binary                             | ``./``              |
+-------------------+------------------------------------------------+---------------------+
| ``nthreads``      | number of threads in the calculations          | ``1``               |
+-------------------+------------------------------------------------+---------------------+
| ``version``       | version of Tinker program                      | ``8.7``             |
+-------------------+------------------------------------------------+---------------------+

