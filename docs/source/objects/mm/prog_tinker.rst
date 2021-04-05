
Tinker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tinker :cite:`Rackers2018` is a complete and general package for molecular mechanics and dynamics, with some special
features for biopolymers. Tinker has the ability to use any of several common parameter sets, such
as Amber (ff94, ff96, ff98, ff99, ff99SB), CHARMM (19, 22, 22/CMAP), Allinger MM (MM2-1991 and
MM3-2000), OPLS (OPLS-UA, OPLS-AA), Merck Molecular Force (MMFF), Liam Dang's polarizable model,
AMOEBA (2004, 2009, 2013, 2017, 2018) polarizable atomic multipole force fields, AMOEBA+ that adds
charge penetration effects, and HIPPO (Hydrogen-like Interatomic Polarizable POtential) force field.

+------------------------+------------------------------------------------+---------------------+
| Parameters             | Work                                           | Default             |
+========================+================================================+=====================+
| **molecule**           | Molecule object                                |                     |  
| (:class:`Molecule`)    |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **scheme**             | Type of QM/MM scheme                           | *None*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **embedding**          | Charge embedding options                       | *None*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **vdw**                | Van der Waals interactions                     | *None*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **l_periodic**         | Use periodicity in the calculations            | *False*             |
| *(boolean)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **cell_par**           | Cell lattice parameters (lengths and angles)   | *6 \* [ 0.0 ]*      |
| *(double, list)*       |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **xyz_file**           | Initial tinker.xyz file                        | *'./tinker.xyz'*    |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **key_file**           | Initial tinker.key file                        | *'./tinker.key'*    |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **mm_path**            | Path for MM binary                             | *'./'*              |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **nthreads**           | Number of threads in the calculations          | *1*                 |
| *(integer)*            |                                                |                     |
+------------------------+------------------------------------------------+---------------------+
| **version**            | Version of Tinker                              | *'8.7'*             |
| *(string)*             |                                                |                     |
+------------------------+------------------------------------------------+---------------------+

Detailed description of parameters
''''''''''''''''''''''''''''''''''''

- **scheme** *(string)* - Default: *None*

  This parameter specifies type of the QM/MM scheme. The current version of PyUNIxMD supports two types of the QM/MM scheme.
  One is the subtractive scheme and the other is the additive scheme. The detailed expressions for
  these schemes are given in the literature. :cite:`Senn2009`

  + *'additive'*: Total three calculations was carried out, which correspond
    to a QM calculation on the inner region, an MM calculation on the outer region,
    and an explicit QM-MM coupling term, respectively. In general, it is
    recommended to use the additive scheme for solution systems.
  + *'subtractive'*: Total three calculations was carried out, which correspond
    to a QM calculation on the inner region, an MM calculation on the entire region,
    and an MM calculation on the inner region, respectively. In this scheme,
    no explicit QM-MM coupling terms are needed since they are implicitly included
    in three calculations. It is recommended to use the subtractive scheme for proteins.

\

- **embedding** *(string)* - Default: *None*

  This parameter specifies type of charge-charge interactions between the inner and outer regions.
  Current version of PyUNIxMD supports two types of charge-charge embedding.
  One is mechanical interaction and the other is electrostatic interaction.
  The **embedding** of the MM object must be same with the **embedding** defined in the QM object.
  If this parameter is *None*, the charge-charge embedding is not included in the QM/MM calculation.

  + *'mechanical'*: The charge-charge interactions are treated at MM level.
    The energies are calculated from the interactions between point charges.
  + *'electrostatic'*: The charge-charge interactions are treated at QM level.
    The point charges of the outer regions are added to the one-electron terms of the
    Hamiltonian in QM calculation. Hence, the polarization effect from the point charges are considered in this embedding.

\

- **vdw** *(string)* - Default: *None*

  This parameter specifies type of van der Waals interactions between the inner and outer regions.
  Current version of PyUNIxMD supports one type of van der Waals interaction,
  which is the Lennard-Jones interaction. The other types of van der Waals
  interactions provided in Tinker are not currently interfaced with PyUNIxMD.
  If this parameter is *None*, the van der Waals interactions are not included in the QM/MM calculation.

  + *'lennardjones'*: The Lennard-Jones interactions are used for van der Waals interactions.

\

- **l_periodic** *(boolean)* - Default: *False*

  When **l_periodic** is set to *True*, a periodicity is considered in the calculation.
  Only :math:`\Gamma`-point sampling is supported with Tinker when the periodicity is considered.

\

- **cell_par** *(double, list)* - Default: *6 \* [ 0.0 ]*

  This parameter specifies cell lattice parameters of the periodic unit cell.
  The list consists of six elements and first three elements correspond to
  the :math:`a`, :math:`b`, and :math:`c` lengths while the last three
  elements correspond to the :math:`\alpha`, :math:`\beta`, and :math:`\gamma` angles, respectively.

\

- **xyz_file** *(string)* - Default: *'./tinker.xyz'*

  This parameter specifies a TINKER Cartesian coordinates file of the initial geometry.
  This file must include correct atom types and bonding information.

\

- **key_file** *(string)* - Default: *'./tinker.key'*

  This parameter specifies initial 'tinker.key' file used in the calculations. The keywords of Tinker
  except **embedding**, **vdw**, and the periodicity can be included in this file.
  For example, if you want to add some constraints to the systems, then
  the related keywords can be added to the 'tinker.key' file.

\

- **mm_path** *(string)* - Default: *'./'*

  Path for Tinker binary. In our interfacing scripts, the executable file,
  'testgrad' (or 'testgrad.x') is used to calculate the energies and forces in MM level.

  This parameter determines a path for Tinker binaries such as testgrad.
  (For example, '/my_disk/my_name/Tinker/bin/').

\

- **nthreads** *(integer)* - Default: *1*

  This parameter specifies number of threads in the calculation. To use this option, you must check
  that your binaries of Tinker supports OpenMP parallelization.

\

- **version** *(string)* - Default: *'8.7'*

  This parameter determines the version of Tinker.
  PyUNIxMD is currently based on version 8.7 of Tinker.

