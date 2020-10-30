
If you want to interface other quantum chemistry programs, you can write your own interfacing scripts.

In [UNIXMDHOME]/src/qm/diy/ directory, there are four template scripts containing essential components for interfacing.

First, you should define a class for the quantum chemistry program that inherits the QM_calculator class. (``program_name.py``)

.. literalinclude:: ../../../src/qm/diy/program_name.py

Then you define classes for theories (ex. DFT, HF, CASSCF, ...) that you want to implement from the program. (``theory1.py``)

.. literalinclude:: ../../../src/qm/diy/theory1.py

Class for theory (Theory1) inherits program class (Program_name) and its variables ``molecule.l_nacme`` and ``self.re_calc`` must be initialized.
If the program can calculate NACVs and you want to use them, then ``molecule.l_nacme`` should be set to ``False``.
In contrast, if the program cannot calculate NACVs, ``molecule.l_nacme`` should be set to ``True`` and time-derivative NAC matrix elements should be provided instead of NACVs.
If the program can calculate gradients of multiple electronic states at a single run, you can set ``self.re_calc`` to ``False``.

Theory class must have ``get_data`` method because ``run`` method in dynamics object calls ``get_data`` method.
The method must get all arguments shown in the example code.

Arguments in ``get_data`` method is explained in modules section.

