Input File Description
======================

Here we provide a description of the input variables related to each of the applications in the code.
Phoebe tries to check the variables provided in input, and stops the code if it detects missing or invalid input parameters.

Note that the executable, ``phoebe`` is the same for all applications. The applications are selected through the :ref:`appName` parameter.

.. note::
   The input file text is case sensitive!

-------------------

QE to Phoebe
------------

**Functionality:** Convert the electron-phonon coupling from Quantum-ESPRESSO format to Phoebe format.
This postprocesses the data for  Wannier interpolation or EPA calculations in Phoebe.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "elPhQeToPhoebe"

* :ref:`phD2FileName`

* :ref:`elPhInterpolation`

* :ref:`electronH0Name`

* :ref:`wannier90Prefix`

* :ref:`quantumEspressoPrefix`

* :ref: `epaMinEnergy`

* :ref: `epaMaxEnergy`

* :ref: `epaNumBins`

* :ref: `epaSmearingEnergy`

.. raw:: html

  <h3>Sample input file (Wannier interpolation)</h3>

::

  appName = "elPhQeToPhoebe"
  elPhInterpolation = "wannier"

  phD2FileName = "./silicon.fc"
  electronH0Name = "./si_tb.dat",
  wannier90Prefix = "si"
  quantumEspressoPrefix = "silicon"

.. raw:: html

  <h3>Sample input file (EPA)</h3>

::

  appName = "elPhQeToPhoebe"
  elPhInterpolation = "epa"

  phD2FileName = "./silicon.fc"
  electronH0Name = "./out/silicon.xml",
  quantumEspressoPrefix = "silicon"

  electronFourierCutoff = 4.
  epaMinEnergy = -4. eV
  epaMaxEnergy = 10. eV
  epaNumBins = 10
  epaSmearingEnergy = 0.05 eV


-----------------------------------

Phonon BTE Solver
----------------------

**Functionality:** Build and solve the phonon Boltzmann Transport Equation (BTE) to compute phonon transport properties such as thermal conductivity, relaxation times and viscosity.

.. raw:: html

  <h3>Input variables</h3>

:ref:`appName` = "phononTransport"

* :ref:`phD2FileName`

* :ref:`phD3FileName`

* :ref:`sumRuleD2`

* :ref:`qMesh`

* :ref:`temperatures`

* :ref:`minTemperature`

* :ref:`maxTemperature`

* :ref:`deltaTemperature`

* :ref:`smearingMethod`

* :ref:`smearingWidth`

* :ref:`solverBTE`

* :ref:`scatteringMatrixInMemory`

* :ref:`windowType`

* :ref:`windowEnergyLimit`

* :ref:`windowPopulationLimit`

* :ref:`maxIterationsBTE`

* :ref:`convergenceThresholdBTE`

* :ref:`dimensionality`

* :ref:`constantRelaxationTime`

* :ref:`withIsotopeScattering`

* :ref:`massVariance`

* :ref:`boundaryLength`

* :ref:`useSymmetries`

* :ref:`phonopyDispFileName`

* :ref:`dispFCFileName`



.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "phononTransport"
  phD2FileName = "./ForceConstants2nd"
  sumRuleD2 = "crystal"
  phD3FileName = "./ForceConstants3rd"
  qMesh = [10,10,10]
  temperatures = [300.]
  smearingMethod = "adaptiveGaussian"
  solverBTE = ["variational","relaxons"]
  scatteringMatrixInMemory = true
  boundaryLength = 10. mum

-------------------------------------

Electron BTE Solver
-----------------------------------

**Functionality:** Build and solve the electronic Boltzmann Transport Equation (BTE) using Wannier interpolation. Output quantites are electrical conductivity, electronic thermal conductivity, Seebeck coefficient, electron viscosity, and electronic lifetimes (on the kMesh used in this calculation).

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronWannierTransport"

* :ref:`phD2FileName`

* :ref:`sumRuleD2`

* :ref:`electronH0Name`

* :ref:`elphFileName`

* :ref:`kMesh`

* :ref:`temperatures`

* :ref:`minTemperature`

* :ref:`maxTemperature`

* :ref:`deltaTemperature`

* :ref:`dopings`

* :ref:`chemicalPotentials`

* :ref:`minChemicalPotential`

* :ref:`maxChemicalPotential`

* :ref:`deltaChemicalPotential`

* :ref:`smearingMethod`

* :ref:`smearingWidth`

* :ref:`windowType`

* :ref:`dimensionality`

* :ref:`constantRelaxationTime`

* :ref:`convergenceThresholdBTE`

* :ref:`maxIterationsBTE`

* :ref:`windowType`

* :ref:`windowEnergyLimit`

* :ref:`windowPopulationLimit`

* :ref:`solverBTE`

* :ref:`scatteringMatrixInMemory`

* :ref:`fermiLevel`

* :ref:`numOccupiedStates`

* :ref:`useSymmetries`

.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "electronWannierTransport"
  phD2FileName = "./silicon.fc"
  sumRuleD2 = "crystal"
  electronH0Name = "./si_tb.dat",
  elphFileName = "silicon.phoebe.elph.dat"
  kMesh = [15,15,15]
  temperatures = [300.]
  dopings = [1.e21]
  smearingMethod = "gaussian"
  smearingWidth = 0.5 eV
  windowType = "population"


-------------------------------------

EPA Transport Calculation
-----------------------------------

**Functionality:** Build and solve the electronic Boltzmann Transport Equation (BTE) using Wannier interpolation. Output quantites are electrical conductivity, electronic thermal conductivity, Seebeck coefficient, electron viscosity, and electronic lifetimes as a function of bin energy.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "transportEPA"

* :ref:`electronH0Name`

* :ref:`epaFileName`

* :ref:`electronFourierCutoff`

* :ref:`epaEnergyStep`

* :ref:`epaEnergyRange`

* :ref:`kMesh`

* :ref:`temperatures`

* :ref:`dopings`


.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "transportEpa"

  electronH0Name = "./out/silicon.xml",
  epaFileName = "./silicon.phoebe.epa.dat"

  electronFourierCutoff = 4.
  epaEnergyStep = 0.01 eV
  epaEnergyRange = 3.0 eV

  kMesh = [10,10,10]
  temperatures = [300.]
  dopings = [1.0e21]

-----------------------------------

Phonon Lifetimes on a Path
--------------------------

**Functionality:** Compute phonon lifetimes/linewidths on a specified path through the Brillouin zone.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "phononLifetimes"

* :ref:`phD2FileName`

* :ref:`sumRuleD2`

* :ref:`phD3FileName`

* :ref:`qMesh`

* :ref:`temperatures`

* :ref:`minTemperature`

* :ref:`maxTemperature`

* :ref:`deltaTemperature`

* :ref:`smearingMethod`

* :ref:`smearingWidth`

* :ref:`constantRelaxationTime`

* :ref:`withIsotopeScattering`

* :ref:`massVariance`

* :ref:`boundaryLength`

* :ref:`deltaPath`

* :ref:`beginEndPointPath`

* :ref:`phonopyDispFileName`

* :ref:`dispFCFileName`


.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "phononLifetimes"
  phD2FileName = "../Silicon/espresso.ifc2",
  sumRuleD2 = "simple"
  phD3FileName = "../Silicon/ShengBTEForceConstants3rd"
  qMesh = [15,15,15]
  temperatures = [600.]
  smearingMethod = "gaussian"
  smearingWidth = 10. cmm1
  deltaPath = 0.01
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path


-----------------------------------

Electron Lifetimes on a Path
-----------------------------

**Functionality:** Compute electron-phonon lifetimes/linewidths on a path through the Brillouin zone using Wannier interpolation.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronLifetimes"

* :ref:`phD2FileName`

* :ref:`electronH0Name`

* :ref:`sumRuleD2`

* :ref:`elphFileName`

* :ref:`kMesh`

* :ref:`temperatures`

* :ref:`minTemperature`

* :ref:`maxTemperature`

* :ref:`deltaTemperature`

* :ref:`dopings`

* :ref:`chemicalPotentials`

* :ref:`minChemicalPotential`

* :ref:`maxChemicalPotential`

* :ref:`deltaChemicalPotential`

* :ref:`smearingMethod`

* :ref:`smearingWidth`

* :ref:`constantRelaxationTime`

* :ref:`numOccupiedStates`

* :ref:`fermiLevel`

* :ref:`deltaPath`

* :ref:`beginEndPointPath`


.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "electronLifetimes"
  phD2FileName = "./silicon.fc",
  sumRuleD2 = "crystal"
  electronH0Name = "./si_tb.dat",
  elphFileName = "./silicon.phoebe.elph.dat"
  kMesh = [15,15,15]
  temperatures = [600.]
  dopings = [1.e22]
  smearingMethod = "gaussian"
  smearingWidth = 0.5 eV
  deltaPath = 0.01
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path

-----------------------------------

Phonon Dos
----------

**Functionality:** Compute the phonon density of states.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "phononDos"

* :ref:`phD2FileName`

* :ref:`sumRuleD2`

* :ref:`qMesh`

* :ref:`dosMinEnergy`

* :ref:`dosMaxEnergy`

* :ref:`dosDeltaEnergy`

* :ref:`phonopyDispFileName`

* :ref:`dispFCFileName`

.. raw:: html

  <h3>Sample input file</h3>

::

  phD2FileName = "qespresso/silicon.fc",
  sumRuleD2 = "simple"
  qMesh = [10,10,10]
  appName = "phononDos"
  dosMinEnergy = 0. cmm1
  dosMaxEnergy = 600. cmm1
  dosDeltaEnergy = 0.5 cmm1

--------------

Phonon Bands
------------

**Functionality:** Compute the phonon band structure on a path through the Brillouin zone.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "phononBands"

* :ref:`phD2FileName`

* :ref:`sumRuleD2`

* :ref:`deltaPath`

* :ref:`beginEndPointPath`

* :ref:`phonopyDispFileName`

* :ref:`dispFCFileName`


.. raw:: html

  <h3>Sample input file</h3>

::

  phD2FileName = "qespresso/silicon.fc",
  sumRuleD2 = "simple"
  appName = "phononBands"
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path

-----------------------------------

Electron DoS
------------

Electron DoS (Wannier interpolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Functionality:** Compute the electronic density of states. Electronic bands are interpolated to a finer mesh using Wannier interpolation.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronWannierDos"

* :ref:`electronH0Name`

* :ref:`fermiLevel`

* :ref:`kMesh`

* :ref:`dosMinEnergy`

* :ref:`dosMaxEnergy`

* :ref:`dosDeltaEnergy`

* :ref:`beginEndCrystal`

.. raw:: html

  <h3>Sample input file</h3>

::

  electronH0Name = "qespresso/si_tb.dat",
  kMesh = [10,10,10]
  appName = "electronWannierDos"
  dosMinEnergy = -6. eV
  dosMaxEnergy = 20. eV
  dosDeltaEnergy = 0.1 eV
  begin crystal
   Si 0.00000   0.00000   0.00000
   Si 1.34940   1.34940   1.34940
  end crystal

-----------------------------------

Electron DoS (Fourier interpolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Functionality:** Compute the electronic density of states. Electronic bands are interpolated to a finer mesh using Fourier interpolation.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronFourierDos"

* :ref:`electronH0Name`

* :ref:`kMesh`

* :ref:`fermiLevel`

* :ref:`dosMinEnergy`

* :ref:`dosMaxEnergy`

* :ref:`dosDeltaEnergy`

* :ref:`electronFourierCutoff`


.. raw:: html

  <h3>Sample input file</h3>

::

  electronH0Name = "qespresso/out/silicon.xml",
  kMesh = [10,10,10]
  appName = "electronFourierDos"
  dosMinEnergy = -6. eV
  dosMaxEnergy = 20. eV
  dosDeltaEnergy = 0.1 eV
  electronFourierCutoff = 4.

----------------------------------

Electron Bands
-----------------------------------

Electron Bands (Wannier interpolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Functionality:** Compute the phonon band structure on a path through the Brillouin zone. Electronic bands are interpolated to a finer mesh using Wannier interpolation.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronWannierBands"

* :ref:`electronH0Name`

* :ref:`fermiLevel`

* :ref:`deltaPath`

* :ref:`beginEndPointPath`

* :ref:`beginEndCrystal`


.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "electronWannierBands"
  electronH0Name = "qespresso/si_tb.dat",
  deltaPath = 0.01
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path
  begin crystal
   Si 0.00000   0.00000   0.00000
   Si 1.34940   1.34940   1.34940
  end crystal

-----------------------------------

Electron Bands (Fourier interpolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Functionality:** Compute the electronic band structure on a path through the Brillouin zone. Electronic bands are interpolated to a finer mesh using Fourier interpolation.

.. raw:: html

  <h3>Input variables</h3>


:ref:`appName` = "electronFourierBands"

* :ref:`electronH0Name`

* :ref:`fermiLevel`

* :ref:`deltaPath`

* :ref:`electronFourierCutoff`

* :ref:`beginEndPointPath`

.. raw:: html

  <h3>Sample input file</h3>

::

  appName = "electronFourierBands"
  electronH0Name = "qespresso/out/silicon.xml",
  deltaPath = 0.01
  electronFourierCutoff = 4.
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path



-----------------------------------

Input Variable Descriptions
----------------------------

Below are descriptions of the input variables available in Phoebe, along with their formats and when applicable, default values. If no default value is listed, this means the default value is nothing. If a required variable is left blank, Phoebe will throw a related error message.

.. _appName:

appName
^^^^^^^

* **Description:** This parameter, which must always be present, identifies which app (functionality) you want to run.

* **Format:** *string*

* **Required:** yes

**Possible values:**
  * "elPhQeToPhoebe": app to convert electron-phonon coupling from QE to Phoebe format (must be run before running any electron transport).

  * "phononTransport": app to solve the phonon BTE and compute phonon transport properties.

  * "electronWannierTransport": app to solve the electron BTE and compute electron transport properties.

  * "phononLifetimes": app to compute the phonon lifetimes on a given Brillouin zone path.

  * "electronLifetimes": app to compute the electron lifetimes on a given Brillouin zone path.

  * "phononDos": app to compute the phonon density of states.

  * "phononBands": app to compute the phonon bands on a Brillouin zone path.

  * "electronWannierDos": app to compute the electron density of states with Wannier interpolation.

  * "electronFourierDos": app to compute the electron density of states with Fourier interpolation.

  * "electronWannierBands": app to compute the electron bands with Wannier interpolation on a Brillouin zone path.

  * "electronFourierBands": app to compute the electron bands with Fourier interpolation on a Brillouin zone path.


.. _phD2FileName:

phD2FileName
^^^^^^^^^^^^

* **Description:** Path to a file containing harmonic force constants. File formats supported are: Quantum-ESPRESSO output of ``q2r.x`` (``prefix.fc``) or phono3py output (``fc2.hdf5``).

* **Format:** *string*

* **Required:** yes (for all phonon apps)


.. _phD3FileName:

phD3FileName
^^^^^^^^^^^^

* **Description:** Path to a file containing anharmonic (3rd order) force constants. File formats supported are: ShengBTE (FORCE_CONSTANTS_3RD) or phono3py (fc3.hdf5).

* **Format:** *string*

* **Required:** yes (for phonon transport and lifetime apps)

.. _phonopyDispFileName:

phonopyDispFileName
^^^^^^^^^^^^^^^^^^^

* **Description:** Path to the phono3py_disp.yaml file output by phono3py. (In the case of running only the harmonic phonons with phonopy, this file is named phonopy_disp.yaml).

* **Format:** *string*

* **Required:** yes (for calculations using phono3py)

.. _dispFCFileName:

dispFCFileName
^^^^^^^^^^^^^^

* **Description:** Path to the disp_fc3.yaml file output by phono3py. (In the case of running only the harmonic phonons with phonopy, this file is named disp.yaml.

* **Format:** *string*

* **Required:** yes (for calculations using phono3py)


.. _sumRuleD2:

sumRuleD2
^^^^^^^^^

* **Description:** If specified, applies an acoustic sum rule to the phonon harmonic force constants. Allowed values are "simple" or "crystal", with the same algorithm and meaning of Quantum-ESPRESSO ``matdyn.x`` program.

* **Format:** *string*

* **Required:** yes (for phonon calculations)


.. _qMesh:

qMesh
^^^^^

* **Description:** Triplet of integers with the fine q-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of phonon properties.

* **Format:** *list of int*

* **Required:** yes (for phonon calculations)


.. _kMesh:

kMesh
^^^^^

* **Description:** Triplet of integers with the fine k-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of electronic properties. In electron-phonon transport calculations, `qMesh` is set to be equal to this value and does not need to be specified by the user.

* **Format:** *list of int*

* **Required:** yes (for electronic calculations)


.. _temperatures:

temperatures
^^^^^^^^^^^^

* **Description:** List with the values of temperatures (in Kelvin) to be used in the calculation. If scatteringMatrixInMemory=true, only one value of temperature is allowed.

* **Format:** *list of doubles*

* **Required:** yes (for transport and lifetime calculations)


.. _smearingMethod:

smearingMethod
^^^^^^^^^^^^^^

* **Description:** Selects the level of approximation for replacing the Dirac-delta approximating energy conservation. Allowed values are "gaussian" and "adaptiveGaussian" (preferred)

* **Format:** *string*

* **Required:** yes (for transport and lifetime calculations)


.. _smearingWidth:

smearingWidth
^^^^^^^^^^^^^

* **Description:** This parameter is required if :ref:`smearingMethod` = "gaussian", where this parameter represents the full-width half-maximum of the Gaussian used to approximate the Dirac-delta conserving energy. Example: smearingWidth = 0.5 eV

* **Format:** *double+units*

* **Required:** yes (when :ref:`smearingMethod` = "gaussian")


.. _solverBTE:

solverBTE
^^^^^^^^^

* **Description:** If specified, solves the Boltzmann equation beyond the relaxation time approximation. Allowed values are: "variational", "iterative", and "relaxons", see the Theory section for a detailed explanation. Example: solverBTE=["variational","relaxons"]

* **Format:** *list of strings*

* **Required:** no


.. _scatteringMatrixInMemory:

scatteringMatrixInMemory
^^^^^^^^^^^^^^^^^^^^^^^^

* **Description:** If true, the scattering matrix is kept in memory, and only one temperature is allowed. In exchange for a larger memory usage, exact BTE solvers are much faster. Disable this flag to reduce the memory footprint, at the cost of slowing down the exact BTE solvers.

* **Format:** *bool*

* **Required:** no

* **Default:** `true`


.. _windowType:

windowType
^^^^^^^^^^

* **Description:** Enables the window used to discard some phonon states that don't contribute to transport. Possible values are "nothing", "population" and "energy". "nothing" means window is not applied; "population" means phonon states are discarded if :math:`\frac{\partial \bar{n}}{\partial T} <` windowPopulationLimit, where :math:`\frac{\partial \bar{n}}{\partial T}` is the Bose--Einstein distribution derivative.

* **Format:** *string*

* **Required:** no

* **Default:** `"nothing"`

.. _windowEnergyLimit:

windowEnergyLimit
^^^^^^^^^^^^^^^^^

* **Description:** Additional parameter for energy :ref:`windowType`. Specify two values :math:`E_{min}` and :math:`E_{max}` (in electronVolts) such that we discard all phonon states  with energy outside of these bounds.

* **Format:** *list of doubles*

* **Required:** yes (if :ref:`windowType` = "energy")


.. _windowPopulationLimit:

windowPopulationLimit
^^^^^^^^^^^^^^^^^^^^^

* **Description:** Required if :ref:`windowType` = "population". Cutoff values for discarding phonon states based on their equilibrium phonon occupation number, such that :math:`\frac{\partial \bar{n}}{\partial T} <` windowPopulationLimit.

* **Format:** *double*

* **Required:** yes (if :ref:`windowType` = "population")


.. _maxIterationsBTE:

maxIterationsBTE
^^^^^^^^^^^^^^^^

* **Description:** Maximum number of iterations for iterative and variational BTE solvers. If the maximum number of iterations is reached, the code will throw an error.

* **Format:** *int*

* **Required:** no

* **Default:** `50`


.. _convergenceThresholdBTE:

convergenceThresholdBTE
^^^^^^^^^^^^^^^^^^^^^^^

* **Description:** Convergence criterion to stop iterative BTE solvers. The calculation is converged if the transport coefficients have a relative change smaller than convergenceThresholdBTE.

* **Format:** *double*

* **Required:** no

* **Default:** `1.0e-5`


.. _dimensionality:

dimensionality
^^^^^^^^^^^^^^

* **Description:** Input the dimensionality of the material. As a result, transport coefficients tensors will be of size (dim x dim), and units will be suitably scaled for the desired dimensionality.

* **Format:** *int*

* **Required:** no

* **Default:** `3`


.. _constantRelaxationTime:

constantRelaxationTime
^^^^^^^^^^^^^^^^^^^^^^

* **Description:** If specified, we solve the BTE with the constant relaxation time approximation, where the electron or phonon lifetime is set to this input value. (Fast but inaccurate!)

* **Format:** *double+units*

* **Required:** no


.. _withIsotopeScattering:

withIsotopeScattering
^^^^^^^^^^^^^^^^^^^^^

* **Description:** Controls whether to include or not phonon-isotope scattering

* **Format:** *bool*

* **Required:** no

* **Default:** `true`


.. _massVariance:

massVariance
^^^^^^^^^^^^

* **Description:** User can specify a list of custom atomic mass variances :math:`g_2^s`. See Theory section for a description. The mass variances must be ordered in the same way that atomic species are specified in the file :ref:`phD2FileName`. Defaults to the mass variance for natural isotopic abundance.

* **Format:** *list of doubles*

* **Required:** no

* **Default:** natural isotopic abundance


.. _boundaryLength:

boundaryLength
^^^^^^^^^^^^^^

* **Description:** If specified, includes the phonon-boundary scattering within the RTA approximation. Example: boundaryLength = 10 mum

* **Format:** *double+units*

* **Required:** no


.. _electronH0Name:

electronH0Name
^^^^^^^^^^^^^^

* **Description:** For Wannier-interpolation-based calculations, electronH0Name must contain the path to the ``{prefix}_tb.dat`` file generated by Wannier90. For Fourier-interpolation-based calculations, electronH0Name must contain the path to the Quantum-ESPRESSO ``{outdir}/{prefix}.xml`` file generated by ``pw.x``.

* **Format:** *string*

* **Required:** yes


.. _dosMinEnergy:

dosMinEnergy
^^^^^^^^^^^^

* **Description:** Used in conjunction with :ref:`dosMaxEnergy` and :ref:`dosDeltaEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.

* **Format:** *double+units*

* **Required:** yes (for dos calculation)


.. _dosMaxEnergy:

dosMaxEnergy
^^^^^^^^^^^^

* **Description:** Used in conjunction with :ref:`dosMinEnergy` and :ref:`dosDeltaEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.

* **Format:** *double+units*

* **Required:** yes (for dos calculation)


.. _dosDeltaEnergy:

dosDeltaEnergy
^^^^^^^^^^^^^^

* **Description:** Used in conjunction with :ref:`dosMinEnergy` and :ref:`dosMaxEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.

* **Format:** *double+units*

* **Required:** yes (for dos calculation)


.. _electronFourierCutoff:

electronFourierCutoff
^^^^^^^^^^^^^^^^^^^^^

* **Description:** A parameter controlling the search of lattice vectors used for the Fourier interpolation of the electronic band structure. In detail, the lattice vectors used for the Fourier interpolation are searched in a supercell of size `electronFourierCutoff`:sup:`3` the primitive unit cell. Set it to at least 2.

* **Format:** *double*

* **Required:** yes (for electron Fourier DoS and bands)


.. _beginEndCrystal:

begin/end crystal
^^^^^^^^^^^^^^^^^

* **Description:** Specify the atomic species and atomic positions inside the crystal. This needs to be specified in some apps like WannierBands or WannierDos, as the output files of Wannier90 doesn't provide all the information about the crystal.

* **Format:** Namelist format, with atomic symbol and position coordinates in units of Angstroms. Example::

    begin crystal
     Si 0.00000   0.00000   0.00000
     Si 1.34940   1.34940   1.34940
    end crystal

* **Required:** yes (for electron Wannier DoS and bands calculations)


.. _beginEndPointPath:

begin/end point path
^^^^^^^^^^^^^^^^^^^^

* **Description:** Specify the path of wavectors in the Brillouin zone used in apps such as `phononBands` or `phononLifetimes`. Use the parameter :ref:`deltaPath` to control the number of wavevectors in each segment.

* **Format:** Namelist format, as pairs of special point symbol and wavevector coordinates. Wavevector coordinates are in fractional coordinates with respect to the primitive reciprocal lattice vectors. Example::

    begin point path
     L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
     G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
     X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
     K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
    end point path

* **Required:** yes (for electron and phonon bands and lifetime apps)

.. _dopings:

dopings
^^^^^^^

* **Description:** Specify a list of doping concentrations, in cm :sup:`-3`, to compute electronic properties at various doping concentrations. The chemical potentials corresponding to this doping concentrations will be computed.

* **Format:** *list of doubles*

* **Required:** yes (for electron transport and lifetime calculations, unless :ref:`chemicalPotentials` is specified)


.. _chemicalPotentials:

chemicalPotentials
^^^^^^^^^^^^^^^^^^

* **Description:** Specify a list of chemical potentials to be used for the calculation of properties as a function of the chemical potential. If used in electron Wannier transport and scatteringMatrixInMemory=true, then only one value of chemical potentials can be specified. Values are in eV.

* **Format:** *list of doubles*

* **Required:** yes (unless `minChemicalPotential, maxChemicalPotential, deltaChemicalPotential` variables are present, or :ref:`dopings` are specified).


.. _elphFileName:

elphFileName
^^^^^^^^^^^^

* **Description:** Path to the file generated by the app `elPhQeToPhoebe` containing the electron-phonon coupling in the Wannier representation (e.g. ``{prefix}.phoebe.elph.dat``)

* **Format:** *string*

* **Required:** yes (for electron transport and lifetime apps)


.. _epaMinEnergy:

epaMinEnergy
^^^^^^^^^^^^

* **Description:** Specifies the minimum the energy bin value over which the electron-phonon coupling will be averaged for the EPA approximation post-processing in the ``elPhQeToPhoebe`` app.

* **Format:** *double*

* **Required:** yes (for EPA ``elPhQeToPhoebe`` runs)


.. _epaMaxEnergy:

epaMaxEnergy
^^^^^^^^^^^^

* **Description:** Specifies the maximum the energy bin value over which the electron-phonon coupling will be averaged for the EPA approximation post-processing in the ``elPhQeToPhoebe`` app.

* **Format:** *double*

* **Required:** yes (for EPA ``elPhQeToPhoebe`` runs)


.. _epaNumBins:

epaNumBins
^^^^^^^^^^^

* **Description:** The number of energy bins, ranging from :ref:`epaMinEnergy` to :ref:`epaMaxEnergy` used to average the electron-phonon matrix elements for EPA calculations in the ``elPhQeToPhoebe`` app.

* **Format:** *int*

* **Required:** yes (for EPA ``elPhQeToPhoebe`` runs)


.. _epaSmearingEnergy:

epaSmearingEnergy
^^^^^^^^^^^^^^^^^

* **Description:** Specifies the Gaussian width used in the moving least squares averaging procedure used to average the electron-phonon matrix elements for an EPA calculation.

* **Format:** *double*

* **Required:** yes (for EPA ``elPhQeToPhoebe`` runs)


.. _epaFileName:

epaFileName
^^^^^^^^^^^

* **Description:** This is the path to the file ``*.phoebe.epa.dat``, which is created by ``elPhQeToPhoebe``.

* **Format:** *string*

* **Required:** yes (for EPA transport app)


.. _epaEnergyStep:

epaEnergyStep
^^^^^^^^^^^^^

* **Description:** The energy interval used to integrate the transport coefficients, i.e. lifetimes will be computed every ``epaEnergyStep`` energies.

* **Format:** *double*

* **Required:** yes (for EPA transport app)


.. _epaEnergyRange:

epaEnergyRange
^^^^^^^^^^^^^^

* **Description:** EPA lifetimes will be computed for all energies in proximity of the chemical potential, i.e. for all energies such that :math:`|\epsilon-\mu|<\text{epaEnergyRange}`. This variable specifies that range.

* **Format:** *double*

* **Required:** yes (for EPA transport app)


.. _deltaPath:

deltaPath
^^^^^^^^^

* **Description:** This variable controls how far apart are the wavevectors when a path in the Brillouin zone is specified, and it represents the distance (in Bohr) between wavevectors. Can be used when a path of wavevectors is specified with the :ref:`beginEndPointPath` key.

* Default: 0.05 Bohr:sup:`-1`

* **Format:** *string*

* **Required:** no


.. _elPhInterpolation:

elPhInterpolation
^^^^^^^^^^^^^^^^^

* **Description:** Can be either "wannier" or "epa". The first prepares the electron-phonon coupling for the transport calculation with Wannier interpolation (i.e. does the transformation from Bloch to Wannier representation). The second, prepares the electron-phonon coupling to be used with the EPA approximation.

* **Format:** *string*

* **Required:** yes (for elPhQeToPhoebe app)


.. _wannier90Prefix:

wannier90Prefix
^^^^^^^^^^^^^^^

* **Description:** Set to the same value of ``prefix`` in Wannier90. It's used to locate the files ``{prefix}.eig`` generated by ``wannier90.x``.

* **Format:** *string*

* **Required:** yes (for any electron app using Wannier interpolation)


.. _quantumEspressoPrefix:

quantumEspressoPrefix
^^^^^^^^^^^^^^^^^^^^^

* **Description:** Set to the same value of ``prefix`` in Quantum-ESPRESSO. It's used to locate the files ``{prefix}.dyn*`` or ``{prefix}.phoebe.*.dat`` generated by ``ph.x``.

* **Format:** *string*

* **Required:** yes


.. _fermiLevel:

fermiLevel
^^^^^^^^^^

* **Description:** Sets the fermi level of the ground state. Can be specified e.g. in Bands or DOS apps to specify an otherwise unknown fermi level. This quantity is read from file for transport calculations: this input parameter overwrites that value, use with caution.

* **Format:** *double+units*

* **Required:** no


.. _numOccupiedStates:

numOccupiedStates
^^^^^^^^^^^^^^^^^

* **Description:** Determines the number of occupied Kohn-Sham states at the ground state. The default value might be read from the :ref:`electronH0Name` (when this is the Quantum-ESPRESSO xml file) or the file with the el-ph interaction (so, the user may not need to specify it for transport calculations). This value controls where the Fermi level is set. The user alternatively can specify the :ref:`fermiLevel` (and :ref:`numOccupiedStates` will be computed from the Fermi level).

* **Format:** *double*

* **Required:** no


.. _minChemicalPotential:

minChemicalPotential
^^^^^^^^^^^^^^^^^^^^

* **Description:** To be used together with :ref:`maxChemicalPotential` and :ref:`deltaChemicalPotential`, sets the code to compute properties at all chemical Potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with :ref:`chemicalPotentials` to instead manually specify the chemical potentials of the calculation.

* **Format:** *double*

* **Required:** yes (either specify :ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, and :ref:`deltaChemicalPotential` **or** :ref:`chemicalPotentials`.


.. _maxChemicalPotential:

maxChemicalPotential
^^^^^^^^^^^^^^^^^^^^

* **Description:** To be used together with :ref:`minChemicalPotential` and :ref:`deltaChemicalPotential`, sets the code to compute properties at all chemical potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with chemicalPotentials to instead manually specify the chemical potentials of the calculation.

* **Format:** *double*

* **Required:** yes (either specify :ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, and :ref:`deltaChemicalPotential` **or** :ref:`chemicalPotentials`)


.. _deltaChemicalPotential:

deltaChemicalPotential
^^^^^^^^^^^^^^^^^^^^^^

* **Description:** To be used together with :ref:`minChemicalPotential` and :ref:`maxChemicalPotential`, sets the code to compute properties at all chemical Potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with :ref:`chemicalPotentials` to instead manually specify the chemical potentials of the calculation.

* **Format:** *double*

* **Required:** yes (either specify :ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, and :ref:`deltaChemicalPotential` **or** :ref:`chemicalPotentials`)


.. _minTemperature:

minTemperature
^^^^^^^^^^^^^^

* **Description:** To be used together with :ref:`maxTemperature` and :ref:`deltaTemperature`, sets the code to compute observables at temperatures (in Kelvin) between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.

* **Format:** *double*

* **Required:** yes (either set :ref:`minTemperature`, :ref:`maxTemperature`, and :ref:`deltaTemperature` **or** :ref:`temperatures`)


.. _maxTemperature:

maxTemperature
^^^^^^^^^^^^^^

* **Description:** To be used together with :ref:`minTemperature` and :ref:`deltaTemperature`, sets the code to compute observables at temperatures (in Kelvin) between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.

* **Format:** *double*

* **Required:** yes (either set :ref:`minTemperature`, :ref:`maxTemperature`, and :ref:`deltaTemperature` **or** :ref:`temperatures`)


.. _deltaTemperature:

deltaTemperature
^^^^^^^^^^^^^^^^

* **Description:** To be used together with minTemperature and maxTemperature, sets the code to compute observables at temperatures (in Kelvin) between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.

* **Format:** *double*

* **Required:** yes (either set :ref:`minTemperature`, :ref:`maxTemperature`, and :ref:`deltaTemperature` **or** :ref:`temperatures`)


.. _useSymmetries:

useSymmetries
^^^^^^^^^^^^^

* **Description:** When set to true, triggers the BTE to be computed only on the irreducible wedge of the Brillouin zone. For systems with several symmetries, this speeds up calculations. On the other hand, it may slow down the code for systems with few symmetries. If set to true, the viscosity is only computed at the RTA level.

* **Format:** *bool*

* **Required:** no

* **Default:** `false`

