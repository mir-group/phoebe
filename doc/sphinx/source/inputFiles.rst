Input files overview
====================

Here we provide a description of the input variables of the various apps in the code.
Note that the text in the user input file is case sensitive!
In the code, we try to check the variables provided in input, and stop the code if we detect missing or invalid input parameters.
Note that the executable is the same for all this applications; the applications are selected through the :ref:`appName` parameter.




Phonon BTE
----------

Target: build and solve the phonon Boltzmann Transport Equation (BTE), and compute phonon transport properties such as thermal conductivity, relaxation times and viscosity.

Input variables:

* :ref:`appName` = "phononTransport"

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


Sample input file::

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




Electron BTE, Wannier interpolation
-----------------------------------

Target: build and solve the electronic Boltzmann Transport Equation (BTE), using Wannier-based interpolation. Output quantites are electrical conductivity, electronic thermal conductivity, Seebeck coefficient, electron viscosity and electronic lifetimes.  

Input variables:

* :ref:`appName` = "electronWannierTransport"
  
* :ref:`phD2FileName`
  
* :ref:`sumRuleD2`
  
* :ref:`electronH0Name`
  
* :ref:`epwFileName`
  
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

Sample input file::

  appName = "electronWannierTransport"
  phD2FileName = "./silicon.fc"
  sumRuleD2 = "crystal"
  electronH0Name = "./si_tb.dat",
  epwFileName = "silicon.phoebe.elph.dat"
  kMesh = [15,15,15]
  temperatures = [300.]
  dopings = [1.e21]
  smearingMethod = "gaussian"
  smearingWidth = 0.5 eV
  windowType = "population"




QE to Phoebe
------------

Target: convert the electron-phonon coupling from a Quantum-ESPRESSO format to a Phoebe format.
In doing so, we postprocess data for the Wannier or EPA interpolations.

Input variables:

* :ref:`appName` = "elPhQeToPhoebe"
  
* :ref:`phD2FileName`
  
* :ref:`elPhInterpolation`
  
* :ref:`electronH0Name`
  
* :ref:`wannier90Prefix`
  
* :ref:`quantumEspressoPrefix`


Sample input file::

  appName = "elPhQeToPhoebe"
  elPhInterpolation = "wannier"
  phD2FileName = "./silicon.fc"
  electronH0Name = "./si_tb.dat",
  wannier90Prefix = "si"
  quantumEspressoPrefix = "silicon"



Phonon Lifetimes
----------------

Target: compute phonon lifetimes/linewidths on a path in the Brillouin zone.

Input variables:

* :ref:`appName` = "phononLifetimes"
  
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


Sample input file::

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



Electron Wannier Lifetimes
--------------------------

Target: compute electron lifetimes/linewidths on a path in the Brillouin zone, using Wannier-based interpolations.

Input variables:

* :ref:`appName` = "electronLifetimes"
  
* :ref:`phD2FileName`
  
* :ref:`electronH0Name`
  
* :ref:`sumRuleD2`
  
* :ref:`epwFileName`
  
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

  
Sample input file::

  appName = "electronLifetimes"
  phD2FileName = "./silicon.fc",
  sumRuleD2 = "crystal"
  electronH0Name = "./si_tb.dat",
  epwFileName = "./silicon.phoebe.elph.dat"
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


Phonon Dos
----------

Target: compute the phonon Density of States.

Input variables:

* :ref:`appName` = "phononDos"
  
* :ref:`phD2FileName`
  
* :ref:`sumRuleD2`
  
* :ref:`qMesh`
  
* :ref:`dosMinEnergy`
  
* :ref:`dosMaxEnergy`
  
* :ref:`dosDeltaEnergy`

Sample input file::

  phD2FileName = "qespresso/silicon.fc",
  sumRuleD2 = "simple"
  qMesh = [10,10,10]
  appName = "phononDos"
  dosMinEnergy = 0. cmm1
  dosMaxEnergy = 600. cmm1
  dosDeltaEnergy = 0.5 cmm1




Electron DoS, Wannier interpolation
-----------------------------------

Target: compute the electronic Density of States. Electronic bands are interpolated to a finer mesh using maximally localized Wannier function interpolation.

Input variables:

* :ref:`appName` = "electronWannierDos"
  
* :ref:`electronH0Name`
  
* :ref:`fermiLevel`
  
* :ref:`kMesh`
  
* :ref:`dosMinEnergy`
  
* :ref:`dosMaxEnergy`
  
* :ref:`dosDeltaEnergy`
  
* :ref:`beginEndCrystal`

Sample input file::

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



Electron DoS, Fourier interpolation
-----------------------------------

Target: compute the electronic Density of States. Electronic bands are interpolated to finer meshes using a Fourier interpolation.

Input variables:

* :ref:`appName` = "electronFourierDos"
  
* :ref:`electronH0Name`
  
* :ref:`kMesh`
  
* :ref:`fermiLevel`
  
* :ref:`dosMinEnergy`
  
* :ref:`dosMaxEnergy`
  
* :ref:`dosDeltaEnergy`
  
* :ref:`electronFourierCutoff`

Sample input file::

  electronH0Name = "qespresso/out/silicon.xml",
  kMesh = [10,10,10]
  appName = "electronFourierDos"
  dosMinEnergy = -6. eV
  dosMaxEnergy = 20. eV
  dosDeltaEnergy = 0.1 eV
  electronFourierCutoff = 4.



Phonon Bands
------------

Target: compute the phonon band structure on a path in the Brillouin zone.

Input variables:

* :ref:`appName` = "phononBands"
  
* :ref:`phD2FileName`
  
* :ref:`sumRuleD2`
  
* :ref:`deltaPath`
  
* :ref:`beginEndPointPath`

Sample input file::

  phD2FileName = "qespresso/silicon.fc",
  sumRuleD2 = "simple"
  appName = "phononBands"
  begin point path
   L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
   G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
   X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000 
   K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
  end point path



Electron Bands, Wannier interpolation
-------------------------------------

Target: compute the phonon band structure on a path in the Brillouin zone. Electronic bands are interpolated using Wannier functions.

Input variables:

* :ref:`appName` = "electronWannierBands"
  
* :ref:`electronH0Name`
  
* :ref:`fermiLevel`
  
* :ref:`deltaPath`
  
* :ref:`beginEndPointPath`
  
* :ref:`beginEndCrystal`


Sample input file::

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



Electron Bands, Fourier interpolation
-------------------------------------

Target: compute the electronic band structure on a path in the Brillouin zone. Electronic bands are interpolated using Wannier functions.

Input variables:

* :ref:`appName` = "electronFourierBands"
  
* :ref:`electronH0Name`
  
* :ref:`fermiLevel`
  
* :ref:`deltaPath`
  
* :ref:`electronFourierCutoff`
  
* :ref:`beginEndPointPath`

Sample input file::

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








Variable descriptions
---------------------


.. _appName:

appName
^^^^^^^

* This parameter, which must always be present, identifies which App (functionality) you want to run. Allowed values are:
  
  * "phononTransport": app to solve the phonon BTE and compute phonon transport properties.
    
  * "phononLifetimes": app to compute the phonon lifetimes.

  * "elPhQeToPhoebe": app to convert electron-phonon coupling from qe2Phoebe (must be run before running any electron Transport).
    
  * "electronWannierTransport": app to solve the electron BTE with Wannier interpolation.
    
  * "electronLifetimes": app to compute the electron lifetimes.
  
  * "phononDos": app to compute the phonon density of states.
    
  * "electronWannierDos": app to compute the electron density of states with Wannier interpolation.
    
  * "electronFourierDos": app to compute the electron Density of States with Fouerier interpolation.
  
  * "phononBands": app to compute the phonon bands on a path.
    
  * "electronWannierBands": app to compute the electron bands with Wannier interpolation on a path in the Brillouin zone.
    
  * "electronFourierBands": app to compute the electron bands with Fourier interpolation on a path in the Brillouin zone. 

* *string*
  
* Required


.. _phD2FileName:

phD2FileName
^^^^^^^^^^^^

* Path to the file with harmonic force constants. File format supported is Quantum-ESPRESSO output of q2r.x.

* *string*
  
* Required


.. _phD3FileName:

phD3FileName
^^^^^^^^^^^^

* Path to the file with anharmonic 3rd order force constants.
  
* *string*
  
* Required


.. _sumRuleD2:

sumRuleD2
^^^^^^^^^

* If specified, applies an acoustic sum rule to the phonon harmonic force constants. Allowed values are "simple" or "crystal", with the same algorithm and meaning of Quantum-ESPRESSO matdyn.x program.
  
* *string*
  
* Required


.. _qMesh:

qMesh
^^^^^

* Triplet of integers with the fine q-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of phonon properties.
  
* *list of int*
  
* Required


.. _kMesh:

kMesh
^^^^^

* Triplet of integers with the fine k-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of electronic properties. In electron-phonon transport calculations, `qMesh` is set to be equal to this value and does not need to be specified by the user.
  
* *list of int*
  
* Required


.. _temperatures:

temperatures
^^^^^^^^^^^^

* List with the values of temperatures to be used in the calculation. If scatteringMatrixInMemory=true, only one value of temperature is allowed.
  
* *list of doubles*
  
* Required


.. _smearingMethod:

smearingMethod
^^^^^^^^^^^^^^

* Selects the level of approximation for replacing the Dirac-delta approximating energy conservation. Allowed values are "gaussian" and "adaptiveGaussian" (preferred)
  
* *string*
  
* Required


.. _smearingWidth:

smearingWidth
^^^^^^^^^^^^^

* This parameter is required if :ref:`smearingMethod` = "gaussian", where this parameter represents the full width at half maximum of the gaussian used to approximate the Dirac-delta conserving energy. Example: smearingWidth = 0.5 eV
  
* *double+units*
  
* Required if :ref:`smearingMethod` = "gaussian"


.. _solverBTE:

solverBTE
^^^^^^^^^

* If specified, solves the Boltzmann equation beyond the relaxation time approximation. Allowed values are: "variational", "iterative", and "relaxons", see the Theory section for a detailed explanation. Example: solverBTE=["variational","relaxons"]
  
* *list of strings*
  
* Optional


.. _scatteringMatrixInMemory:

scatteringMatrixInMemory
^^^^^^^^^^^^^^^^^^^^^^^^

* If true, the scattering matrix is kept in memory, and only one temperature is allowed. In exchange for a larger memory usage, exact BTE solvers are much faster; disable this flag to reduce the memory footprint but slowing down the exact BTE solvers.
  
* *bool*
  
* Default=`true`
  
* Optional


.. _windowType:

windowType
^^^^^^^^^^

* Enables the window used to discard some phonon states that don't contribute to transport. Possible values are "nothing", "population" and "energy". "nothing" means window is not applied; "population" means phonon states are discarded if :math:`\frac{\partial \bar{n}}{\partial T} <` windowPopulationLimit, where :math:`\frac{\partial \bar{n}}{\partial T}` is the Bose--Einstein distribution derivative.
  
* *string*
  
* Default: `"nothing"`
  
* Optional


.. _windowEnergyLimit:

windowEnergyLimit
^^^^^^^^^^^^^^^^^

* Additional parameter for energy :ref:`windowType`. Specify two values :math:`E_{min}` and :math:`E_{max}` (in electronVolts) such that we discard all phonon states  with energy outside of these bounds.
  
* *list of doubles*
  
* Required if :ref:`windowType` = "energy"


.. _windowPopulationLimit:

windowPopulationLimit
^^^^^^^^^^^^^^^^^^^^^

* Required if :ref:`windowType` = "population". Cutoff values for discarding phonon states based on their equilibrium phonon occupation number, such that :math:`\frac{\partial \bar{n}}{\partial T} <` windowPopulationLimit.
  
* *double*
  
* Required if :ref:`windowType` = "population"


.. _maxIterationsBTE:

maxIterationsBTE
^^^^^^^^^^^^^^^^

* Maximum number of iterations for iterative and variational BTE solvers. If the maximum number of iterations is reached, the code will throw an error.
  
* *int*
  
* Default: `50`
  
* Optional


.. _convergenceThresholdBTE:

convergenceThresholdBTE
^^^^^^^^^^^^^^^^^^^^^^^

* Convergence criterion to stop iterative BTE solvers. The calculation is converged if the transport coefficients have a relative change smaller than convergenceThresholdBTE.
  
* *double*
  
* Default: `1.0e-5`
  
* Optional


.. _dimensionality:

dimensionality
^^^^^^^^^^^^^^

* Input the dimensionality of the material. As a result, transport coefficients tensors will be of size (dim x dim), and units will be suitably scaled for the desired dimensionality.
  
* *int*
  
* Default: `3`
  
* Optional


.. _constantRelaxationTime:

constantRelaxationTime
^^^^^^^^^^^^^^^^^^^^^^

* If specified, we solve the BTE with the constant relaxation time approximation, where the phonon lifetime is set to this input value. (Fast but inaccurate!)
  
* *double+units*
  
* Optional


.. _withIsotopeScattering:

withIsotopeScattering
^^^^^^^^^^^^^^^^^^^^^

* Controls whether to include or not phonon-isotope scattering
  
* *bool*
  
* Default: `true`
  
* Optional


.. _massVariance:

massVariance
^^^^^^^^^^^^

* User can specify a list of custom atomic mass variances :math:`g_2^s`. See Theory section for a description. The mass variances must be ordered in the same way that atomic species are specified in the file :ref:`phD2FileName`. Defaults to the mass variance for natural isotopic abundance.
  
* *list of doubles*

* Default: natural isotopic abundance
  
* Optional


.. _boundaryLength:

boundaryLength
^^^^^^^^^^^^^^

* If specified, includes the phonon-boundary scattering within the RTA approximation. Example: boundaryLength = 10 mum
  
* *double+units*
  
* Optional


.. _electronH0Name:

electronH0Name
^^^^^^^^^^^^^^

* For Wannier-interpolation-based calculations, `electronH0Name` must contain the path to the `{prefix}_tb.dat` file generated by Wannier90. For Fourier-interpolation-based calculations, `electronH0Name` must contain the path to the Quantum-ESPRESSO {outdir}/{prefix}.xml file generated by `pw.x`.
  
* *string*
  
* Required


.. _dosMinEnergy:

dosMinEnergy
^^^^^^^^^^^^

* Used in conjunction with :ref:`dosMaxEnergy` and :ref:`dosDeltaEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.
  
* *double+units*
  
* Required


.. _dosMaxEnergy:

dosMaxEnergy
^^^^^^^^^^^^

* Used in conjunction with :ref:`dosMinEnergy` and :ref:`dosDeltaEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.
  
* *double+units*
  
* Required


.. _dosDeltaEnergy:

dosDeltaEnergy
^^^^^^^^^^^^^^

* Used in conjunction with :ref:`dosMinEnergy` and :ref:`dosMaxEnergy` to compute the Density of States every :ref:`dosDeltaEnergy` increments between :ref:`dosMinEnergy` and :ref:`dosMaxEnergy`.
  
* *double+units*
  
* Required


.. _electronFourierCutoff:

electronFourierCutoff
^^^^^^^^^^^^^^^^^^^^^

* A parameter controlling the search of lattice vectors used for the Fourier interpolation of the electronic band structure. In detail, the lattice vectors used for the Fourier interpolation are searched in a supercell of size `electronFourierCutoff`:sup:`3` the primitive unit cell. Set it to at least 2.
  
* *double*
  
* Required


.. _beginEndCrystal:

begin/end crystal
^^^^^^^^^^^^^^^^^

* Specify the atomic species and atomic positions inside the crystal. This needs to be specified in some apps like WannierBands or WannierDos, as the output files of Wannier90 doesn't provide all the information about the crystal.
  
* Namelist format, with atomic symbol and position coordinates in units of Angstroms. Example::

    begin crystal
     Si 0.00000   0.00000   0.00000
     Si 1.34940   1.34940   1.34940
    end crystal

* Required


.. _beginEndPointPath:

begin/end point path
^^^^^^^^^^^^^^^^^^^^

* Specify the path of wavectors in the Brillouin zone used in apps such as `phononBands` or `phononLifetimes`. Use the parameter :ref:`deltaPath` to control the number of wavevectors in each segment.
  
* Namelist format, as pairs of special point symbol and wavevector coordinates. Wavevector coordinates are in fractional coordinates with respect to the primitive reciprocal lattice vectors. Example::

    begin point path
     L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
     G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
     X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000 
     K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
    end point path

* Required


.. _dopings:

dopings
^^^^^^^

* Specify a list of doping concentrations, in cm :sup:`-3`, to compute electronic properties at various doping concentrations. The chemical potentials corresponding to this doping concentrations will be computed.
  
* *list of doubles*
  
* Required; alternatively one must specify :ref:`chemicalPotentials`.


.. _chemicalPotentials:

chemicalPotentials
^^^^^^^^^^^^^^^^^^

* Specify a list of chemical potentials to be used for the calculation of properties as a function of the chemical potential. If used in electron Wannier transport and scatteringMatrixInMemory=true, then only one value of chemical potentials can be specified. Values are in eV.
  
* *list of doubles*
  
* Required. The user can substitute this parameter by specifying `(minChemicalPotential,maxChemicalPotential,deltaChemicalPotential)`.


.. _epwFileName:

epwFileName
^^^^^^^^^^^

* Path to the file generated by the app `elPhQeToPhoebe` containing the electron-phonon coupling in the Wannier representation (e.g. `{prefix}.phoebe.elph.dat`)
  
* *string*
  
* Required


.. _deltaPath:

deltaPath
^^^^^^^^^

* This variable controls how far apart are the wavevectors when a path in the Brillouin zone is specified, and it represents the distance (in Bohr) between wavevectors. Can be used when a path of wavevectors is specified with the :ref:`beginEndPointPath` key.
  
* Default: 0.05 Bohr:sup:`-1`
  
* *string*
  
* Optional


.. _elPhInterpolation:

elPhInterpolation
^^^^^^^^^^^^^^^^^

* Can be either "wannier" or "epa". The first, prepares the electron-phonon coupling for the transport calculation with Wannier interpolation (i.e. does the transformation from Bloch to Wannier representation). The second, prepares the electron-phonon coupling to be used with the EPA approximation.
  
* *string*
  
* Required


.. _wannier90Prefix:

wannier90Prefix
^^^^^^^^^^^^^^^

* Set to the same value of `prefix` in Wannier90. It's used to locate the files `{prefix}.eig` generated by `wannier90.x`.
  
* *string*
  
* Required


.. _quantumEspressoPrefix:

quantumEspressoPrefix
^^^^^^^^^^^^^^^^^^^^^

* Set to the same value of `prefix` in Quantum-ESPRESSO. It's used to locate the files `{prefix}.dyn*` or `{prefix}.phoebe.*.dat` generated by `ph.x`.
  
* *string*
  
* Required


.. _fermiLevel:

fermiLevel
^^^^^^^^^^

* Sets the fermi level of the ground state. Can be specified e.g. in Bands or DOS apps to specify an otherwise unknown fermi level. This quantity is read from file for transport calculations: this input parameter overwrites that value, use with caution.
  
* *double+units*
  
* Optional


.. _numOccupiedStates:

numOccupiedStates
^^^^^^^^^^^^^^^^^

* Determines the number of occupied Kohn-Sham states at the ground state. The default value might be read from the :ref:`electronH0Name` (when this is the Quantum-ESPRESSO xml file) or the file with the el-ph interaction (so, the user may not need to specify it for transport calculations). This value controls where the Fermi level is set. The user alternatively can specify the :ref:`fermiLevel` (and :ref:`numOccupiedStates` will be computed from the Fermi level).
  
* *double*
  
* Optional.


.. _minChemicalPotential:

minChemicalPotential
^^^^^^^^^^^^^^^^^^^^

* To be used together with :ref:`maxChemicalPotential` and :ref:`deltaChemicalPotential`, sets the code to compute properties at all chemical Potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with :ref:`chemicalPotentials` to instead manually specify the chemical potentials of the calculation.
  
* *double*
  
* (Required) either specify (:ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, :ref:`deltaChemicalPotential`) or :ref:`chemicalPotentials`.


.. _maxChemicalPotential:

maxChemicalPotential
^^^^^^^^^^^^^^^^^^^^

* To be used together with :ref:`minChemicalPotential` and :ref:`deltaChemicalPotential`, sets the code to compute properties at all chemical potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with chemicalPotentials to instead manually specify the chemical potentials of the calculation.
  
* *double*
  
* (Required) either specify (:ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, :ref:`deltaChemicalPotential`) or :ref:`chemicalPotentials`.


.. _deltaChemicalPotential:

deltaChemicalPotential
^^^^^^^^^^^^^^^^^^^^^^

* To be used together with :ref:`minChemicalPotential` and :ref:`maxChemicalPotential`, sets the code to compute properties at all chemical Potentials between :ref:`minChemicalPotential` and :ref:`maxChemicalPotential` in steps of :ref:`deltaChemicalPotential`. Can be exchanged with :ref:`chemicalPotentials` to instead manually specify the chemical potentials of the calculation.
  
* *double*
  
* (Required) either specify (:ref:`minChemicalPotential`, :ref:`maxChemicalPotential`, :ref:`deltaChemicalPotential`) or :ref:`chemicalPotentials`.


.. _minTemperature:

minTemperature
^^^^^^^^^^^^^^

* To be used together with :ref:`maxTemperature` and :ref:`deltaTemperature`, sets the code to compute observables at temperatures between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.
  
* *double*
  
* (Required): either set (:ref:`minTemperature`, :ref:`maxTemperature`, :ref:`deltaTemperature`) or :ref:`temperatures`.


.. _maxTemperature:

maxTemperature
^^^^^^^^^^^^^^

* To be used together with :ref:`minTemperature` and :ref:`deltaTemperature`, sets the code to compute observables at temperatures between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.
  
* *double*
  
* (Required): either set (:ref:`minTemperature`, :ref:`maxTemperature`, :ref:`deltaTemperature`) or :ref:`temperatures`.


.. _deltaTemperature:

deltaTemperature
^^^^^^^^^^^^^^^^

* To be used together with minTemperature and maxTemperature, sets the code to compute observables at temperatures between :ref:`minTemperature` and :ref:`maxTemperature` in steps of :ref:`deltaTemperature`.
  
* *double*

* (Required): either set (:ref:`minTemperature`, :ref:`maxTemperature`, :ref:`deltaTemperature`) or :ref:`temperatures`.

  
.. _useSymmetries:

useSymmetries
^^^^^^^^^^^^^

* When set to true, triggers the BTE to be computed only on the irreducible wedge of the Brillouin zone. For systems with several symmetries, this speeds up calculations. On the other hand, it may slow down the code for systems with few symmetries. If set to true, the viscosity is only computed at the RTA level.
  
* *bool*

* Default = `false`

* Optional

