@page Input Input files overview

Here we provide a description of the input variables of the various apps in the code.
Note that the text in the user input file is case sensitive!
In the code, we try to check the variables provided in input, and stop the code if we detect missing or invalid input parameters.
Note that the executable is the same for all this applications; the applications are selected through the @ref appName parameter.




@section phtr Phonon BTE

Target: build and solve the phonon Boltzmann Transport Equation (BTE), and compute phonon transport properties such as thermal conductivity, relaxation times and viscosity.

Input variables:
<ul>
<li> @ref appName = "phononTransport"</li>
<li> @ref phD2FileName</li>
<li> @ref phD3FileName</li>
<li> @ref sumRuleD2</li>
<li> @ref qMesh</li>
<li> @ref temperatures</li>
<li> @ref minTemperature</li>
<li> @ref maxTemperature</li>
<li> @ref deltaTemperature</li>
<li> @ref smearingMethod</li>
<li> @ref smearingWidth</li>
<li> @ref solverBTE</li>
<li> @ref scatteringMatrixInMemory</li>
<li> @ref windowType</li>
<li> @ref windowEnergyLimit</li>
<li> @ref windowPopulationLimit</li>
<li> @ref maxIterationsBTE</li>
<li> @ref convergenceThresholdBTE</li>
<li> @ref dimensionality</li>
<li> @ref constantRelaxationTime</li>
<li> @ref withIsotopeScattering</li>
<li> @ref massVariance</li>
<li> @ref boundaryLength</li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~



@section elwTransport Electron BTE, Wannier interpolation

Target: build and solve the electronic Boltzmann Transport Equation (BTE), using Wannier-based interpolation. Output quantites are electrical conductivity, electronic thermal conductivity, Seebeck coefficient, electron viscosity and electronic lifetimes.  

Input variables:
<ul>
<li> @ref appName = "electronWannierTransport" </li>
<li> @ref phD2FileName </li>
<li> @ref sumRuleD2 </li>
<li> @ref electronH0Name </li>
<li> @ref epwFileName </li>
<li> @ref kMesh </li>
<li> @ref temperatures </li>
<li> @ref minTemperature</li>
<li> @ref maxTemperature</li>
<li> @ref deltaTemperature</li>
<li> @ref dopings </li>
<li> @ref chemicalPotentials </li>
<li> @ref minChemicalPotential </li>
<li> @ref maxChemicalPotential </li>
<li> @ref deltaChemicalPotential </li>
<li> @ref smearingMethod </li>
<li> @ref smearingWidth </li>
<li> @ref windowType </li>
<li> @ref dimensionality </li>
<li> @ref constantRelaxationTime </li>
<li> @ref convergenceThresholdBTE </li>
<li> @ref maxIterationsBTE</li>
<li> @ref windowType </li>
<li> @ref windowEnergyLimit </li>
<li> @ref windowPopulationLimit </li>
<li> @ref solverBTE </li>
<li> @ref scatteringMatrixInMemory </li>
<li> @ref fermiLevel </li>
<li> @ref numOccupiedStates </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~



@section qe2Phoebe QE to Phoebe

Target: convert the electron-phonon coupling from a Quantum-ESPRESSO format to a Phoebe format.
In doing so, we postprocess data for the Wannier or EPA interpolations.

Input variables:
<ul>
<li> @ref appName = "elPhQeToPhoebe" </li>
<li> @ref phD2FileName </li>
<li> @ref elPhInterpolation </li>
<li> @ref electronH0Name </li>
<li> @ref wannier90Prefix </li>
<li> @ref quantumEspressoPrefix </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
appName = "elPhQeToPhoebe"
elPhInterpolation = "wannier"
phD2FileName = "./silicon.fc"
electronH0Name = "./si_tb.dat",
wannier90Prefix = "si"
quantumEspressoPrefix = "silicon"
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section phTau Phonon Lifetimes

Target: compute phonon lifetimes/linewidths on a path in the Brillouin zone.

Input variables:
<ul>
<li> @ref appName = "phononLifetimes" </li>
<li> @ref phD2FileName </li>
<li> @ref sumRuleD2 </li>
<li> @ref phD3FileName </li>
<li> @ref qMesh </li>
<li> @ref temperatures </li>
<li> @ref minTemperature</li>
<li> @ref maxTemperature</li>
<li> @ref deltaTemperature</li>
<li> @ref smearingMethod </li>
<li> @ref smearingWidth </li>
<li> @ref constantRelaxationTime</li>
<li> @ref withIsotopeScattering</li>
<li> @ref massVariance</li>
<li> @ref boundaryLength</li>
<li> @ref deltaPath </li>
<li> @ref beginEndPointPath </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section elwTau Electron Wannier Lifetimes

Target: compute electron lifetimes/linewidths on a path in the Brillouin zone, using Wannier-based interpolations.

Input variables:
<ul>
<li> @ref appName = "electronLifetimes" </li>
<li> @ref phD2FileName </li>
<li> @ref electronH0Name </li>
<li> @ref sumRuleD2 </li>
<li> @ref epwFileName </li>
<li> @ref kMesh </li>
<li> @ref temperatures </li>
<li> @ref minTemperature</li>
<li> @ref maxTemperature</li>
<li> @ref deltaTemperature</li>
<li> @ref dopings </li>
<li> @ref chemicalPotentials </li>
<li> @ref minChemicalPotential </li>
<li> @ref maxChemicalPotential </li>
<li> @ref deltaChemicalPotential </li>
<li> @ref smearingMethod </li>
<li> @ref smearingWidth </li>
<li> @ref constantRelaxationTime </li>
<li> @ref numOccupiedStates </li>
<li> @ref fermiLevel </li>
<li> @ref deltaPath </li>
<li> @ref beginEndPointPath </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section phDos Phonon Dos

Target: compute the phonon Density of States.

Input variables:
<ul>
<li> @ref appName = "phononDos" </li>
<li> @ref phD2FileName </li>
<li> @ref sumRuleD2 </li>
<li> @ref qMesh </li>
<li> @ref dosMinEnergy </li>
<li> @ref dosMaxEnergy </li>
<li> @ref dosDeltaEnergy </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
phD2FileName = "qespresso/silicon.fc",
sumRuleD2 = "simple"
qMesh = [10,10,10]
appName = "phononDos"
dosMinEnergy = 0. cmm1
dosMaxEnergy = 600. cmm1
dosDeltaEnergy = 0.5 cmm1
~~~~~~~~~~~~~~~~~~~~~~~~~~~



@section elwDos Electron DoS, Wannier interpolation

Target: compute the electronic Density of States. Electronic bands are interpolated to a finer mesh using maximally localized Wannier function interpolation.

Input variables:
<ul>
<li> @ref appName = "electronWannierDos" </li>
<li> @ref electronH0Name </li>
<li> @ref fermiLevel </li>
<li> @ref kMesh </li>
<li> @ref dosMinEnergy </li>
<li> @ref dosMaxEnergy </li>
<li> @ref dosDeltaEnergy </li>
<li> @ref beginEndCrystal </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section elfDos Electron DoS, Fourier interpolation

Target: compute the electronic Density of States. Electronic bands are interpolated to finer meshes using a Fourier interpolation.

Input variables:
<ul>
<li> @ref appName = "electronFourierDos" </li>
<li> @ref electronH0Name </li>
<li> @ref kMesh </li>
<li> @ref fermiLevel </li>
<li> @ref dosMinEnergy </li>
<li> @ref dosMaxEnergy </li>
<li> @ref dosDeltaEnergy </li>
<li> @ref electronFourierCutoff </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
electronH0Name = "qespresso/out/silicon.xml",
kMesh = [10,10,10]
appName = "electronFourierDos"
dosMinEnergy = -6. eV
dosMaxEnergy = 20. eV
dosDeltaEnergy = 0.1 eV
electronFourierCutoff = 4.
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section phBands Phonon Bands

Target: compute the phonon band structure on a path in the Brillouin zone.

Input variables:
<ul>
<li> @ref appName = "phononBands" </li>
<li> @ref phD2FileName </li>
<li> @ref sumRuleD2 </li>
<li> @ref deltaPath </li>
<li> @ref beginEndPointPath </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
phD2FileName = "qespresso/silicon.fc",
sumRuleD2 = "simple"
appName = "phononBands"
begin point path
L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000 
K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
end point path
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section elwBands Electron Bands, Wannier interpolation

Target: compute the phonon band structure on a path in the Brillouin zone. Electronic bands are interpolated using Wannier functions.

Input variables:
<ul>
<li> @ref appName = "electronWannierBands" </li>
<li> @ref electronH0Name </li>
<li> @ref fermiLevel </li>
<li> @ref deltaPath </li>
<li> @ref beginEndPointPath </li>
<li> @ref beginEndCrystal </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~


@section elfBands Electron Bands, Fourier interpolation

Target: compute the electronic band structure on a path in the Brillouin zone. Electronic bands are interpolated using Wannier functions.

Input variables:
<ul>
<li> @ref appName = "electronFourierBands" </li>
<li> @ref electronH0Name </li>
<li> @ref fermiLevel </li>
<li> @ref deltaPath </li>
<li> @ref electronFourierCutoff </li>
<li> @ref beginEndPointPath </li>
</ul>

Sample input file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~







@section VarDesc Variable descriptions

@subsubsection appName appName
<ul>
<li> This parameter, which must always be present, identifies which App (functionality) you want to run. Allowed values are:
  <ul>
  <li> "phononTransport": app to solve the @ref phtr and compute phonon transport properties.
  <li> "phononLifetimes": app to compute the @ref phTau.

  <li> "elPhQeToPhoebe": app to convert electron-phonon coupling from @ref qe2Phoebe (must be run before running any electron Transport).
  <li> "electronWannierTransport": app to solve the @ref elwTransport.
  <li> "electronLifetimes": app to compute the @ref elwTau.
  
  <li> "phononDos": app to compute the @ref phDos.
  <li> "electronWannierDos": app to compute the @ref elwDows.
  <li> "electronFourierDos": app to compute the @ref elfDos.
  
  <li> "phononBands": app to compute the @ref phBands on a path.
  <li> "electronWannierBands": app to compute the @ref elwBands on a path in the Brillouin zone.
  <li> "electronFourierBands": app to compute the @ref elfBands on a path in the Brillouin zone. 
</ul>
<li> *string*
<li> Required
</ul>


@subsubsection phD2FileName phD2FileName
<ul>
<li>Path to the file with harmonic force constants. File format supported is Quantum-ESPRESSO output of q2r.x.
<li> *string*
<li> Required
</ul>


@subsubsection phD3FileName phD3FileName
<ul>
<li> Path to the file with anharmonic 3rd order force constants.
<li> *string*
<li> Required
</ul>


@subsubsection sumRuleD2 sumRuleD2
<ul>
<li> If specified, applies an acoustic sum rule to the phonon harmonic force constants. Allowed values are "simple" or "crystal", with the same algorithm and meaning of Quantum-ESPRESSO matdyn.x program.
<li> *string*
<li> Required
</ul>

@subsubsection qMesh qMesh
<ul>
<li> Triplet of integers with the fine q-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of phonon properties.
<li> *list of int*
<li> Required
</ul>

@subsubsection kMesh kMesh
<ul>
<li> Triplet of integers with the fine k-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations of electronic properties. In electron-phonon transport calculations, `qMesh` is set to be equal to this value and does not need to be specified by the user.
<li> *list of int*
<li> Required
</ul>

@subsubsection temperatures temperatures
<ul>
<li> List with the values of temperatures to be used in the calculation. If scatteringMatrixInMemory=true, only one value of temperature is allowed.
<li> *list of doubles*
<li> Required
</ul>


@subsubsection smearingMethod smearingMethod
<ul>
<li> Selects the level of approximation for replacing the Dirac-delta approximating energy conservation. Allowed values are "gaussian" and "adaptiveGaussian" (preferred)
<li> *string*
<li> Required
</ul>


@subsubsection smearingWidth smearingWidth
<ul>
<li> This parameter is required if @ref smearingMethod="gaussian", where this parameter represents the full width at half maximum of the gaussian used to approximate the Dirac-delta conserving energy. Example: smearingWidth = 0.5 eV
<li> *double+units*
<li> Required if @ref smearingMethod="gaussian"
</ul>


@subsubsection solverBTE solverBTE
<ul>
<li> If specified, solves the Boltzmann equation beyond the relaxation time approximation. Allowed values are: "variational", "iterative", and "relaxons", see the Theory section for a detailed explanation. Example: solverBTE=["variational","relaxons"]
<li> *list of strings*
<li> Optional
</ul>


@subsubsection scatteringMatrixInMemory scatteringMatrixInMemory
<ul>
<li> If true, the scattering matrix is kept in memory, and only one temperature is allowed. In exchange for a larger memory usage, exact BTE solvers are much faster; disable this flag to reduce the memory footprint but slowing down the exact BTE solvers.
<li> *bool*
<li> Default=`true`
<li> Optional
</ul>


@subsubsection windowType windowType
<ul>
<li> Enables the window used to discard some phonon states that don't contribute to transport. Possible values are "nothing", "population" and "energy". "nothing" means window is not applied; "population" means phonon states are discarded if \f$ \frac{\partial \bar{n}}{\partial T} <\f$ windowPopulationLimit, where \f$ \frac{\partial \bar{n}}{\partial T}\f$ is the Bose--Einstein distribution derivative.
<li> *string*
<li> Default: `"nothing"`
<li> Optional
</ul>


@subsubsection windowEnergyLimit windowEnergyLimit
<ul>
<li> Additional parameter for energy @ref windowType. Specify two values \f$E_{min}\f$ and \f$E_{max}\f$(in electronVolts) such that we discard all phonon states  with energy outside of these bounds.
<li> *list of doubles*
<li> Required if @ref windowType="energy"
</ul>


@subsubsection windowPopulationLimit windowPopulationLimit
<ul>
<li> Required if @ref windowType="population". Cutoff values for discarding phonon states based on their equilibrium phonon occupation number, such that \f$ \frac{\partial \bar{n}}{\partial T} <\f$ windowPopulationLimit.
<li> *double*
<li> Required if @ref windowType="population"
</ul>


@subsubsection maxIterationsBTE maxIterationsBTE
<ul>
<li> Maximum number of iterations for iterative and variational BTE solvers. If the maximum number of iterations is reached, the code will throw an error.
<li> *int*
<li> Default: `50`
<li> Optional
</ul>


@subsubsection convergenceThresholdBTE convergenceThresholdBTE
<ul>
<li> Convergence criterion to stop iterative BTE solvers. The calculation is converged if the transport coefficients have a relative change smaller than convergenceThresholdBTE. 
<li> *double*
<li> Default: `1.0e-5`
<li> Optional
</ul>
 

@subsubsection dimensionality dimensionality
<ul>
<li> Input the dimensionality of the material. As a result, transport coefficients tensors will be of size (dim x dim), and units will be suitably scaled for the desired dimensionality.
<li> *int*
<li> Default: `3`
<li> Optional
</ul>


@subsubsection constantRelaxationTime constantRelaxationTime
<ul>
<li> If specified, we solve the BTE with the constant relaxation time approximation, where the phonon lifetime is set to this input value. (Fast but inaccurate!)
<li> *double+units*
<li> Optional
</ul>


@subsubsection withIsotopeScattering withIsotopeScattering
<ul>
<li> Controls whether to include or not phonon-isotope scattering
<li> *bool*
<li> Default: `true`
<li> Optional
</ul>


@subsubsection massVariance massVariance
<ul>
<li> User can specify a list of custom atomic mass variances \f$ g_2^s \f$. See Theory section for a description. The mass variances must be ordered in the same way that atomic species are specified in the file @ref phD2FileName. Defaults to the mass variance for natural isotopic abundance.
<li> *list of doubles*
<li> Default: natural isotopic abundance
<li> Optional
</ul>


@subsubsection boundaryLength boundaryLength
<ul>
<li> If specified, includes the phonon-boundary scattering within the RTA approximation. Example: boundaryLength = 10 mum
<li> *double+units*
<li> Optional
</ul>

@subsubsection electronH0Name electronH0Name
<ul>
<li> For Wannier-interpolation-based calculations, `electronH0Name` must contain the path to the `{prefix}_tb.dat` file generated by Wannier90. For Fourier-interpolation-based calculations, `electronH0Name` must contain the path to the Quantum-ESPRESSO {outdir}/{prefix}.xml file generated by `pw.x`.
<li> *string*
<li> Required
</ul>

@subsubsection dosMinEnergy
<ul>
<li> Used in conjunction with @ref dosMaxEnergy and @ref dosDeltaEnergy to compute the Density of States every @ref dosDeltaEnergy increments between @ref dosMinEnergy and @ref dosMaxEnergy.
<li> *double+uints*
<li> Required
</ul>

@subsubsection dosMaxEnergy dosMaxEnergy
<ul>
<li> Used in conjunction with @ref dosMinEnergy and @ref dosDeltaEnergy to compute the Density of States every @ref dosDeltaEnergy increments between @ref dosMinEnergy and @ref dosMaxEnergy.
<li> *double+uints*
<li> Required
</ul>

@subsubsection dosDeltaEnergy dosDeltaEnergy
<ul>
<li> Used in conjunction with @ref dosMinEnergy and @ref dosMaxEnergy to compute the Density of States every @ref dosDeltaEnergy increments between @ref dosMinEnergy and @ref dosMaxEnergy.
<li> *double+uints*
<li> Required
</ul>

@subsubsection electronFourierCutoff electronFourierCutoff
<ul>
<li> A parameter controlling the search of lattice vectors used for the Fourier interpolation of the electronic band structure. In detail, the lattice vectors used for the Fourier interpolation are searched in a supercell of size `electronFourierCutoff`\f$^3\f$ the primitive unit cell. Set it to at least 2.
<li> *double*
<li> Required
</ul>

@subsubsection beginEndCrystal begin/end crystal 
<ul>
<li> Specify the atomic species and atomic positions inside the crystal. This needs to be specified in some apps like WannierBands or WannierDos, as the output files of Wannier90 doesn't provide all the information about the crystal.
<li> Namelist format, with atomic symbol and position coordinates in units of Angstroms. Example:
~~~~~~~~{.c}
begin crystal
Si 0.00000   0.00000   0.00000
Si 1.34940   1.34940   1.34940
end crystal
~~~~~~~~
<li> Required
</ul>

@subsubsection beginEndPointPath begin/end point path
<ul>
<li> Specify the path of wavectors in the Brillouin zone used in apps such as `phononBands` or `phononLifetimes`. Use the parameter @ref deltaPath to control the number of wavevectors in each segment.
<li> Namelist format, as pairs of special point symbol and wavevector coordinates. Wavevector coordinates are in fractional coordinates with respect to the primitive reciprocal lattice vectors. Example:
~~~~~~~~{.c}
begin point path
L 0.50000  0.50000 0.5000 G 0.00000  0.00000 0.0000
G 0.00000  0.00000 0.0000 X 0.50000  0.00000 0.5000
X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000 
K 0.37500 -0.37500 0.0000 G 0.00000  0.00000 0.0000
end point path
~~~~~~~~
<li> Required
</ul>

@subsubsection dopings dopings
<ul>
<li> Specify a list of doping concentrations, in cm\f$^{-3}\f$, to compute electronic properties at various doping concentrations. The chemical potentials corresponding to this doping concentrations will be computed.
<li> *list of doubles*
<li> Required; alternatively one must specify @ref chemicalPotentials.
</ul> 

@subsubsection chemicalPotentials chemicalPotentials
<ul>
<li> Specify a list of chemical potentials to be used for the calculation of properties as a function of the chemical potential. If used in electron Wannier transport and scatteringMatrixInMemory=true, then only one value of chemical potentials can be specified. Values are in eV.
<li> *list of doubles*
<li> Required. The user can substitute this parameter by specifying `(minChemicalPotential,maxChemicalPotential,deltaChemicalPotential)`.
</ul>

@subsubsection epwFileName epwFileName
<ul>
<li> Path to the file generated by the app `elPhQeToPhoebe` containing the electron-phonon coupling in the Wannier representation (e.g. `{prefix}.phoebe.elph.dat`)
<li> *string*
<li> Required
</ul>

@subsubsection deltaPath deltaPath
<ul>
<li> This variable controls how far apart are the wavevectors when a path in the Brillouin zone is specified, and it represents the distance (in Bohr) between wavevectors. Can be used when a path of wavevectors is specified with the @ref beginEndPointPath key. 
<li> Default: 0.05 Bohr\f$^{-1}\f$
<li> *string*
<li> Optional
</ul>

@subsubsection elPhInterpolation elPhInterpolation
<ul>
<li> Can be either "wannier" or "epa". The first, prepares the electron-phonon coupling for the transport calculation with Wannier interpolation (i.e. does the transformation from Bloch to Wannier representation). The second, prepares the electron-phonon coupling to be used with the EPA approximation.
<li> *string*
<li> Required
</ul>

@subsubsection wannier90Prefix wannier90Prefix
<ul>
<li> Set to the same value of `prefix` in Wannier90. It's used to locate the files `{prefix}.eig` generated by `wannier90.x`.
<li> *string*
<li> Required
</ul>

@subsubsection quantumEspressoPrefix quantumEspressoPrefix
<ul>
<li> Set to the same value of `prefix` in Quantum-ESPRESSO. It's used to locate the files `{prefix}.dyn*` or `{prefix}.phoebe.*.dat` generated by `ph.x`.
<li> *string*
<li> Required
</ul>

@subsubsection fermiLevel fermiLevel
<ul>
<li> Sets the fermi level of the ground state. Can be specified e.g. in Bands or DOS apps to specify an otherwise unknown fermi level. This quantity is read from file for transport calculations: this input parameter overwrites that value, use with caution.
<li> *double+units*
<li> Optional
</ul>

@subsubsection numOccupiedStates numOccupiedStates
<ul>
<li> Determines the number of occupied Kohn-Sham states at the ground state. The default value can be read from either the @ref electronH0FileName (when this is the Quantum-ESPRESSO xml file) or the file with the el-ph interaction. The user choose instead to specify the @ref fermiLevel (@ref numOccupiedStates can be computed knowing that).
<li> *double*
<li> Optional.
</ul>

@subsubsection minChemicalPotential minChemicalPotential
<ul>
<li> To be used together with @ref maxChemicalPotential and @ref deltaChemicalPotential, sets the code to compute properties at all chemical Potentials between @ref minChemicalPotential and @ref maxChemicalPotential in steps of @ref deltaChemicalPotential. Can be exchanged with @ref chemicalPotentials to instead manually specify the chemical potentials of the calculation.
<li> *double*
<li> (Required) either specify (@ref minChemicalPotential, @ref maxChemicalPotential, @ref deltaChemicalPotential) or @ref chemicalPotentials.
</ul>

@subsubsection maxChemicalPotential maxChemicalPotential
<ul>
<li> To be used together with @ref minChemicalPotential and @ref deltaChemicalPotential, sets the code to compute properties at all chemical potentials between @ref minChemicalPotential and @ref maxChemicalPotential in steps of @ref deltaChemicalPotential. Can be exchanged with chemicalPotentials to instead manually specify the chemical potentials of the calculation.
<li> *double*
<li> (Required) either specify (@ref minChemicalPotential, @ref maxChemicalPotential, @ref deltaChemicalPotential) or @ref chemicalPotentials.
</ul>

@subsubsection deltaChemicalPotential deltaChemicalPotential
<ul>
<li> To be used together with @ref minChemicalPotential and @ref maxChemicalPotential, sets the code to compute properties at all chemical Potentials between @ref minChemicalPotential and @ref maxChemicalPotential in steps of @ref deltaChemicalPotential. Can be exchanged with @ref chemicalPotentials to instead manually specify the chemical potentials of the calculation.
<li> *double*
<li> (Required) either specify (@ref minChemicalPotential, @ref maxChemicalPotential, @ref deltaChemicalPotential) or @ref chemicalPotentials.
</ul>

@subsubsection minTemperature minTemperature
<ul>
<li> To be used together with @ref maxTemperature and @ref deltaTemperature, sets the code to compute observables at temperatures between @ref minTemperature and @ref maxTemperature in steps of @ref deltaTemperature.
<li> *double*
<li> (Required): either set (@ref minTemperature, @ref maxTemperature, @ref deltaTemperature) or @ref temperatures.
</ul>

@subsubsection maxTemperature maxTemperature
<ul>
<li> To be used together with @ref minTemperature and @ref deltaTemperature, sets the code to compute observables at temperatures between @ref minTemperature and @ref maxTemperature in steps of @ref deltaTemperature.
<li> *double*
<li> (Required): either set (@ref minTemperature, @ref maxTemperature, @ref deltaTemperature) or @ref temperatures.
</ul>

@subsubsection deltaTemperature deltaTemperature
<ul>
<li> To be used together with minTemperature and maxTemperature, sets the code to compute observables at temperatures between @ref minTemperature and @ref maxTemperature in steps of @ref deltaTemperature.
<li> *double*
<li> (Required): either set (@ref minTemperature, @ref maxTemperature, @ref deltaTemperature) or @ref temperatures.
</ul>
