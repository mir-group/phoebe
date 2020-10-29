@page Input Input files overview

@section phtr Phonon transport

<table>
<caption id="phTransport">Input variables for phonon transport</caption>
<tr><th> Name <th> Kind <th> Default <th> Required <th> Description

<tr><td> <b>appName</b> <td> *string* <td> `None` <td> Required <td> Set to "phononTransport" to run this application.

<tr><td> <b>phD2FileName</b> <td> *string*      <td> `None` <td> Required <td> Path to the file with harmonic force constants. File format supported is Quantum-ESPRESSO output of q2r.x.

<tr><td> <b>phD3FileName</b> <td> *string*      <td> `None` <td> Required <td> Path to the file with anharmonic 3rd order force constants.

<tr><td> <b>sumRuleD2</b> <td> *string*      <td> `None` <td> Required <td> If specified, applies an acoustic sum rule to the phonon harmonic force constants. Allowed values are "simple" or "crystal", with the same algorithm and meaning of Quantum-ESPRESSO matdyn.x program.

<tr><td> <b>qMesh</b> <td> *int list*    <td> `None` <td> Required <td> Triplet of integers with the fine q-point Monkhorst-Pack mesh that will be used for Brillouin zone integrations.

<tr><td> <b>temperatures</b> <td> *double list* <td> `None` <td> Required <td> List with the values of temperatures to be used in the calculation. If scatteringMatrixInMemory=true, only one value of temperature is allowed.

<tr><td> <b>smearingMethod</b> <td> *string*      <td> `None` <td> Required <td> Selects the level of approximation for replacing the Dirac-delta approximating energy conservation. Allowed values are "gaussian" and "adaptiveGaussian" (preferred)

<tr><td> <b>smearingWidth</b>  <td> *double+unit* <td> `None` <td> (Required) <td> This parameter is required if smearingMethod="gaussian", where this parameter represents the full width at half maximum of the gaussian used to approximate the Dirac-delta conserving energy. Example: smearingWidth = 0.5 eV

<tr><td> <b>solverBTE</b> <td> *string list* <td> `None` <td> Optional <td> If specified, solves the Boltzmann equation beyond the relaxation time approximation. Allowed values are: "variational", "iterative", and "relaxons", see the Theory section for a detailed explanation. Example: solverBTE=["variational","relaxons"]

<tr><td> <b>scatteringMatrixInMemory</b> <td> *bool* <td> `true` <td> Optional <td> If true, the scattering matrix is kept in memory, and only one temperature is allowed. In exchange for a larger memory usage, exact BTE solvers are much faster; disable this flag to reduce the memory footprint but slowing down the exact BTE solvers.

<tr><td> <b>windowType</b> <td> *string* <td> `"nothing"` <td> Optional <td> Enables the window used to discard some phonon states that don't contribute to transport. Possible values are "nothing", "population" and "energy". "nothing" means window is not applied; "population" means phonon states are discarded if dn/dT<windowPopulationLimit, where dn/dT is the Bose--Einstein distribution derivative.

<tr><td> <b>windowEnergyLimit</b> <td> *double list* <td> `None` <td> (Required) <td> Required if windowType="energy". Specify two values \f$E_{min}\f$ and \f$E_{max}\f$(in electronVolts) such that we discard all phonon states  with energy outside of these bounds.

<tr><td> <b>windowPopulationLimit</b> <td> *double* <td> `None` <td> (Required) <td> Required if windowType="population". Cutoff values for discarding phonon states based on their equilibrium phonon occupation number, such that dn/dT<windowPopulationLimit. 

<tr><td> <b>maxIterationsBTE</b> <td> *int* <td> `50` <td> Optional <td> Maximum number of iterations for iterative and variational BTE solvers. If the maximum number of iterations is reached, the code will throw an error.
 
<tr><td> <b>convergenceThresholdBTE</b> <td> *double* <td> `1.0e-5` <td> Optional <td> Convergence criterion to stop iterative BTE solvers. The calculation is converged if the transport coefficients have a relative change smaller than convergenceThresholdBTE. 
 
<tr><td> <b>dimensionality</b> <td> *int* <td> `3` <td> Optional <td> Input the dimensionality of the material. As a result, transport coefficients tensors will be of size dimxdim, and units will be suitably scaled for the desired dimensionality.

<tr><td> <b>constantRelaxationTime</b> <td> *double+unit* <td> `None` <td> Optional <td> If specified, we solve the BTE with the constant relaxation time approximation, where the phonon lifetime is set to this input value. (Fast but inaccurate!)

<tr><td> <b>withIsotopeScattering</b> <td> *bool* <td> `true` <td> Optional <td> Controls whether to include or not phonon-isotope scattering

<tr><td> <b>massVariance</b> <td> *double list* <td> Natural abundance <td> Optional <td> User can specify a list of custom atomic mass variances \f$ g_2^s \f$. See Theory section for a description. The mass variances must be ordered in the same way that atomic species are specified in the file "phD2FileName" . Defaults to the mass variance for natural isotopic abundance.
  
<tr><td> <b>boundaryLength</b> <td> *double+unit* <td> `None` <td> Optional <td> If specified, includes the phonon-boundary scattering within the RTA approximation. Example: boundaryLength = 10 mum

</tr>
</table>

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



@section eltr Electron transport with Wannier interpolation
@section epa Electron transport - EPA approximation
@section phdos Phonon density of states
@section elwdos Electron density of states with Wannier interpolation
@section elfdos Electron density of states with Fourier interpolation
@section phbs Phonon band structure
@section elwbs Electron band structure with Wannier interpolation
@section elfbs Electron band structure with Fourier interpolation

