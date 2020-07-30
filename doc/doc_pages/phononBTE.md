@page Theory Theory
@section PHBTE Phonon Boltzmann transport equation

We follow the theory of Ref: www.dx.doi.org/10.1103/PhysRevB.88.045430
Here we report the most relevant points of the manuscript. 

When a gradient of temperature \f$\nabla T\f$ is established in a system, a subsequent heat flux will start propagating in the medium.
Without loss of generality we assume the gradient of temperature to be along the direction \f$x\f$.
 The flux of heat, collinear to the temperature gradient, can be written in terms of phonon energies \f$\hbar\omega_{\mathbf{q}j}\f$, phonon group velocities \f$v_{q j}\f$ in the \f$x\f$ direction, and the perturbed phonon population \f$n_{q j}\f$:
 \begin{equation}
\frac{1}{N_0 \Omega} \sum_{q j} \hbar \omega_{q j} v_{q j} n_{q j} = - k \frac{\partial T}{ \partial x}
\end{equation}
On the l.h.s \f$\omega_{\mathbf{q}j }\f$ is the angular frequency of the phonon mode with vector \f$\mathbf{q}\f$ and branch index \f$j\f$, \f$\Omega\f$ is the volume of 
the unit cell and the sum runs over a uniform mesh of \f$N_0\f$ \f$\mathbf{q}\f$ points. 
On the r.h.s. \f$k\f$ is the diagonal component of the thermal conductivity in the temperature-gradient direction. % with \f$\alpha\f$ and $\beta$ the Cartesian indeces.
 Knowledge of the perturbed phonon population allows heat flux and subsequently thermal conductivity to be evaluated.
Unlike phonon scattering by defects, impurities and boundaries, anharmonic scattering represents an intrinsic resistive 
process and in high quality samples, at room temperature, it dominates the behaviour of lattice thermal conductivity balancing the perturbation due to the gradient of temperature.
The balance equation, namely the Boltzmann Transport Equation (BTE), formulated in 1929 by Peierls is:
\begin{equation}
-v_{q j}\frac {\partial T} {\partial x} \frac{\partial n_{q j}}{\partial T} + \frac{\partial n_{q j}}{\partial t}\bigg|_{scatt} = 0
\end{equation}
with the first term indicating the phonon diffusion due to the temperature gradient and the second term the scattering rate due to all the scattering processes.
This equation has to be solved self consistently.
In the general approach, for small perturbation from the equilibrium, the temperature gradient of the perturbed phonon population is replaced with the temperature gradient of the equilibrium phonon population \f$\partial n_{q j} / \partial T = \partial \bar{n}_{q j} / \partial T \f$ where \f$\bar{n}_{q j} = (e^{\hbar \omega_{q j} /k_BT} - 1)^{-1}\f$; while for the scattering term it can be expanded about its equilibrium value in terms of a first order perturbation \f$f^{\mathrm{EX}}\f$:
 \begin{equation}
 n_{q j} \simeq \bar{n}_{q j}+\bar{n}_{q j}(\bar{n}_{q j}+1) \frac{\partial T}{\partial x}f^{\mathrm{EX}}_{q j}
 \end{equation}
The linearized BTE can then be written in the following form:
\begin{equation}
-v_{q j}\left(\frac{\partial \bar{n}_{q j}}{\partial T}\right) =  
  \sum_{q' j',q'' j''}\Big[ P_{q j,q' j'}^{q'' j''}(f^{\mathrm{EX}}_{q j}+f^{\mathrm{EX}}_{q' j'}-f^{\mathrm{EX}}_{q'' j''}) + \frac{1}{2} P^{q' j',q'' j''}_{q j} (f^{\mathrm{EX}}_{q j}-f^{\mathrm{EX}}_{q' j'}-f^{\mathrm{EX}}_{q'' j''} )\Big] + \sum_{q' j'}  P^{\mathrm{isot}}_{q j,q' j'}  (f^{\mathrm{EX}}_{q j} - f^{\mathrm{EX}}_{q' j'}) + P^{\mathrm{be}}_{q j} f^{\mathrm{EX}}_{q j}
\end{equation}
where the sum on \f$q'\f$ and \f$q''\f$ is performed in the Brillouin Zone (BZ).
The \f$\mathrm{EX}\f$ superscript of the first order perturbation \f$f^{\mathrm{EX}}\f$ denotes the exact solution of the BTE, to be distinguished from the approximated solutions that we will discuss later.
In this last equation the anharmonic scattering processes as well as the scattering with the isotopic impurities and the border effect are considered. 
More specifically \f$P_{q j,q' j'}^{q'' j''}\f$ is the scattering rate at the equilibrium  of a process where a phonon mode $q j$ scatters by absorbing another mode \f$qp j'\f$ to generate a third phonon mode \f$q'' j''\f$.
While \f$P^{q' j',q'' j''}_{q j}\f$ is the scattering rate at the equilibrium of a process where a phonon mode \f$\mathbf{q}j\f$ decays in two modes \f$\mathbf{q}'j'\f$ and \f$\mathbf{q}''j'' \f$. 

The two scattering rates have the forms:
\begin{eqnarray}
P^{qpp j''}_{q j,qp j'}&{=}& \frac{2 \pi}{N_0 \hbar^2} \sum_{\mathbf{G}}
	  |V^{(3)}(q j,q' j',{-}q'' j'')|^2  \nonumber \\
	&&    \bar{n}_{q j}\bar{n}_{q' j'}(\bar{n}_{q'' j''}+1) \delta_{q{+}q' {-}q'', \mathbf{G}}\nonumber \\
	&&  \delta(\hbar \omega_{q j} +\hbar \omega_{q' j'}-\hbar \omega_{q'' j''}) \label{coal}  
\end{eqnarray}
\begin{eqnarray}
P^{q' j',q'' j''}_{q j}&{=}& \frac{2 \pi}{N_0 \hbar^2 } \sum_{\mathbf{G}}
	    |V^{(3)}(q j,{-}q' j',{-}q'' j'')|^2 \nonumber \\
        &&       \bar{n}_{q j}(\bar{n}_{q' j'}{+}1)(\bar{n}_{q'' j''}{+}1)\delta_{q{-}q' {-}q'', \mathbf{G}} \nonumber \\ 
        &&   \delta(\hbar \omega_{q j}-\hbar \omega_{q' j'}-\hbar \omega_{q'' j''} )\label{dec} 
\end{eqnarray}
 with \f$ {\mathbf G}\f$ the reciprocal lattice vectors.
In order to evaluate them it is necessary to compute the third derivative \f$V^{(3)}\f$ of  the total energy of the crystal \f$\mathcal{E}^{tot}(\{u_{s \alpha} (\mathbf{R}_l) \})\f$, with respect to the atomic displacement \f$u_{s \alpha} (\mathbf{R}_l)\f$, from the equilibrium position, of the s-th atom, 
along the \f$\alpha\f$ Cartesian coordinate in the crystal cell identified by the lattice vector \f$\mathbf{R}_l\f$ :
\begin{equation}
V^{(3)}(\mathbf{q} j,\mathbf{q}' j',\mathbf{q}'' j'')= \frac{\partial^3 \mathcal{E}^{cell}}
                                                    {\partial X_{\mathbf{q} j},\partial X_{\mathbf{q}' j'},\partial X_{\mathbf{q}'' j''}}
\end{equation}
 where \f$\mathcal{E}^{cell}\f$ is the energy per unit cell.
 The non-dimensional quantity \f$X_{\mathbf{q} j}\f$ is defined by 
\begin{equation}
X_{\mathbf{q} j}= \frac{1}{N_0}\sum_{l,s,\alpha} \sqrt{\frac{2 M_s \omega_{q j}} {\hbar}} z^{s \alpha^*}_{q j}  u_{s \alpha }(\mathbf{R}_l) e^{-i\mathbf{q}\cdot \mathbf{R}_l}
\end{equation}
with \f$z^{s \alpha}_{\mathbf{q}j} \f$ being the orthogonal phonon eigenmodes normalized on the unit cell and $M_s$ the atomic masses.
This expression of X can be used to transform the output of a density-functional code, i.e. the matrix of energy derivatives in real space \f$\mathcal{E}(R_l s\alpha,R'_{l'} s' \alpha',R''_{l''}s''\alpha'')\f$ to the Fourier space.
The matrix is actually a periodic function, so it can be possible to neglect one of the Bravais lattice vector indices of such a tensor.
Note that, Quantum Espresso provides the matrix \f$\mathcal{E}(0 s\alpha,R'_{l'} s' \alpha',R''_{l''}s''\alpha'')\f$ while frozen phonon codes such as phonopy and related use \f$\mathcal{E}(R_l s\alpha,R'_{l'} s' \alpha',0 s''\alpha'')\f$, i.e. set a different bravais lattice vector to zero.
Phoebe uses the latter convention at the time being.

The rate of the elastic scattering with isotopic impurities has the form:
\begin{eqnarray}
  P_{q j,q' j'}^{\mathrm{isot}} & = & \frac{\pi}{2 N_0} \omega_{q j}\omega_{q' j'}  
                   \left[ \bar{n}_{q j} \bar{n}_{q' j'} + \frac{\bar{n}_{q j} + \bar{n}_{q' j'}} {2} \right ] \nonumber \\
                 & &\sum_{s} g^{s}_{2}   |  \sum_{\alpha} z^{s \alpha^*}_{\mathbf{q}j} \cdot z^{s \alpha}_{q' j'} |^2 \delta (\omega_{q j}- \omega_{q' j'})
\end{eqnarray}
where \f$g^s_2 = \frac{(M_s - \langle  M_s\rangle)^2}{ \langle M_s \rangle^2 }\f$ is the average over the mass distribution of the atom of type \f$s\f$.
In presence of two isotopes \f$M_s\f$ and \f$M_{s'}\f$ it can be written in terms of the concentration \f$\epsilon\f$ and mass change \f$\Delta M_s= M_{s'} - M_s\f$ :
\begin{equation}
 g^s_2=  \epsilon(1-\epsilon)  \frac{ | \Delta M_s |}{ \langle M_s \rangle} 
\end{equation}
with \f$\langle M_s \rangle = M_s + \epsilon \Delta M_s\f$.\\
Eventually, in a system of finite size, \f$P_{q j}^{\mathrm{be}} \f$ describes the reflection of a phonon from the border:
\begin{equation}
P_{q j}^{\mathrm{be}} = \frac{v_{q j}}{L}\bar{n}_{q j}(\bar{n}_{q j}+1) 
\end{equation}
where \f$L\f$ is the Casimir length of the sample.
The border scattering is treated in the relaxation time approximation and it results in a process in which a phonon from a specific state(\f$\mathbf{q} j\f$) is reemitted from the surface contributing only to the equilibrium distribution.

For the sake of clarity we will contract from here on the vector \f$\mathbf{q}\f$ and branch index \f$j\f$ in a single mode index \f$\nu\f$.
The BTE of Eq. \ref{BTE2} can be written as  a linear system in matrix form:
\begin{equation}
\mathbf{A} \mathbf{f}^{\mathrm{EX}}=\mathbf{b}
\label{linearsyst}
\end{equation}
with the vector \f$b_{\nu'} =-v_{\nu'}\hbar \omega_{\nu'} \bar{n}_{\nu'}(\bar{n}_{\nu'}+1) \f$ and the matrix
\begin{equation}
A_{\nu,\nu'} = [{\sum_{\nu'',\nu'''}} (P^{\nu''}_{\nu,\nu'''} + \frac{ P_{\nu''',\nu''}^{\nu}}{2} ) + \sum_{\nu''} P^{\mathrm{isot}}_{\nu,\nu''} + P^{\mathrm{be}}_{\nu} ] \delta_{\nu,\nu'} - {\sum_{\nu''}} (  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu}  ) + P^{\mathrm{isot}}_{\nu,\nu'} 
\end{equation}
where we have used \f$P^{\nu', \nu''}_{\nu}=P_{\nu', \nu''}^{\nu}\f$ from the detailed balance condition \f$\bar{n}_{\nu}(\bar{n}_{\nu'}+1)(\bar{n}_{\nu''}+1) = (\bar{n}_{\nu}+1)\bar{n}_{\nu'}\bar{n}_{\nu''}\f$ (valid under the assumption \f$\hbar \omega = \hbar \omega' + \hbar \omega''\f$).
In this form the matrix is symmetric and positive semi-definite and it can be decomposed in \f$\mathbf{A} = \mathbf{A}^{\mathrm{out}} +\mathbf{A}^{\mathrm{in}} \f$,
where
\begin{eqnarray}
A^{\mathrm{out}}_{\nu,\nu'} &=& \frac{\bar{n}_{\nu}(\bar{n}_{\nu} +1)} {\tau^{\mathrm{T}}_{\nu}}\delta_{\nu,\nu'} \\
A^{\mathrm{in}}_{\nu,\nu'} &=&  -  \sum_{\nu''} \left(  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu} \right )    + P^{\mathrm{isot}}_{\nu,\nu'} 
\end{eqnarray}
with \f$\tau^{\mathrm{T}}_{\nu}\f$ being the phonon relaxation time.
The \f$\mathbf{A}^{\mathrm{out}}\f$ diagonal matrix describes the depopulation of phonon states due to the scattering mechanisms while the \f$\mathbf{A}^{\mathrm{in}}\f$ matrix describes their repopulation due to the incoming scattered phonons.

The solution of the linear system in Eq. \ref{linearsyst} is obtained formally by inverting the matrix \f${\mathbf A}\f$.
\begin{equation}
{\mathbf f}^{\mathrm{EX}} =   \frac{1}{\mathbf{A}}  {\mathbf b}
\end{equation}
and subsequently the thermal conductivity will be evaluated as:
\begin{equation}
k =  \lambda {\mathbf b} \cdot {\mathbf f}^{\mathrm{EX}}
= - \frac{\hbar}{N_0\Omega  k_B T^2}\sum_{\nu}v_{\nu}
				    \omega_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1) f_{\nu}^{\mathrm{EX}}
\end{equation}
with \f$\lambda= 1 /(N_0\Omega k_B T^2)\f$.
 

@section PHRTA RTA solution of the phonon BTE
In the relaxation time approximation (RTA), we set \f$\mathbf{A}^{\mathrm{in}}\f$ to zero
\begin{equation}
{\mathbf f}^{\mathrm{SMA}} =\frac{1}{ \mathbf{A}^{\mathrm{out}}}  {\mathbf b}
\end{equation}
Inverting \f$\mathbf{A}^{\mathrm{out}}\f$ is trivial due to its diagonal form.
The lattice thermal conductivity in RTA is then 
\begin{equation}
k^{\mathrm{RTA}}=\lambda \mathbf{b} \cdot \mathbf{f}^{\mathrm{SMA}}=\frac{\hbar^2}{N_0\Omega k_B T^2}\sum_{\nu}v^2_{\nu} \omega^2_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1)\tau^{\mathrm{T}}_{\nu}.
\end{equation}


@section PHITER Iterative solution of the phonon BTE - Omini-Sparavigna method

An exact solution of the BTE that does not imply either storing or the explicit inversion of matrix \f$\mathbf{A}\f$ has been proposed by Omini and Sparavigna by converging with respect to the iteration \f$i\f$ the following:
\begin{equation}
\mathbf{f}_{ i+1} =\frac{1} {\mathbf{A}^{\mathrm{out} } } \mathbf{b} - \frac{1} {\mathbf{A}^{\mathrm{out} } } \mathbf{A}^{\mathrm{in}}  \mathbf{f}_{i}
\end{equation}
with the iteration zero consisting in the RTA \f$\mathbf{f}_0=\mathbf{f}^{\mathrm{RTA}}\f$.
Instead of storing and inverting \f$\mathbf{A}\f$, it just requires the evaluation of \f$\mathbf{A}^{\mathrm{in}}\:\mathbf{f}_{i}\f$, at each iteration \f$i\f$ of the OS method, which is an operation computationally much less demanding.
Once the convergence is obtained the thermal conductivity is evaluated by:
\begin{equation}
k^{\mathrm{NV}}(\mathbf{f}_i)=\lambda \mathbf{b}\cdot \mathbf{f}_{i}
\label{kOS}
\end{equation}
From a mathematical point of views the OS iterative procedure 
can be written as a geometric series:
 \begin{equation}
\mathbf{f}_{ i} = \sum_{j=0,i} \left(-\frac{1}{\mathbf{A}^{\mathrm{out}}}  \mathbf{A}^{\mathrm{in}}\right)^{j} \frac{1}{\mathbf{A}^{\mathrm{out}}} \:  \mathbf{b} \;.
\end{equation}

@section PHVAR Iterative solution of the phonon BTE - Variational method

An alternative approach consists in using the properties of the matrix \f${\mathbf A} \f$ to find the exact solution of the linearized BTE, via the variational principle.
Indeed the solution  of the BTE is the vector \f$\mathbf{f}^{\mathrm{EX}}\f$ which makes  stationary the quadratic form
\begin{equation}
\mathcal{F}(\mathbf{f}) =\frac{1}{2} {\mathbf f} \cdot{\mathbf A} {\mathbf f}- {\mathbf b} \cdot {\mathbf f}
\end{equation}
for a generic vector \f$\mathbf{f}\f$.
Since $\mathbf{A}$ is positive the stationary point is the global and single minimum of this functional.
One can then define a variational conductivity functional: 
\begin{equation} 
k^\mathrm{V}(\mathbf{f}) = - 2 \lambda \mathcal{F}({\mathbf f})
\label{quadratic}
\end{equation}
that has the property \f$k^\mathrm{V}(\mathbf{f}^{\mathrm{EX}})=k\f$ while any other value of \f$k^{\mathrm{V}}(\mathbf{f})\f$  underestimates \f$k\f$.
In other words, finding the minimum of the quadratic form is equivalent to maximizing the thermal conductivity functional. 
As a consequence an error \f$\delta \mathbf{f}= \mathbf{f} - \mathbf{f}^{\mathrm{EX}}\f$  results in an error in conductivity, linear in \f$\delta \mathbf{f}\f$ when using the non-variational estimator, and quadratic in the variational form.

Here we solve the BTE on a grid (as in OS procedure) by using the conjugate gradient method, to obtain the exact solution of the BTE equation.
In order to speed up the convergence of the conjugate gradient we take advantage of the diagonal and dominant role of \f$\mathbf{A}^{\mathrm{out}}\f$ and we use a preconditioned conjugate gradient.
Formally, this corresponds to use in the minimization the rescaled variable:
\begin{equation}
\tilde{{\mathbf f}} = \sqrt{{\mathbf A^{\mathrm{out}}}} {\mathbf f}
\label{prec1}
\end{equation}
and then, with respect to this new variable, minimize the quadratic form \f$\tilde{\mathcal{F}}(\tilde{\mathbf{f}}) = \mathcal{F}(\mathbf{f})\f$ where:
\begin{equation}
\tilde{\mathcal{F}}( \tilde{\mathbf{f}}) =\frac{1}{2} \tilde{\mathbf{f}}\cdot \tilde{\mathbf{A}} \tilde{\mathbf{f}}- \tilde{\mathbf{ b}}\cdot\tilde{\mathbf {f}}
\end{equation}
and  
\begin{equation}
\tilde{{\mathbf A}} =\frac{1}{ \sqrt{{\mathbf A^{\mathrm{out}}}}} {\mathbf A}\frac{1}{ \sqrt{{\mathbf A^{\mathrm{out}}}}}
\end{equation}
\begin{equation}
\tilde{{\mathbf b}} =\frac{1}{ \sqrt{{\mathbf A^{\mathrm{out}}}}} {\mathbf b} \label{prec3}
\end{equation}

Notice that \f$\tilde{\mathbf{f}}^{\mathrm{RTA}}=\tilde{\mathbf{b}}\f$.
The square root evaluation of \f$\mathbf{A}^{\mathrm{out}}\f$ is trivial due to its diagonal form.
The computational cost per iteration of the conjugate gradient scheme is equivalent to the OS one, but it always converges and requires a smaller number of iterations.


The conjugate gradient minimization requires the evaluation of the gradient \f$\mathbf{g}_i= \mathbf{A} \mathbf{f}_i - \mathbf{b}\f$ and a line minimization.
Since the form is quadratic the line minimization can be done analytically and exactly.
Moreover the information required by the line minimization at  iteration \f$i\f$ can be recycled to compute the gradient at the next iteration \f$i+1\f$.
Starting with an the initial vector \f$\mathbf{f}_0= \mathbf{f}^{\mathrm{RTA}}\f$, initial gradient \f$\mathbf{g}_0=\mathbf{A}\mathbf{f}_0 -\mathbf{f}^{\mathrm{RTA}}\f$ and letting \f$\mathbf{h}_0= -\mathbf{g}_0\f$, the conjugate gradient method can be summarized with the
recurrence:
\begin{equation}
\mathbf{t}_i =\mathbf{A} \mathbf{h}_i
\end{equation}
\begin{equation}
  {\mathbf f}_{i+1} = {\mathbf f}_{i} - \frac{\mathbf {g}_{i} \cdot {\mathbf{h}_{i}} } {\mathbf{h}_{i} \cdot \mathbf{t}_i } \mathbf{h}_{i}
  \end{equation}
  \begin{equation}
\mathbf{g}_{i+1} = \mathbf{g}_{i}-\frac{\mathbf {g}_{i} \cdot {\mathbf{h}_{i}} } {\mathbf{h}_{i} \cdot \mathbf{t}_i }\mathbf{t}_i
\end{equation}
\begin{equation}
 \mathbf{h}_{i+1} = -\mathbf{g}_{i+1} + \frac{\mathbf{g}_{i+1} \cdot \mathbf{g}_{i+1}}{{\mathbf{g}_{i}} \cdot {\mathbf{g}_{i}} }  {\mathbf h}_{i} 
\end{equation}
where \f$\mathbf{h}_i\f$ is the search direction and \f$\mathbf{t}_i\f$ is an auxiliary vector.
Notice that each iteration requires only one application of the matrix \f$\mathbf{A}\f$ on the vector \f$\mathbf{h}_i\f$ as in the OS method.




@section PHRELAXONS Relaxons solution to the BTE
TBD after deciding how to handle symmetries


@section SYMM BTE with symmetries
To be written



@section VELOCITY Phonon velocity operator
The velocity operator matrix elements (e.g. along the x direction) can be computed using the Hellmann-Feynman theorem from the Dynamical matrix \f$\mathbf{\mathcal{D}}\f$:
\begin{equation}
V_{j j'}(q) = \sum_{\alpha \alpha' s s'} \frac{1}{2 \sqrt{M_s M_{s'}} \omega_{q j} }  z^{s \alpha^*}_{q j}  \frac{\partial \mathcal{D}^{\alpha \alpha'}_{s s'}(q)}{ \partial q_x}   z^{s' \alpha'}_{q j'}
\end{equation}
In the non-degenerate case, the group velocity is \f$v_{q j}=V_{j j}(q)\f$ while in the degenerate one we use the phonon polarization vectors that diagonalize the matrix in the degenerate subspace.



@section SMEARING Dirac-delta approximations
The delta function for the energy conservation can be replaced by a Gaussian
\begin{equation}
 \delta(\hbar \omega)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar \omega/ \sigma )^2)} \;,
\end{equation}
where \f$\sigma\f$ is a constant decided by user input.
It is important to note that when the delta function is substituted with a Gaussian the detailed balance condition is only valid under approximation.
The definition used above guarantees that the scattering matrix is symmetric and non-negative. 

Another method is the adaptive-gaussian smearing scheme [https://link.aps.org/doi/10.1103/PhysRevB.75.195121].
Specifically, we want to approximate a dirac-delta function of the form:
\begin{equation}
 \delta(\hbar (\omega_1+\omega_2-\omega_3)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar (\omega_1+\omega_2-\omega_3)/ \sigma )^2)}.
\end{equation}
This time, we allow \f$\sigma\f$ to be a value dependent on the energies.
Specifically, we build it as:
\begin{equation}
\sigma = \frac{1}{\sqrt{12}} \sqrt{ \sum_{\beta} (\sum_{\alpha} (v_2-v_3) \frac{M_{\alpha \beta}}{N_{\beta}}  )^2 }
\end{equation}
where \f$M\f$ is a matrix comprised of the primitive cell lattice vectors (each column is a lattice vector) ,\f$v_2\f$ and \f$v_3\f$ are phonon group velocities, and \f$N_{\beta}\f$ is the number of wavevectors sampled along direction \f$\beta\f$.
Note that the adaptive scheme may be critical in the case where the velocity sum to zero: in that case, we skip the scattering event, unless we have an exact energy conservation taking place.





