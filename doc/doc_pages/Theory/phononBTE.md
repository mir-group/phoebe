@page thPhononBTE Phonon Boltzmann transport equation

We follow the theory of the reference available [at this link] (https://arxiv.org/abs/1212.0470).
Here we report the most relevant points of the manuscript. 

When a gradient of temperature \f$\nabla T\f$ is established in a system, a subsequent heat flux will start propagating in the medium.
Without loss of generality we assume the gradient of temperature to be along the direction \f$x\f$.
 The flux of heat, collinear to the temperature gradient, can be written in terms of phonon energies \f$\hbar\omega_{\boldsymbol{q}j}\f$, phonon group velocities \f$v_{q j}\f$ in the \f$x\f$ direction, and the perturbed phonon population \f$n_{q j}\f$:
 \begin{equation}
\frac{1}{N_0 \Omega} \sum_{q j} \hbar \omega_{q j} \boldsymbol{v}_{\boldsymbol{q} j} n_{\boldsymbol{q} j} = - k \frac{\partial T}{ \partial \boldsymbol{x}}
\end{equation}
On the l.h.s \f$\omega_{\boldsymbol{q}j }\f$ is the angular frequency of the phonon mode with vector \f$\boldsymbol{q}\f$ and branch index \f$j\f$, \f$\Omega\f$ is the volume of 
the unit cell and the sum runs over a uniform mesh of \f$N_0\f$ \f$\boldsymbol{q}\f$ points. 
On the r.h.s. \f$k\f$ is the diagonal component of the thermal conductivity in the temperature-gradient direction. % with \f$\alpha\f$ and $\beta$ the Cartesian indeces.
 Knowledge of the perturbed phonon population allows heat flux and subsequently thermal conductivity to be evaluated.
Unlike phonon scattering by defects, impurities and boundaries, anharmonic scattering represents an intrinsic resistive 
process and in high quality samples, at room temperature, it dominates the behaviour of lattice thermal conductivity balancing the perturbation due to the gradient of temperature.
The balance equation, namely the Boltzmann Transport Equation (BTE), formulated in 1929 by Peierls is:
\begin{equation}
-\boldsymbol{v}_{q j}\cdot \frac {\partial T} {\partial \boldsymbol{x}} \frac{\partial n_{\boldsymbol{q} j}}{\partial T} + \frac{\partial n_{\boldsymbol{q} j}}{\partial t}\bigg|_{scatt} = 0
\end{equation}
with the first term indicating the phonon diffusion due to the temperature gradient and the second term the scattering rate due to all the scattering processes.
This equation has to be solved self consistently.
In the general approach, for small perturbation from the equilibrium, the temperature gradient of the perturbed phonon population is replaced with the temperature gradient of the equilibrium phonon population \f$\partial n_{\boldsymbol{q} j} / \partial T = \partial \bar{n}_{\boldsymbol{q} j} / \partial T \f$ where \f$\bar{n}_{\boldsymbol{q} j} = (e^{\hbar \omega_{\boldsymbol{q} j} /k_BT} - 1)^{-1}\f$; while for the scattering term it can be expanded about its equilibrium value in terms of a first order perturbation \f$f^{\mathrm{EX}}\f$:
 \begin{equation}
 n_{\boldsymbol{q} j} \simeq \bar{n}_{\boldsymbol{q} j}+\bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q} j}+1) \frac{\partial T}{\partial \boldsymbol{x}}\cdot f^{\mathrm{EX}}_{\boldsymbol{q} j}
 \end{equation}
The linearized BTE can then be written in the following form:
\begin{equation}
-v_{\boldsymbol{q} j}\left(\frac{\partial \bar{n}_{\boldsymbol{q} j}}{\partial T}\right) =  
  \sum_{\boldsymbol{q}' j',\boldsymbol{q}'' j''}\Big[ P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\boldsymbol{q}'' j''}(f^{\mathrm{EX}}_{\boldsymbol{q} j}+f^{\mathrm{EX}}_{\boldsymbol{q}' j'}-f^{\mathrm{EX}}_{\boldsymbol{q}'' j''}) + \frac{1}{2} P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j} (f^{\mathrm{EX}}_{\boldsymbol{q} j}-f^{\mathrm{EX}}_{\boldsymbol{q}' j'}-f^{\mathrm{EX}}_{\boldsymbol{q}'' j''} )\Big] + \sum_{\boldsymbol{q}' j'}  P^{\mathrm{isot}}_{\boldsymbol{q} j,\boldsymbol{q}' j'}  (f^{\mathrm{EX}}_{\boldsymbol{q} j} - f^{\mathrm{EX}}_{\boldsymbol{q}' j'}) + P^{\mathrm{be}}_{\boldsymbol{q} j} f^{\mathrm{EX}}_{\boldsymbol{q} j}
\end{equation}
where the sum on \f$\boldsymbol{q}'\f$ and \f$\boldsymbol{q}''\f$ is performed in the Brillouin Zone (BZ).
The \f$\mathrm{EX}\f$ superscript of the first order perturbation \f$f^{\mathrm{EX}}\f$ denotes the exact solution of the BTE, to be distinguished from the approximated solutions that we will discuss later.
In this last equation the anharmonic scattering processes as well as the scattering with the isotopic impurities and the border effect are considered. 
More specifically \f$P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\boldsymbol{q}'' j''}\f$ is the scattering rate at the equilibrium  of a process where a phonon mode $\boldsymbol{q} j$ scatters by absorbing another mode \f$\boldsymbol{q}' j'\f$ to generate a third phonon mode \f$\boldsymbol{q}'' j''\f$.
While \f$P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j}\f$ is the scattering rate at the equilibrium of a process where a phonon mode \f$\boldsymbol{q}j\f$ decays in two modes \f$\boldsymbol{q}'j'\f$ and \f$\boldsymbol{q}''j'' \f$. 

The two scattering rates have the forms:
\begin{eqnarray}
P^{\boldsymbol{q}'' j''}_{\boldsymbol{q} j,\boldsymbol{q}' j'}&{=}& \frac{2 \pi}{N_0 \hbar^2} \sum_{\boldsymbol{G}}
	  |V^{(3)}(\boldsymbol{q} j,\boldsymbol{q}' j',-\boldsymbol{q}'' j'')|^2  \nonumber \\
	&&    \bar{n}_{\boldsymbol{q} j}\bar{n}_{\boldsymbol{q}' j'}(\bar{n}_{\boldsymbol{q}'' j''}+1) \delta_{\boldsymbol{q}+\boldsymbol{q}' -\boldsymbol{q}'', \boldsymbol{G}}\nonumber \\
	&&  \delta(\hbar \omega_{\boldsymbol{q} j} +\hbar \omega_{\boldsymbol{q}' j'}-\hbar \omega_{\boldsymbol{q}'' j''}) \label{coal}  
\end{eqnarray}
\begin{eqnarray}
P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j}&{=}& \frac{2 \pi}{N_0 \hbar^2 } \sum_{\boldsymbol{G}}
	    |V^{(3)}(\boldsymbol{q} j,-\boldsymbol{q}' j',-\boldsymbol{q}'' j'')|^2 \nonumber \\
        &&       \bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q}' j'}+1)(\bar{n}_{\boldsymbol{q}'' j''}+1)\delta_{\boldsymbol{q}-\boldsymbol{q}' -\boldsymbol{q}'', \boldsymbol{G}} \nonumber \\ 
        &&   \delta(\hbar \omega_{\boldsymbol{q} j}-\hbar \omega_{\boldsymbol{q}' j'}-\hbar \omega_{\boldsymbol{q}'' j''} )\label{dec} 
\end{eqnarray}
 with \f$\boldsymbol{G}\f$ the reciprocal lattice vectors.
In order to evaluate them it is necessary to compute the third derivative \f$V^{(3)}\f$ of  the total energy of the crystal \f$\mathcal{E}^{tot}(\{u_{s \alpha} (\boldsymbol{R}_l) \})\f$, with respect to the atomic displacement \f$u_{s \alpha} (\boldsymbol{R}_l)\f$, from the equilibrium position, of the s-th atom, 
along the \f$\alpha\f$ Cartesian coordinate in the crystal cell identified by the lattice vector \f$\boldsymbol{R}_l\f$ :
\begin{equation}
V^{(3)}(\boldsymbol{q} j,\boldsymbol{q}' j',\boldsymbol{q}'' j'')= \frac{\partial^3 \mathcal{E}^{cell}}
                                                    {\partial X_{\boldsymbol{q} j},\partial X_{\boldsymbol{q}' j'},\partial X_{\boldsymbol{q}'' j''}}
\end{equation}
 where \f$\mathcal{E}^{cell}\f$ is the energy per unit cell.
 The non-dimensional quantity \f$X_{\boldsymbol{q} j}\f$ is defined by 
\begin{equation}
X_{\boldsymbol{q} j}= \frac{1}{N_0}\sum_{l,s,\alpha} \sqrt{\frac{2 M_s \omega_{\boldsymbol{q} j}} {\hbar}} z^{s \alpha^*}_{\boldsymbol{q} j}  u_{s \alpha }(\boldsymbol{R}_l) e^{-i\boldsymbol{q}\cdot \boldsymbol{R}_l}
\end{equation}
with \f$z^{s \alpha}_{\boldsymbol{q}j} \f$ being the orthogonal phonon eigenmodes normalized on the unit cell and $M_s$ the atomic masses.
This expression of X can be used to transform the output of a density-functional code, i.e. the matrix of energy derivatives in real space \f$\mathcal{E}(\boldsymbol{R}_l s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{R}''_{l''}s''\alpha'')\f$ to the Fourier space.
The matrix is actually a periodic function, so it can be possible to neglect one of the Bravais lattice vector indices of such a tensor.
Note that, Quantum Espresso provides the matrix \f$\mathcal{E}(\boldsymbol{0} s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{R}''_{l''}s''\alpha'')\f$ while frozen phonon codes such as phonopy and related use \f$\mathcal{E}(\boldsymbol{R}_l s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{0} s''\alpha'')\f$, i.e. set a different bravais lattice vector to zero.
Phoebe uses the latter convention at the time being.

The rate of the elastic scattering with isotopic impurities has the form:
\begin{eqnarray}
  P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\mathrm{isot}} & = & \frac{\pi}{2 N_0} \omega_{\boldsymbol{q} j}\omega_{\boldsymbol{q}' j'}  
                   \left[ \bar{n}_{\boldsymbol{q} j} \bar{n}_{\boldsymbol{q}' j'} + \frac{\bar{n}_{\boldsymbol{q} j} + \bar{n}_{\boldsymbol{q}' j'}} {2} \right ] \nonumber \\
                 & &\sum_{s} g^{s}_{2}   |  \sum_{\alpha} z^{s \alpha^*}_{\boldsymbol{q}j} \cdot z^{s \alpha}_{\boldsymbol{q}' j'} |^2 \delta (\omega_{\boldsymbol{q} j}- \omega_{\boldsymbol{q}' j'})
\end{eqnarray}
where \f$g^s_2 = \frac{(M_s - \langle  M_s\rangle)^2}{ \langle M_s \rangle^2 }\f$ is the average over the mass distribution of the atom of type \f$s\f$.
In presence of two isotopes \f$M_s\f$ and \f$M_{s'}\f$ it can be written in terms of the concentration \f$\epsilon\f$ and mass change \f$\Delta M_s= M_{s'} - M_s\f$ :
\begin{equation}
 g^s_2=  \epsilon(1-\epsilon)  \frac{ | \Delta M_s |}{ \langle M_s \rangle} 
\end{equation}
with \f$\langle M_s \rangle = M_s + \epsilon \Delta M_s\f$.\\
Eventually, in a system of finite size, \f$P_{q j}^{\mathrm{be}} \f$ describes the reflection of a phonon from the border:
\begin{equation}
P_{\boldsymbol{q} j}^{\mathrm{be}} = \frac{v_{\boldsymbol{q} j}}{L}\bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q} j}+1) 
\end{equation}
where \f$L\f$ is the Casimir length of the sample.
The border scattering is treated in the relaxation time approximation and it results in a process in which a phonon from a specific state(\f$\boldsymbol{q} j\f$) is reemitted from the surface contributing only to the equilibrium distribution.

For the sake of clarity we will contract from here on the vector \f$\boldsymbol{q}\f$ and branch index \f$j\f$ in a single mode index \f$\nu\f$.
The BTE of Eq. \ref{BTE2} can be written as  a linear system in matrix form:
\begin{equation}
\boldsymbol{A} \boldsymbol{f}^{\mathrm{EX}}=\boldsymbol{b}
\label{linearsyst}
\end{equation}
with the vector \f$b_{\nu'} =-v_{\nu'}\hbar \omega_{\nu'} \bar{n}_{\nu'}(\bar{n}_{\nu'}+1) \f$ and the matrix
\begin{equation}
A_{\nu,\nu'} = [{\sum_{\nu'',\nu'''}} (P^{\nu''}_{\nu,\nu'''} + \frac{ P_{\nu''',\nu''}^{\nu}}{2} ) + \sum_{\nu''} P^{\mathrm{isot}}_{\nu,\nu''} + P^{\mathrm{be}}_{\nu} ] \delta_{\nu,\nu'} - {\sum_{\nu''}} (  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu}  ) + P^{\mathrm{isot}}_{\nu,\nu'} 
\end{equation}
where we have used \f$P^{\nu', \nu''}_{\nu}=P_{\nu', \nu''}^{\nu}\f$ from the detailed balance condition \f$\bar{n}_{\nu}(\bar{n}_{\nu'}+1)(\bar{n}_{\nu''}+1) = (\bar{n}_{\nu}+1)\bar{n}_{\nu'}\bar{n}_{\nu''}\f$ (valid under the assumption \f$\hbar \omega = \hbar \omega' + \hbar \omega''\f$).
In this form the matrix is symmetric and positive semi-definite and it can be decomposed in \f$\boldsymbol{A} = \boldsymbol{A}^{\mathrm{out}} +\boldsymbol{A}^{\mathrm{in}} \f$,
where
\begin{eqnarray}
A^{\mathrm{out}}_{\nu,\nu'} &=& \frac{\bar{n}_{\nu}(\bar{n}_{\nu} +1)} {\tau^{\mathrm{T}}_{\nu}}\delta_{\nu,\nu'} \\
A^{\mathrm{in}}_{\nu,\nu'} &=&  -  \sum_{\nu''} \left(  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu} \right )    + P^{\mathrm{isot}}_{\nu,\nu'} 
\end{eqnarray}
with \f$\tau^{\mathrm{T}}_{\nu}\f$ being the phonon relaxation time.
The \f$\boldsymbol{A}^{\mathrm{out}}\f$ diagonal matrix describes the depopulation of phonon states due to the scattering mechanisms while the \f$\boldsymbol{A}^{\mathrm{in}}\f$ matrix describes their repopulation due to the incoming scattered phonons.

The solution of the linear system in Eq. \ref{linearsyst} is obtained formally by inverting the matrix \f${\boldsymbol A}\f$.
\begin{equation}
{\boldsymbol f}^{\mathrm{EX}} =   \frac{1}{\boldsymbol{A}}  {\boldsymbol b}
\end{equation}
and subsequently the thermal conductivity will be evaluated as:
\begin{equation}
k =  \lambda {\boldsymbol b} \cdot {\boldsymbol f}^{\mathrm{EX}}
= - \frac{\hbar}{N_0\Omega  k_B T^2}\sum_{\nu}v_{\nu}
				    \omega_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1) f_{\nu}^{\mathrm{EX}}
\end{equation}
with \f$\lambda= 1 /(N_0\Omega k_B T^2)\f$.
 

@section PHRTA RTA solution of the phonon BTE
In the relaxation time approximation (RTA), we set \f$\boldsymbol{A}^{\mathrm{in}}\f$ to zero
\begin{equation}
{\boldsymbol f}^{\mathrm{SMA}} =\frac{1}{ \boldsymbol{A}^{\mathrm{out}}}  {\boldsymbol b}
\end{equation}
Inverting \f$\boldsymbol{A}^{\mathrm{out}}\f$ is trivial due to its diagonal form.
The lattice thermal conductivity in RTA is then 
\begin{equation}
k^{\mathrm{RTA}}=\lambda \boldsymbol{b} \cdot \boldsymbol{f}^{\mathrm{SMA}}=\frac{\hbar^2}{N_0\Omega k_B T^2}\sum_{\nu}v^2_{\nu} \omega^2_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1)\tau^{\mathrm{T}}_{\nu}.
\end{equation}


@section thPHITER Iterative solution of the phonon BTE - Omini-Sparavigna method

Note: generally, we recommend the variational method over this. 

An exact solution of the BTE that does not imply either storing or the explicit inversion of matrix \f$\boldsymbol{A}\f$ has been proposed by Omini and Sparavigna by converging with respect to the iteration \f$i\f$ the following:
\begin{equation}
\boldsymbol{f}_{ i+1} =\frac{1} {\boldsymbol{A}^{\mathrm{out} } } \boldsymbol{b} - \frac{1} {\boldsymbol{A}^{\mathrm{out} } } \boldsymbol{A}^{\mathrm{in}}  \boldsymbol{f}_{i}
\end{equation}
with the iteration zero consisting in the RTA \f$\boldsymbol{f}_0=\boldsymbol{f}^{\mathrm{RTA}}\f$.
Instead of storing and inverting \f$\boldsymbol{A}\f$, it just requires the evaluation of \f$\boldsymbol{A}^{\mathrm{in}}\:\boldsymbol{f}_{i}\f$, at each iteration \f$i\f$ of the OS method, which is an operation computationally much less demanding.
Once the convergence is obtained the thermal conductivity is evaluated by:
\begin{equation}
k^{\mathrm{NV}}(\boldsymbol{f}_i)=\lambda \boldsymbol{b}\cdot \boldsymbol{f}_{i}
\label{kOS}
\end{equation}
From a mathematical point of views the OS iterative procedure 
can be written as a geometric series:
 \begin{equation}
\boldsymbol{f}_{ i} = \sum_{j=0,i} \left(-\frac{1}{\boldsymbol{A}^{\mathrm{out}}}  \boldsymbol{A}^{\mathrm{in}}\right)^{j} \frac{1}{\boldsymbol{A}^{\mathrm{out}}} \:  \boldsymbol{b} \;.
\end{equation}

@section PHVAR Iterative solution of the phonon BTE - Variational method

An alternative approach consists in using the properties of the matrix \f${\boldsymbol A} \f$ to find the exact solution of the linearized BTE, via the variational principle.
Indeed the solution  of the BTE is the vector \f$\boldsymbol{f}^{\mathrm{EX}}\f$ which makes  stationary the quadratic form
\begin{equation}
\mathcal{F}(\boldsymbol{f}) =\frac{1}{2} {\boldsymbol f} \cdot{\boldsymbol A} {\boldsymbol f}- {\boldsymbol b} \cdot {\boldsymbol f}
\end{equation}
for a generic vector \f$\boldsymbol{f}\f$.
Since $\boldsymbol{A}$ is positive the stationary point is the global and single minimum of this functional.
One can then define a variational conductivity functional: 
\begin{equation} 
k^\mathrm{V}(\boldsymbol{f}) = - 2 \lambda \mathcal{F}({\boldsymbol f})
\label{quadratic}
\end{equation}
that has the property \f$k^\mathrm{V}(\boldsymbol{f}^{\mathrm{EX}})=k\f$ while any other value of \f$k^{\mathrm{V}}(\boldsymbol{f})\f$  underestimates \f$k\f$.
In other words, finding the minimum of the quadratic form is equivalent to maximizing the thermal conductivity functional. 
As a consequence an error \f$\delta \boldsymbol{f}= \boldsymbol{f} - \boldsymbol{f}^{\mathrm{EX}}\f$  results in an error in conductivity, linear in \f$\delta \boldsymbol{f}\f$ when using the non-variational estimator, and quadratic in the variational form.

Here we solve the BTE on a grid (as in OS procedure) by using the conjugate gradient method, to obtain the exact solution of the BTE equation.
In order to speed up the convergence of the conjugate gradient we take advantage of the diagonal and dominant role of \f$\boldsymbol{A}^{\mathrm{out}}\f$ and we use a preconditioned conjugate gradient.
Formally, this corresponds to use in the minimization the rescaled variable:
\begin{equation}
\tilde{{\boldsymbol f}} = \sqrt{{\boldsymbol A^{\mathrm{out}}}} {\boldsymbol f}
\label{prec1}
\end{equation}
and then, with respect to this new variable, minimize the quadratic form \f$\tilde{\mathcal{F}}(\tilde{\boldsymbol{f}}) = \mathcal{F}(\boldsymbol{f})\f$ where:
\begin{equation}
\tilde{\mathcal{F}}( \tilde{\boldsymbol{f}}) =\frac{1}{2} \tilde{\boldsymbol{f}}\cdot \tilde{\boldsymbol{A}} \tilde{\boldsymbol{f}}- \tilde{\boldsymbol{ b}}\cdot\tilde{\boldsymbol {f}}
\end{equation}
and  
\begin{equation}
\tilde{{\boldsymbol A}} =\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}} {\boldsymbol A}\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}}
\end{equation}
\begin{equation}
\tilde{{\boldsymbol b}} =\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}} {\boldsymbol b} \label{prec3}
\end{equation}

Notice that \f$\tilde{\boldsymbol{f}}^{\mathrm{RTA}}=\tilde{\boldsymbol{b}}\f$.
The square root evaluation of \f$\boldsymbol{A}^{\mathrm{out}}\f$ is trivial due to its diagonal form.
The computational cost per iteration of the conjugate gradient scheme is equivalent to the OS one, but it always converges and requires a smaller number of iterations.


The conjugate gradient minimization requires the evaluation of the gradient \f$\boldsymbol{g}_i= \boldsymbol{A} \boldsymbol{f}_i - \boldsymbol{b}\f$ and a line minimization.
Since the form is quadratic the line minimization can be done analytically and exactly.
Moreover the information required by the line minimization at  iteration \f$i\f$ can be recycled to compute the gradient at the next iteration \f$i+1\f$.
Starting with an the initial vector \f$\boldsymbol{f}_0= \boldsymbol{f}^{\mathrm{RTA}}\f$, initial gradient \f$\boldsymbol{g}_0=\boldsymbol{A}\boldsymbol{f}_0 -\boldsymbol{f}^{\mathrm{RTA}}\f$ and letting \f$\boldsymbol{h}_0= -\boldsymbol{g}_0\f$, the conjugate gradient method can be summarized with the
recurrence:
\begin{equation}
\boldsymbol{t}_i =\boldsymbol{A} \boldsymbol{h}_i
\end{equation}
\begin{equation}
  {\boldsymbol f}_{i+1} = {\boldsymbol f}_{i} - \frac{\boldsymbol {g}_{i} \cdot {\boldsymbol{h}_{i}} } {\boldsymbol{h}_{i} \cdot \boldsymbol{t}_i } \boldsymbol{h}_{i}
  \end{equation}
  \begin{equation}
\boldsymbol{g}_{i+1} = \boldsymbol{g}_{i}-\frac{\boldsymbol {g}_{i} \cdot {\boldsymbol{h}_{i}} } {\boldsymbol{h}_{i} \cdot \boldsymbol{t}_i }\boldsymbol{t}_i
\end{equation}
\begin{equation}
 \boldsymbol{h}_{i+1} = -\boldsymbol{g}_{i+1} + \frac{\boldsymbol{g}_{i+1} \cdot \boldsymbol{g}_{i+1}}{{\boldsymbol{g}_{i}} \cdot {\boldsymbol{g}_{i}} }  {\boldsymbol h}_{i} 
\end{equation}
where \f$\boldsymbol{h}_i\f$ is the search direction and \f$\boldsymbol{t}_i\f$ is an auxiliary vector.
Notice that each iteration requires only one application of the matrix \f$\boldsymbol{A}\f$ on the vector \f$\boldsymbol{h}_i\f$ as in the OS method.




@section PHRELAXONS Relaxons solution to the BTE
We first diagonalize the scattering matrix:
\begin{equation}
\frac{1}{N_k} \sum_{\nu} \Omega_{\nu\nu'} \theta_{\nu'\alpha} = \frac{1}{\tau_{\alpha}} \theta_{\nu\alpha}
\end{equation}
where \f$ \theta \f$ are eigenvectors, \f$ \alpha \f$ are eigenvalue indices, \f$ \frac{1}{\tau_{\alpha}} \f$ are eigenvalues, and the scattering matrix is:
\begin{equation}
\Omega_{\nu\nu'} = \frac{ A_{\nu\nu'} } { \sqrt{ \bar{n}_{\nu}(\bar{n}_{\nu}+1) \bar{n}_{\nu'}(\bar{n}_{\nu'}+1)  } }
\end{equation}
Next, we compute the velocities:
\begin{equation}
\boldsymbol{V}_{\alpha} = \frac{1}{N_k} \sum_{\nu} \theta_{\nu0} \boldsymbol{v}_{\nu} \theta_{\nu\alpha}
\end{equation}
where
\begin{equation}
\theta_{\nu0} = \sqrt { \frac{ \frac{\partial \bar{n}_{\nu}}{\partial \epsilon} }{C T} } \hbar \omega_{\nu}
\end{equation}
Finally, the thermal conductivity is:
\begin{equation}
k^{ij} = \sum_{\alpha} C V_{\alpha}^i V_{\alpha}^j \tau_{\alpha}
\end{equation}






@section VELOCITY Phonon velocity operator
The velocity operator matrix elements (e.g. along the x direction) can be computed using the Hellmann-Feynman theorem from the Dynamical matrix \f$\boldsymbol{\mathcal{D}}\f$:
\begin{equation}
V^x_{j j'}(\boldsymbol{q}) = \sum_{\alpha \alpha' s s'} \frac{1}{2 \sqrt{M_s M_{s'}} \omega_{\boldsymbol{q} j} }  z^{s \alpha^*}_{\boldsymbol{q} j}  \frac{\partial \mathcal{D}^{\alpha \alpha'}_{s s'}(\boldsymbol{q})}{ \partial q_x}   z^{s' \alpha'}_{\boldsymbol{q} j'}
\end{equation}
In the non-degenerate case, the group velocity is \f$\boldsymbol{v}_{\boldsymbol{q} j}=\boldsymbol{V}_{j j}(\boldsymbol{q})\f$ while in the degenerate one we use the phonon polarization vectors that diagonalize the matrix in the degenerate subspace.



@section thSMEARING Dirac-delta approximations
The delta function for the energy conservation can be replaced by a Gaussian
\begin{equation}
 \delta(\hbar \omega)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar \omega/ \sigma )^2)} \;,
\end{equation}
where \f$\sigma\f$ is a constant decided by user input.
It is important to note that when the delta function is substituted with a Gaussian the detailed balance condition is only valid under approximation.
The definition used above guarantees that the scattering matrix is symmetric and non-negative. 

Another method is the [adaptive-gaussian smearing scheme](https://link.aps.org/doi/10.1103/PhysRevB.75.195121).
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





@section WIGNERPH Wigner correction to phonon thermal conductivity.

The theory is fully described in the Reference available at this [link](https://www.nature.com/articles/s41567-019-0520-x).

In extreme synthesis, the thermal conductivity is estimated as:
\begin{equation}
k_{\alpha\beta} = k^{BTE}_{\alpha\beta} +  \frac{k_BT^2}{\Omega N_k} \sum_{\boldsymbol{q}} \sum_{s\neq s'} \frac{\omega_{\boldsymbol{q}j}+\omega_{\boldsymbol{q}j'}}{2}   V_{jj'}^{\alpha}(\boldsymbol{q}) V_{j'j}^{\beta}(\boldsymbol{q}) \frac{ ( \frac{\partial n_{\boldsymbol{q}j}}{\partial T} + \frac{\partial n_{\boldsymbol{q}j'}}{\partial T})(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'}) }{4(\omega_{\boldsymbol{q}j}-\omega_{\boldsymbol{q}j'})^2 + (\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})^2} 
\end{equation}

where \f$k^{BTE}_{\alpha\beta}\f$ is the thermal conductivity estimated by the Boltzmann transport equation discussed above, and \f$\Gamma_{\boldsymbol{q}j} = \frac{1}{\tau_{\boldsymbol{q}j}}\f$ is the phonon linewidth, i.e. a diagonal element of the scattering matrix.


@section thPhViscosity Thermal Viscosity
The theory is described to far greater extent in this <a href="https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011019">PRX</a>.
The equilibrium of a system of bosonic particles that conserves energy and momentum is the drifting distribution:
\begin{equation}
n_{\nu}^{D}
=
\frac{1}{e^{\beta(\hbar \omega_\nu - \hbar \boldsymbol{q} \cdot \boldsymbol{u})}-1} \;,
\end{equation}
where \f$\boldsymbol{q} \f$ is the phonon wavevector (proportional to the phonon crystal momentum, and \f$\boldsymbol{u}\f$ is the phonon drift velocity.
The thermal viscosity is defined as the coefficient of proportionality between the crystal momentum flux \f$\Pi\f$ and a local perturbation in the drift velocity \f$\boldsymbol{u}\f$.
\begin{equation}
\Pi^{ij} = - \sum_{kl} \eta^{ijkl} \frac{\partial u^k}{\partial r^l}
\end{equation}
and the momentum flux (at least, the component relevant to our case) is defined as:
\begin{equation}
\Pi^{ij} = \frac{1}{V N_q} \sum_{\nu} \hbar q^i v_{\nu}^j n_{\nu}
\end{equation}
The population in response to the perturbation is fixed by the phonon BTE.
At the RTA level, we simply need to solve
\begin{equation}
\boldsymbol{v}_{\nu} \cdot (\frac{\partial n^{D}_{\nu}}{\partial \boldsymbol{u}} \cdot \nabla \boldsymbol{u} )
= - \frac{n_{\nu}}{\tau_{\nu}}
\end{equation}
We linearize the solution, stating \f$ n_{\nu} = n_{\nu} \nabla \boldsymbol{u} \f$, and the equation is readily solved.
Beyond the RTA, we must solve the equation:
\begin{equation}
\frac{\boldsymbol{v}_{\nu}}{\sqrt{\bar{n}_{\nu}(\bar{n}_{\nu}+1)}} \cdot (\frac{\partial n^{D}_{\nu}}{\partial \boldsymbol{u}} \cdot \nabla \boldsymbol{u} )
= - \frac{1}{V N_q} \sum_{\nu'} \Omega_{\nu\nu'} n_{\nu'}
\end{equation}
which we do with the eigenvector formalism.
Using the eigenvectors of the scattering matrix, we expand the phonon population as:
\begin{equation}
n_{\nu} = \sum_{kl} f^{kl}_{\alpha} \theta_{\nu\alpha} \nabla_l u^k
\end{equation}
We find the solution as:
\begin{equation}
f^{kl}_{\alpha} = - \tau_{\alpha} \sum_{\nu} \theta_{\nu\alpha} \frac{\boldsymbol{v}^l_{\nu}}{\sqrt{\bar{n}_{\nu}(\bar{n}_{\nu}+1)}} \frac{\partial n^{D}_{\nu}}{\partial u^k} 
\end{equation}
which can be used to reconstruct the phonon population response.
Finally, the viscosity tensor is symmetrized, finding the thermal viscosity:
\begin{equation}
\mu^{ijkl} = \frac{1}{2} ( \eta^{ijkl} + \eta^{ilkj} )
\end{equation}
The code also prints other quantities needed to write the viscous heat equations derived <a href="https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011019">in this reference</a>.


@section thPhBTESymms Symmetries of the BTE
We exploit the symmetries of the crystal to speed-up the calculation of thermal conductivity.
We took inspiration from <a href='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.265506'>this reference</a>.
Let \f$q\f$ indicate any wavevector in the Brillouin zone.
The symmetries of a crystal identify an irreducible set of wavevectors \f$ q^* \f$, such that any other wavevector \f$q\f$ can be obtained from a rotation of these irreducible wavevectors \f$ q = R q^* \f$.
The basic idea is to restrict the calculation to the irreducible set of wavevectors.
The conductivity for example, is:
\begin{equation}
k^{ij}
= \frac{1}{V N_k} \sum_{\nu} \hbar \omega_{\nu} v^i_{\nu} n^{j}_{\nu}
= \frac{1}{V N_k} \sum_{\nu^*} \sum_{R} \hbar \omega_{\nu^*} (R v_{\nu})_{i} (R n_{\nu})_{j}
\end{equation}
where \f$ R \f$ is the set of rotations used to reconstruct all the symmetry-equivalent wavevectors of \f$q^*\f$, and the summation over \f$\nu^*\f$ is only done in the irreducible set of wavevectors.

The BTE too can be restricted to the irreducible wedge.
\begin{equation}
v^i_{\nu^*} \frac{\partial \bar{n}_{\nu}}{\partial T}
= - \frac{1}{V N_q} \sum_{\nu'} A_{\nu^*\nu'} f^i_{\nu'}
= - \frac{1}{V N_q} \sum_{\nu'^*} \sum_{R} \sum_{j} A_{\nu^*\nu'^*} R_{ij} f^j_{\nu'}
= - \frac{1}{V N_q} \sum_{\nu'^* j} A^{ij}_{\nu^*\nu'^*} f^j_{\nu'}
\end{equation}
Hence, one can work with the same techniques detailed above, provided that we work with an enlarged matrix \f$ A^{ij}_{\nu^*\nu'^*} \f$.

By default, we make use of symmetries.
Some comments:
<ul>
<li> Advantage: for a system with a lot of symmetries, the matrix \f$ A^{ij}_{\nu^*\nu'^*} \f$ is generally smaller than \f$ A_{\nu\nu'} \f$, and thus calculations will be much faster.
<li> Disadvantage 1: we cannot compute viscosity beyond the RTA using symmetries. To do so, one must disable symmetries.
<li> Disadvantage 2: note that the symmetric matrix gains two indices on cartesian coordinates. As a result, in the limit case that there are no symmetries in the system (only the identity), the matrix \f$ A^{ij}_{\nu^*\nu'^*} \f$ will still be computed on the same number of wavevectors of  \f$ A_{\nu\nu'} \f$, but occupies 3x3 times more memory without adding any information. Therefore, for low-symmetry systems, consider disabling symmetries.
</ul>
