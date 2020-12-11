Phonon BTE
==========

We follow the theory of the reference available at this link (https://arxiv.org/abs/1212.0470). 
Here we report the most relevant points of the manuscript. 

When a gradient of temperature :math:`\nabla T` is established in a system, a subsequent heat flux will start propagating in the medium.
Without loss of generality we assume the gradient of temperature to be along the direction :math:`x`.
The flux of heat, collinear to the temperature gradient, can be written in terms of phonon energies :math:`\hbar\omega_{\boldsymbol{q}j}`, phonon group velocities :math:`v_{q j}` in the :math:`x` direction, and the perturbed phonon population :math:`n_{q j}`
 
.. math::
   \frac{1}{N_0 \Omega} \sum_{q j} \hbar \omega_{q j} \boldsymbol{v}_{\boldsymbol{q} j} n_{\boldsymbol{q} j} = - k \frac{\partial T}{ \partial \boldsymbol{x}}

On the l.h.s :math:`\omega_{\boldsymbol{q}j }` is the angular frequency of the phonon mode with vector :math:`\boldsymbol{q}` and branch index :math:`j`, :math:`\Omega` is the volume of the unit cell and the sum runs over a uniform mesh of :math:`N_0 \boldsymbol{q}` points. 
On the r.h.s. :math:`k` is the diagonal component of the thermal conductivity in the temperature-gradient direction, with :math:`\alpha` and :math:`\beta` the Cartesian indeces.
Knowledge of the perturbed phonon population allows heat flux and subsequently thermal conductivity to be evaluated.
Unlike phonon scattering by defects, impurities and boundaries, anharmonic scattering represents an intrinsic resistive process and in high quality samples, at room temperature, it dominates the behaviour of lattice thermal conductivity balancing the perturbation due to the gradient of temperature.
The balance equation, namely the Boltzmann Transport Equation (BTE), formulated in 1929 by Peierls is:

.. math::
   -\boldsymbol{v}_{q j}\cdot \frac {\partial T} {\partial \boldsymbol{x}} \frac{\partial n_{\boldsymbol{q} j}}{\partial T} + \frac{\partial n_{\boldsymbol{q} j}}{\partial t}\bigg|_{scatt} = 0

with the first term indicating the phonon diffusion due to the temperature gradient and the second term the scattering rate due to all the scattering processes.
This equation has to be solved self consistently.
In the general approach, for small perturbation from the equilibrium, the temperature gradient of the perturbed phonon population is replaced with the temperature gradient of the equilibrium phonon population :math:`\partial n_{\boldsymbol{q} j} / \partial T = \partial \bar{n}_{\boldsymbol{q} j} / \partial T ` where `\bar{n}_{\boldsymbol{q} j} = (e^{\hbar \omega_{\boldsymbol{q} j} /k_BT} - 1)^{-1}`; while for the scattering term it can be expanded about its equilibrium value in terms of a first order perturbation :math:`f^{\mathrm{EX}}`:

.. math::  
   n_{\boldsymbol{q} j} \simeq \bar{n}_{\boldsymbol{q} j}+\bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q} j}+1) \frac{\partial T}{\partial \boldsymbol{x}}\cdot f^{\mathrm{EX}}_{\boldsymbol{q} j}
 
The linearized BTE can then be written in the following form:

.. math::
   -v_{\boldsymbol{q} j}\left(\frac{\partial \bar{n}_{\boldsymbol{q} j}}{\partial T}\right) =
   \sum_{\boldsymbol{q}' j',\boldsymbol{q}'' j''}\Big[ P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\boldsymbol{q}'' j''}(f^{\mathrm{EX}}_{\boldsymbol{q} j}+f^{\mathrm{EX}}_{\boldsymbol{q}' j'}-f^{\mathrm{EX}}_{\boldsymbol{q}'' j''})
   + \frac{1}{2} P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j} (f^{\mathrm{EX}}_{\boldsymbol{q} j}-f^{\mathrm{EX}}_{\boldsymbol{q}' j'}-f^{\mathrm{EX}}_{\boldsymbol{q}'' j''} )\Big] \\\\
   + \sum_{\boldsymbol{q}' j'}  P^{\mathrm{isot}}_{\boldsymbol{q} j,\boldsymbol{q}' j'}  (f^{\mathrm{EX}}_{\boldsymbol{q} j} - f^{\mathrm{EX}}_{\boldsymbol{q}' j'}) + P^{\mathrm{be}}_{\boldsymbol{q} j} f^{\mathrm{EX}}_{\boldsymbol{q} j}

where the sum on :math:`\boldsymbol{q}'` and :math:`\boldsymbol{q}"` is performed in the Brillouin Zone.
The :math:`\mathrm{EX}` superscript of the first order perturbation :math:`f^{\mathrm{EX}}` denotes the exact solution of the BTE, to be distinguished from the approximated solutions that we will discuss later.
In this last equation the anharmonic scattering processes as well as the scattering with the isotopic impurities and the border effect are considered. 
More specifically :math:`P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\boldsymbol{q}'' j''}` is the scattering rate at the equilibrium  of a process where a phonon mode :math:`\boldsymbol{q} j` scatters by absorbing another mode :math:`\boldsymbol{q}' j'` to generate a third phonon mode :math:`\boldsymbol{q}'' j''`.
While :math:`P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j}` is the scattering rate at the equilibrium of a process where a phonon mode :math:`\boldsymbol{q}j` decays in two modes :math:`\boldsymbol{q}'j'` and :math:`\boldsymbol{q}''j''`.

The two scattering rates have the forms:

.. math::
   P^{\boldsymbol{q}'' j''}_{\boldsymbol{q} j,\boldsymbol{q}' j'} = \frac{2 \pi}{N_0 \hbar^2} \sum_{\boldsymbol{G}}
   |V^{(3)}(\boldsymbol{q} j,\boldsymbol{q}' j',-\boldsymbol{q}'' j'')|^2
   \bar{n}_{\boldsymbol{q} j}\bar{n}_{\boldsymbol{q}' j'}(\bar{n}_{\boldsymbol{q}'' j''}+1) \delta_{\boldsymbol{q}+\boldsymbol{q}' -\boldsymbol{q}'', \boldsymbol{G}}
   \delta(\hbar \omega_{\boldsymbol{q} j} +\hbar \omega_{\boldsymbol{q}' j'}-\hbar \omega_{\boldsymbol{q}'' j''})

.. math::
   P^{\boldsymbol{q}' j',\boldsymbol{q}'' j''}_{\boldsymbol{q} j} = \frac{2 \pi}{N_0 \hbar^2 } \sum_{\boldsymbol{G}}
   |V^{(3)}(\boldsymbol{q} j,-\boldsymbol{q}' j',-\boldsymbol{q}'' j'')|^2 
   \bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q}' j'}+1)(\bar{n}_{\boldsymbol{q}'' j''}+1)\delta_{\boldsymbol{q}-\boldsymbol{q}' -\boldsymbol{q}'', \boldsymbol{G}}
   \delta(\hbar \omega_{\boldsymbol{q} j}-\hbar \omega_{\boldsymbol{q}' j'}-\hbar \omega_{\boldsymbol{q}'' j''} )

with :math:`\boldsymbol{G}` the reciprocal lattice vectors.
In order to evaluate them it is necessary to compute the third derivative :math:`V^{(3)}` of  the total energy of the crystal :math:`\mathcal{E}^{tot}(\{u_{s \alpha} (\boldsymbol{R}_l) \})`, with respect to the atomic displacement :math:`u_{s \alpha} (\boldsymbol{R}_l)`, from the equilibrium position, of the s-th atom, along the :math:`\alpha` Cartesian coordinate in the crystal cell identified by the lattice vector :math:`\boldsymbol{R}_l`

.. math::
   V^{(3)}(\boldsymbol{q} j,\boldsymbol{q}' j',\boldsymbol{q}'' j'')= \frac{\partial^3 \mathcal{E}^{cell}}
   {\partial X_{\boldsymbol{q} j},\partial X_{\boldsymbol{q}' j'},\partial X_{\boldsymbol{q}'' j''}}

where :math:`\mathcal{E}^{cell}` is the energy per unit cell.
The non-dimensional quantity :math:`X_{\boldsymbol{q} j}` is defined by

.. math::
   X_{\boldsymbol{q} j}= \frac{1}{N_0}\sum_{l,s,\alpha} \sqrt{\frac{2 M_s \omega_{\boldsymbol{q} j}} {\hbar}} z^{s \alpha^*}_{\boldsymbol{q} j}  u_{s \alpha }(\boldsymbol{R}_l) e^{-i\boldsymbol{q}\cdot \boldsymbol{R}_l}

with :math:`z^{s \alpha}_{\boldsymbol{q}j} ` being the orthogonal phonon eigenmodes normalized on the unit cell and :math:`M_s` the atomic masses.
This expression of X can be used to transform the output of a density-functional code, i.e. the matrix of energy derivatives in real space :math:`\mathcal{E}(\boldsymbol{R}_l s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{R}''_{l''}s''\alpha'')` to the Fourier space.
The matrix is actually a periodic function, so it can be possible to neglect one of the Bravais lattice vector indices of such a tensor.
Note that, Quantum Espresso provides the matrix :math:`\mathcal{E}(\boldsymbol{0} s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{R}''_{l''}s''\alpha'')` while frozen phonon codes such as phonopy and related use :math:`\mathcal{E}(\boldsymbol{R}_l s\alpha,\boldsymbol{R}'_{l'} s' \alpha',\boldsymbol{0} s''\alpha'')`, i.e. set a different bravais lattice vector to zero.
Phoebe uses the latter convention at the time being.

The rate of the elastic scattering with isotopic impurities has the form:

.. math::
   P_{\boldsymbol{q} j,\boldsymbol{q}' j'}^{\mathrm{isot}} = \frac{\pi}{2 N_0} \omega_{\boldsymbol{q} j}\omega_{\boldsymbol{q}' j'}  
   \left[ \bar{n}_{\boldsymbol{q} j} \bar{n}_{\boldsymbol{q}' j'} + \frac{\bar{n}_{\boldsymbol{q} j} + \bar{n}_{\boldsymbol{q}' j'}} {2} \right ]
   \sum_{s} g^{s}_{2}   |  \sum_{\alpha} z^{s \alpha^*}_{\boldsymbol{q}j} \cdot z^{s \alpha}_{\boldsymbol{q}' j'} |^2 \delta (\omega_{\boldsymbol{q} j}- \omega_{\boldsymbol{q}' j'})

where :math:`g^s_2 = \frac{(M_s - \langle  M_s\rangle)^2}{ \langle M_s \rangle^2 }` is the average over the mass distribution of the atom of type :math:`s`.
In presence of two isotopes :math:`M_s` and :math:`M_{s'}` it can be written in terms of the concentration :math:`\epsilon` and mass change :math:`\Delta M_s= M_{s'} - M_s` :

.. math::
   g^s_2=  \epsilon(1-\epsilon)  \frac{ | \Delta M_s |}{ \langle M_s \rangle} 

with :math:`\langle M_s \rangle = M_s + \epsilon \Delta M_s`.
Eventually, in a system of finite size, :math:`P_{q j}^{\mathrm{be}}` describes the reflection of a phonon from the border:

.. math::
   P_{\boldsymbol{q} j}^{\mathrm{be}} = \frac{v_{\boldsymbol{q} j}}{L}\bar{n}_{\boldsymbol{q} j}(\bar{n}_{\boldsymbol{q} j}+1) 

where :math:`L` is the Casimir length of the sample.
The border scattering is treated in the relaxation time approximation and it results in a process in which a phonon from a specific state(:math:`\boldsymbol{q} j`) is reemitted from the surface contributing only to the equilibrium distribution.

For the sake of clarity we will contract from here on the vector :math:`\boldsymbol{q}` and branch index :math:`j` in a single mode index :math:`\nu`.
The BTE of Eq. \ref{BTE2} can be written as  a linear system in matrix form:

.. math::
   \boldsymbol{A} \boldsymbol{f}^{\mathrm{EX}}=\boldsymbol{b}

with the vector :math:`b_{\nu'} =-v_{\nu'}\hbar \omega_{\nu'} \bar{n}_{\nu'}(\bar{n}_{\nu'}+1)` and the matrix

.. math::
   A_{\nu,\nu'} = [{\sum_{\nu'',\nu'''}} (P^{\nu''}_{\nu,\nu'''} + \frac{ P_{\nu''',\nu''}^{\nu}}{2} ) + \sum_{\nu''} P^{\mathrm{isot}}_{\nu,\nu''} + P^{\mathrm{be}}_{\nu} ] \delta_{\nu,\nu'} - {\sum_{\nu''}} (  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu}  ) + P^{\mathrm{isot}}_{\nu,\nu'} 

where we have used :math:`P^{\nu', \nu''}_{\nu}=P_{\nu', \nu''}^{\nu}` from the detailed balance condition :math:`\bar{n}_{\nu}(\bar{n}_{\nu'}+1)(\bar{n}_{\nu''}+1) = (\bar{n}_{\nu}+1)\bar{n}_{\nu'}\bar{n}_{\nu''}` (valid under the assumption :math:`\hbar \omega = \hbar \omega' + \hbar \omega''`).
In this form the matrix is symmetric and positive semi-definite and it can be decomposed in :math:`\boldsymbol{A} = \boldsymbol{A}^{\mathrm{out}} +\boldsymbol{A}^{\mathrm{in}}`,
where

.. math::
   A^{\mathrm{out}}_{\nu,\nu'} = \frac{\bar{n}_{\nu}(\bar{n}_{\nu} +1)} {\tau^{\mathrm{T}}_{\nu}}\delta_{\nu,\nu'}

.. math::   
   A^{\mathrm{in}}_{\nu,\nu'} =  -  \sum_{\nu''} \left(  P^{\nu'}_{\nu,\nu''} -P^{\nu''}_{\nu,\nu'}+ P_{\nu',\nu''}^{\nu} \right )    + P^{\mathrm{isot}}_{\nu,\nu'} 

with :math:`\tau^{\mathrm{T}}_{\nu}` being the phonon relaxation time.
The :math:`\boldsymbol{A}^{\mathrm{out}}` diagonal matrix describes the depopulation of phonon states due to the scattering mechanisms while the :math:`\boldsymbol{A}^{\mathrm{in}}` matrix describes their repopulation due to the incoming scattered phonons.

The solution of the linear system in Eq. \ref{linearsyst} is obtained formally by inverting the matrix :math:`{\boldsymbol A}`.

.. math::
   {\boldsymbol f}^{\mathrm{EX}} =   \frac{1}{\boldsymbol{A}}  {\boldsymbol b}

and subsequently the thermal conductivity will be evaluated as:

.. math::
   k =  \lambda {\boldsymbol b} \cdot {\boldsymbol f}^{\mathrm{EX}}
   = - \frac{\hbar}{N_0\Omega  k_B T^2}\sum_{\nu}v_{\nu}
   \omega_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1) f_{\nu}^{\mathrm{EX}}

with :math:`\lambda= 1 /(N_0\Omega k_B T^2)`.
 

RTA solution of the phonon BTE
------------------------------

In the relaxation time approximation (RTA), we set :math:`\boldsymbol{A}^{\mathrm{in}}` to zero

.. math::
   {\boldsymbol f}^{\mathrm{SMA}} =\frac{1}{ \boldsymbol{A}^{\mathrm{out}}}  {\boldsymbol b}

Inverting :math:`\boldsymbol{A}^{\mathrm{out}}` is trivial due to its diagonal form.
The lattice thermal conductivity in RTA is then 

.. math::
   k^{\mathrm{RTA}}=\lambda \boldsymbol{b} \cdot \boldsymbol{f}^{\mathrm{SMA}}=\frac{\hbar^2}{N_0\Omega k_B T^2}\sum_{\nu}v^2_{\nu} \omega^2_{\nu} \bar{n}_{\nu}(\bar{n}_{\nu}+1)\tau^{\mathrm{T}}_{\nu}.



Iterative solution of the phonon BTE - Omini-Sparavigna method
--------------------------------------------------------------

Note: generally, we recommend the variational method over this. 

An exact solution of the BTE that does not imply either storing or the explicit inversion of matrix :math:`\boldsymbol{A}` has been proposed by Omini and Sparavigna by converging with respect to the iteration :math:`i` the following:

.. math::
   \boldsymbol{f}_{ i+1} =\frac{1} {\boldsymbol{A}^{\mathrm{out} } } \boldsymbol{b} - \frac{1} {\boldsymbol{A}^{\mathrm{out} } } \boldsymbol{A}^{\mathrm{in}}  \boldsymbol{f}_{i}

with the iteration zero consisting in the RTA :math:`\boldsymbol{f}_0=\boldsymbol{f}^{\mathrm{RTA}}`.
Instead of storing and inverting :math:`\boldsymbol{A}`, it just requires the evaluation of :math:`\boldsymbol{A}^{\mathrm{in}}\:\boldsymbol{f}_{i}`, at each iteration :math:`i` of the OS method, which is an operation computationally much less demanding.
Once the convergence is obtained the thermal conductivity is evaluated by:

.. math::
   k^{\mathrm{NV}}(\boldsymbol{f}_i)=\lambda \boldsymbol{b}\cdot \boldsymbol{f}_{i}

From a mathematical point of views the OS iterative procedure 
can be written as a geometric series:

.. math::
   \boldsymbol{f}_{ i} = \sum_{j=0,i} \left(-\frac{1}{\boldsymbol{A}^{\mathrm{out}}}  \boldsymbol{A}^{\mathrm{in}}\right)^{j} \frac{1}{\boldsymbol{A}^{\mathrm{out}}} \:  \boldsymbol{b} \;.


Iterative solution of the phonon BTE - Variational method
---------------------------------------------------------

An alternative approach consists in using the properties of the matrix :math:`{\boldsymbol A}` to find the exact solution of the linearized BTE, via the variational principle.
Indeed the solution  of the BTE is the vector :math:`\boldsymbol{f}^{\mathrm{EX}}` which makes  stationary the quadratic form

.. math::
   \mathcal{F}(\boldsymbol{f}) =\frac{1}{2} {\boldsymbol f} \cdot{\boldsymbol A} {\boldsymbol f}- {\boldsymbol b} \cdot {\boldsymbol f}

for a generic vector :math:`\boldsymbol{f}`.
Since :math:`\boldsymbol{A}` is positive the stationary point is the global and single minimum of this functional.
One can then define a variational conductivity functional:

.. math:: 
   k^\mathrm{V}(\boldsymbol{f}) = - 2 \lambda \mathcal{F}({\boldsymbol f})

that has the property :math:`k^\mathrm{V}(\boldsymbol{f}^{\mathrm{EX}})=k` while any other value of :math:`k^{\mathrm{V}}(\boldsymbol{f})`  underestimates :math:`k`.
In other words, finding the minimum of the quadratic form is equivalent to maximizing the thermal conductivity functional. 
As a consequence an error :math:`\delta \boldsymbol{f}= \boldsymbol{f} - \boldsymbol{f}^{\mathrm{EX}}`  results in an error in conductivity, linear in :math:`\delta \boldsymbol{f}` when using the non-variational estimator, and quadratic in the variational form.

Here we solve the BTE on a grid (as in OS procedure) by using the conjugate gradient method, to obtain the exact solution of the BTE equation.
In order to speed up the convergence of the conjugate gradient we take advantage of the diagonal and dominant role of :math:`\boldsymbol{A}^{\mathrm{out}}` and we use a preconditioned conjugate gradient.
Formally, this corresponds to use in the minimization the rescaled variable:

.. math::
   \tilde{{\boldsymbol f}} = \sqrt{{\boldsymbol A^{\mathrm{out}}}} {\boldsymbol f}


and then, with respect to this new variable, minimize the quadratic form :math:`\tilde{\mathcal{F}}(\tilde{\boldsymbol{f}}) = \mathcal{F}(\boldsymbol{f})` where:

.. math::
   \tilde{\mathcal{F}}( \tilde{\boldsymbol{f}}) =\frac{1}{2} \tilde{\boldsymbol{f}}\cdot \tilde{\boldsymbol{A}} \tilde{\boldsymbol{f}}- \tilde{\boldsymbol{ b}}\cdot\tilde{\boldsymbol {f}}

and  

.. math::
   \tilde{{\boldsymbol A}} =\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}} {\boldsymbol A}\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}}

.. math::
   \tilde{{\boldsymbol b}} =\frac{1}{ \sqrt{{\boldsymbol A^{\mathrm{out}}}}} {\boldsymbol b} \label{prec3}


Notice that :math:`\tilde{\boldsymbol{f}}^{\mathrm{RTA}}=\tilde{\boldsymbol{b}}`.
The square root evaluation of :math:`\boldsymbol{A}^{\mathrm{out}}` is trivial due to its diagonal form.
The computational cost per iteration of the conjugate gradient scheme is equivalent to the OS one, but it always converges and requires a smaller number of iterations.


The conjugate gradient minimization requires the evaluation of the gradient :math:`\boldsymbol{g}_i= \boldsymbol{A} \boldsymbol{f}_i - \boldsymbol{b}` and a line minimization.
Since the form is quadratic the line minimization can be done analytically and exactly.
Moreover the information required by the line minimization at  iteration :math:`i` can be recycled to compute the gradient at the next iteration :math:`i+1`.
Starting with an the initial vector :math:`\boldsymbol{f}_0= \boldsymbol{f}^{\mathrm{RTA}}`, initial gradient :math:`\boldsymbol{g}_0=\boldsymbol{A}\boldsymbol{f}_0 -\boldsymbol{f}^{\mathrm{RTA}}` and letting :math:`\boldsymbol{h}_0= -\boldsymbol{g}_0`, the conjugate gradient method can be summarized with the
recurrence:

.. math::
   \boldsymbol{t}_i =\boldsymbol{A} \boldsymbol{h}_i

.. math::
   {\boldsymbol f}_{i+1} = {\boldsymbol f}_{i} - \frac{\boldsymbol {g}_{i} \cdot {\boldsymbol{h}_{i}} } {\boldsymbol{h}_{i} \cdot \boldsymbol{t}_i } \boldsymbol{h}_{i}
  
.. math::
   \boldsymbol{g}_{i+1} = \boldsymbol{g}_{i}-\frac{\boldsymbol {g}_{i} \cdot {\boldsymbol{h}_{i}} } {\boldsymbol{h}_{i} \cdot \boldsymbol{t}_i }\boldsymbol{t}_i

.. math::
   \boldsymbol{h}_{i+1} = -\boldsymbol{g}_{i+1} + \frac{\boldsymbol{g}_{i+1} \cdot \boldsymbol{g}_{i+1}}{{\boldsymbol{g}_{i}} \cdot {\boldsymbol{g}_{i}} }  {\boldsymbol h}_{i} 

where :math:`\boldsymbol{h}_i` is the search direction and :math:`\boldsymbol{t}_i` is an auxiliary vector.
Notice that each iteration requires only one application of the matrix :math:`\boldsymbol{A}` on the vector :math:`\boldsymbol{h}_i` as in the OS method.




Relaxons solution to the BTE
----------------------------

We first diagonalize the scattering matrix:

.. math::
   \frac{1}{N_k} \sum_{\nu} \Omega_{\nu\nu'} \theta_{\nu'\alpha} = \frac{1}{\tau_{\alpha}} \theta_{\nu\alpha}

where :math:`\theta` are eigenvectors, :math:`\alpha` are eigenvalue indices, :math:`\frac{1}{\tau_{\alpha}}` are eigenvalues, and the scattering matrix is:

.. math::
   \Omega_{\nu\nu'} = \frac{ A_{\nu\nu'} } { \sqrt{ \bar{n}_{\nu}(\bar{n}_{\nu}+1) \bar{n}_{\nu'}(\bar{n}_{\nu'}+1)  } }

Next, we compute the velocities:

.. math::
   \boldsymbol{V}_{\alpha} = \frac{1}{N_k} \sum_{\nu} \theta_{\nu0} \boldsymbol{v}_{\nu} \theta_{\nu\alpha}

where

.. math::
   \theta_{\nu0} = \sqrt { \frac{ \frac{\partial \bar{n}_{\nu}}{\partial \epsilon} }{C T} } \hbar \omega_{\nu}

Finally, the thermal conductivity is:

.. math::
   k^{ij} = \sum_{\alpha} C V_{\alpha}^i V_{\alpha}^j \tau_{\alpha}







Phonon velocity operator
------------------------

The velocity operator matrix elements (e.g. along the x direction) can be computed using the Hellmann-Feynman theorem from the Dynamical matrix :math:`\boldsymbol{\mathcal{D}}`:

.. math::
   V^x_{j j'}(\boldsymbol{q}) = \sum_{\alpha \alpha' s s'} \frac{1}{2 \sqrt{M_s M_{s'}} \omega_{\boldsymbol{q} j} }  z^{s \alpha^*}_{\boldsymbol{q} j}  \frac{\partial \mathcal{D}^{\alpha \alpha'}_{s s'}(\boldsymbol{q})}{ \partial q_x}   z^{s' \alpha'}_{\boldsymbol{q} j'}

In the non-degenerate case, the group velocity is :math:`\boldsymbol{v}_{\boldsymbol{q} j}=\boldsymbol{V}_{j j}(\boldsymbol{q})` while in the degenerate one we use the phonon polarization vectors that diagonalize the matrix in the degenerate subspace.



Dirac-delta approximations
--------------------------

The delta function for the energy conservation can be replaced by a Gaussian

.. math::
   \delta(\hbar \omega)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar \omega/ \sigma )^2)} \;,

where :math:`\sigma` is a constant decided by user input.
It is important to note that when the delta function is substituted with a Gaussian the detailed balance condition is only valid under approximation.
The definition used above guarantees that the scattering matrix is symmetric and non-negative. 

Another method is the adaptive-gaussian smearing scheme(https://link.aps.org/doi/10.1103/PhysRevB.75.195121).
Specifically, we want to approximate a dirac-delta function of the form:

.. math::
   \delta(\hbar (\omega_1+\omega_2-\omega_3)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar (\omega_1+\omega_2-\omega_3)/ \sigma )^2)}.

This time, we allow :math:`\sigma` to be a value dependent on the energies.
Specifically, we build it as:

.. math::
   \sigma = \frac{1}{\sqrt{12}} \sqrt{ \sum_{\beta} (\sum_{\alpha} (v_2-v_3) \frac{M_{\alpha \beta}}{N_{\beta}}  )^2 }

where :math:`M` is a matrix comprised of the primitive cell lattice vectors (each column is a lattice vector) ,:math:`v_2` and :math:`v_3` are phonon group velocities, and :math:`N_{\beta}` is the number of wavevectors sampled along direction :math:`\beta`.
Note that the adaptive scheme may be critical in the case where the velocity sum to zero: in that case, we skip the scattering event, unless we have an exact energy conservation taking place.





Wigner correction to phonon thermal conductivity
------------------------------------------------

The theory is fully described in the Reference available at this link (https://www.nature.com/articles/s41567-019-0520-x).

In extreme synthesis, the thermal conductivity is estimated as:

.. math::
   k_{\alpha\beta} = k^{BTE}_{\alpha\beta} +  \frac{k_BT^2}{\Omega N_k} \sum_{\boldsymbol{q}} \sum_{s\neq s'} \frac{\omega_{\boldsymbol{q}j}+\omega_{\boldsymbol{q}j'}}{2}   V_{jj'}^{\alpha}(\boldsymbol{q}) V_{j'j}^{\beta}(\boldsymbol{q}) \frac{ ( \frac{\partial n_{\boldsymbol{q}j}}{\partial T} + \frac{\partial n_{\boldsymbol{q}j'}}{\partial T})(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'}) }{4(\omega_{\boldsymbol{q}j}-\omega_{\boldsymbol{q}j'})^2 + (\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})^2} 


where :math:`k^{BTE}_{\alpha\beta}` is the thermal conductivity estimated by the Boltzmann transport equation discussed above, and :math:`\Gamma_{\boldsymbol{q}j} = \frac{1}{\tau_{\boldsymbol{q}j}}` is the phonon linewidth, i.e. a diagonal element of the scattering matrix.


Thermal Viscosity
-----------------

The theory is described to far greater extent in this reference https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011019.
The equilibrium of a system of bosonic particles that conserves energy and momentum is the drifting distribution:

.. math::
   n_{\nu}^{D}
   =
   \frac{1}{e^{\beta(\hbar \omega_\nu - \hbar \boldsymbol{q} \cdot \boldsymbol{u})}-1} \;,

where :math:`\boldsymbol{q}` is the phonon wavevector (proportional to the phonon crystal momentum, and :math:`\boldsymbol{u}` is the phonon drift velocity.
The thermal viscosity is defined as the coefficient of proportionality between the crystal momentum flux :math:`\Pi` and a local perturbation in the drift velocity :math:`\boldsymbol{u}`.

.. math::
   \Pi^{ij} = - \sum_{kl} \eta^{ijkl} \frac{\partial u^k}{\partial r^l}

and the momentum flux (at least, the component relevant to our case) is defined as:

.. math::
   \Pi^{ij} = \frac{1}{V N_q} \sum_{\nu} \hbar q^i v_{\nu}^j n_{\nu}

The population in response to the perturbation is fixed by the phonon BTE.
At the RTA level, we simply need to solve

.. math::
   \boldsymbol{v}_{\nu} \cdot (\frac{\partial n^{D}_{\nu}}{\partial \boldsymbol{u}} \cdot \nabla \boldsymbol{u} )
   = - \frac{n_{\nu}}{\tau_{\nu}}

We linearize the solution, stating :math:`n_{\nu} = n_{\nu} \nabla \boldsymbol{u}`, and the equation is readily solved.
Beyond the RTA, we must solve the equation:

.. math::
   \frac{\boldsymbol{v}_{\nu}}{\sqrt{\bar{n}_{\nu}(\bar{n}_{\nu}+1)}} \cdot (\frac{\partial n^{D}_{\nu}}{\partial \boldsymbol{u}} \cdot \nabla \boldsymbol{u} )
   = - \frac{1}{V N_q} \sum_{\nu'} \Omega_{\nu\nu'} n_{\nu'}

which we do with the eigenvector formalism.
Using the eigenvectors of the scattering matrix, we expand the phonon population as:

.. math::
   n_{\nu} = \sum_{kl} f^{kl}_{\alpha} \theta_{\nu\alpha} \nabla_l u^k

We find the solution as:

.. math::
   f^{kl}_{\alpha} = - \tau_{\alpha} \sum_{\nu} \theta_{\nu\alpha} \frac{\boldsymbol{v}^l_{\nu}}{\sqrt{\bar{n}_{\nu}(\bar{n}_{\nu}+1)}} \frac{\partial n^{D}_{\nu}}{\partial u^k} 

which can be used to reconstruct the phonon population response.
Finally, the viscosity tensor is symmetrized, finding the thermal viscosity:

.. math::
   \mu^{ijkl} = \frac{1}{2} ( \eta^{ijkl} + \eta^{ilkj} )

The code also prints other quantities needed to write the viscous heat equations derived in this reference https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011019.


Symmetries of the BTE
---------------------

We exploit the symmetries of the crystal to speed-up the calculation of thermal conductivity.
We took inspiration from this reference https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.265506.
Let :math:`q` indicate any wavevector in the Brillouin zone.
The symmetries of a crystal identify an irreducible set of wavevectors :math:`q^*`, such that any other wavevector :math:`q` can be obtained from a rotation of these irreducible wavevectors :math:`q = R q^*`.
The basic idea is to restrict the calculation to the irreducible set of wavevectors.
The conductivity for example, is:

.. math::
   k^{ij}
   = \frac{1}{V N_k} \sum_{\nu} \hbar \omega_{\nu} v^i_{\nu} n^{j}_{\nu}
   = \frac{1}{V N_k} \sum_{\nu^*} \sum_{R} \hbar \omega_{\nu^*} (R v_{\nu})_{i} (R n_{\nu})_{j}

where :math:`R` is the set of rotations used to reconstruct all the symmetry-equivalent wavevectors of :math:`q^*`, and the summation over :math:`\nu^*` is only done in the irreducible set of wavevectors.

The BTE too can be restricted to the irreducible wedge.

.. math::
   v^i_{\nu^*} \frac{\partial \bar{n}_{\nu}}{\partial T}
   = - \frac{1}{V N_q} \sum_{\nu'} A_{\nu^*\nu'} f^i_{\nu'}
   = - \frac{1}{V N_q} \sum_{\nu'^*} \sum_{R} \sum_{j} A_{\nu^*\nu'^*} R_{ij} f^j_{\nu'}
   = - \frac{1}{V N_q} \sum_{\nu'^* j} A^{ij}_{\nu^*\nu'^*} f^j_{\nu'}

Hence, one can work with the same techniques detailed above, provided that we work with an enlarged matrix :math:`A^{ij}_{\nu^*\nu'^*}`.

By default, we make use of symmetries.
Some comments:

* Advantage: for a system with a lot of symmetries, the matrix :math:`A^{ij}_{\nu^*\nu'^*}` is generally smaller than :math:`A_{\nu\nu'}`, and thus calculations will be much faster.
  
* Disadvantage 1: we cannot compute viscosity beyond the RTA using symmetries. To do so, one must disable symmetries.

* Disadvantage 2: note that the symmetric matrix gains two indices on cartesian coordinates. As a result, in the limit case that there are no symmetries in the system (only the identity), the matrix :math:`A^{ij}_{\nu^*\nu'^*}` will still be computed on the same number of wavevectors of  :math:`A_{\nu\nu'}`, but occupies 3x3 times more memory without adding any information. Therefore, for low-symmetry systems, consider disabling symmetries.

