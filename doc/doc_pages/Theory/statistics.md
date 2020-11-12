@page thStat Independent particle properties

@section thPartStat Particle statistics

Equations available in class Particle.
The Fermi--Dirac distribution (for electrons) for a Bloch state (\f$\boldsymbol{k}\f$,\f$b\f$) (\f$\boldsymbol{k}\f$,\f$b\f$ are Bloch indeces) is

\begin{equation}
f_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}-\mu}{k_BT})}+1}
\end{equation}

The Bose--Einstein distribution (for phonons) is:
\begin{equation}
n_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}}{k_BT})}-1}
\end{equation}

A computationally stable way to evaluate derivatives is:

\begin{equation}
\frac{\partial n_{\boldsymbol{k},b}}{\partial T} = \frac{\epsilon_{\boldsymbol{k},b}}{4k_BT^2} \sinh( \frac{\epsilon_{\boldsymbol{k},b}}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial f_{\boldsymbol{k},b}}{\partial T} = \frac{\epsilon_{\boldsymbol{k},b}-\mu}{4k_BT^2} \cosh( \frac{\epsilon_{\boldsymbol{k},b}-\mu}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial n_{\boldsymbol{k},b}}{\partial \epsilon} = - \frac{1}{4k_BT} \sinh( \frac{\epsilon_{\boldsymbol{k},b}}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial f_{\boldsymbol{k},b}}{\partial \epsilon} = - \frac{1}{4k_BT} \cosh( \frac{\epsilon_{\boldsymbol{k},b}-\mu}{2T} ) 
\end{equation}


@section thCV Specific heat

The specific heat at constant volume for phonons (electrons) is evaluated as:

\begin{equation}
C_v = \frac{1}{V N_k} \sum_{\boldsymbol{k},b} \epsilon_{\boldsymbol{k},b} \frac{\partial n_{\boldsymbol{k},b}}{\partial T}
\end{equation}

where \f$V\f$ is the volume of the primitive crystal unit cell, \f$N_k\f$ is the number of wavevectors used, \f$\boldsymbol{k}\f$ is a wavevector index, \f$b\f$ is a branch/band index, \f$\epsilon\f$ is the phonon(electron, w.r.t. the chemical potential) energy, \f$T\f$ the temperature and \f$n\f$ is the Bose--Einstein (Fermi--Dirac) distribution function.





@section thDos Density of States

The density of states is defined as the number of states that are available to a particle at a certain energy \f$ E \f$:

\begin{equation}
DOS(E) = \frac{1}{(2\pi)^d V} \sum_b \int_{BZ} d\boldsymbol{k} \delta(E-\epsilon_{\boldsymbol{k}b})
\end{equation}
where \f$d\f$ is the dimensionality, \f$ V \f$ is the crystal unit cell volume, \f$ b \f$ is a Bloch index over bands, \f$ \boldsymbol{k} \f$ is a Bloch index over wavevectors, \f$ \epsilon \f$ is the particle energy and the integration is carried over the Brillouin zone.

The definition holds for both phonon and electrons.

The integral is sampled with a uniform grid of points in the Brillouin zone.
In the DoSApps, the integral of the Dirac-delta is implemented using the tetrahedron method.