@page thStat Particle statistics

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