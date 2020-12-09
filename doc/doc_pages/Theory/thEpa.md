@page thEpa Electron-phonon averaging (EPA)

\section thEppa Electron-phonon averaging

In phoebe, we implemented the electron-phonon average method discuss in ref. https://doi.org/10.1016/j.mtphys.2018.07.001.
The purpose of this method is to have a fast tool for the computation of electronic transport coefficients limited by electron-phonon scattering.
The benefits of this methodology are that EPA:
1. avoids Wannierization: while accurate, finding Wannier functions is tricky and is not yet fully automatizable; and
2. gives priority to speed rather than accuracy, allowing a bigger computational throughput.

We only consider transport coefficients within the relaxation time approximation.
In such case, the electronic lifetimes due to the electron-phonon scattering are:

\f{eqnarray}
\tau_{n\boldsymbol{k}}^{-1}(\mu,T)
&=&
\frac{V}{(2\pi)^2 \hbar} \sum_{m\nu}
\int_{BZ} d\boldsymbol{q}
|g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\bigg[ \big(n(\omega_{\boldsymbol{q}\nu},T) + f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} + \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})  \\
&& \big(n(\omega_{\boldsymbol{q}\nu},T) + 1 - f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} - \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m}) \bigg]
\f}
While accurate, this expression is still too expensive for some purposes.

The main idea of EPA is to replace integrals over the Brillouin zone with some averaged quantities.
In particular, rather than working with the full phonon dispersion, we average the frquencies of each phonon mode: 
\begin{equation}
\omega_{\boldsymbol{q}\nu}
\to
\bar{\omega}_{\nu}
\end{equation}
So that we have \f$ 3 N_{atoms} \f$ modes per crystal.
Similarly, the make an energy binning on the electron-phonon coupling strength, as
\begin{equation}
|g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\to
g^2_{\nu} (\epsilon_{\boldsymbol{k}n}, \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})
\end{equation}
There can be various ways of doing this averaging procedure on \f$g\f$.
In Phoebe, we implemented a gaussian quadrature procedure.
In this approach, the electron-phonon is approximated as a weighted sum of the electron-phonon coupling computed from an ab-initio code:
\begin{equation}
g^2_{\nu} (\epsilon_1,\epsilon_2)
=
\frac{1}{W}
\sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}} |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\end{equation}
where 
\begin{equation}
W = \sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}}
\end{equation}
and the weights \f$ w \f$ are found by minimizing
\begin{equation}
\sum_{mn} \int d\boldsymbol{q} d\boldsymbol{k} ( g^2_{\nu} (\epsilon_1,\epsilon_2) - |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2 )^2
\exp\bigg( -\frac{(\epsilon_{\boldsymbol{k}n}-\epsilon_1)^2+(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m}-\epsilon_2)^2}{2\sigma^2_{gauss}} \bigg)
\end{equation}

After the electron-phonon coupling is approximated in this way, the electron lifetimes can be integrated as:
\f{eqnarray}
\tau^{-1}(\epsilon,\mu,T)
=
\frac{2\pi V}{g_s \hbar} \sum_{\nu}
&& g^2_{\nu}(\epsilon,\epsilon+\bar{\omega}_{\nu})
\big(n(\bar{\omega}_{\nu},T) + f(\epsilon + \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon + \bar{\omega}_{\nu})  +  \\
&& g^2_{\nu}(\epsilon,\epsilon-\bar{\omega}_{\nu})
\big(n(\bar{\omega}_{\nu},T) + f(\epsilon - \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon - \bar{\omega}_{\nu})
\f}
where
\begin{equation}
v^2_{\alpha\beta} (\epsilon) \rho (\epsilon)
=
\sum_n \int_{BZ} d\boldsymbol{k} v_{n\boldsymbol{k}\alpha} v_{n\boldsymbol{k}\beta} \delta(\epsilon-\epsilon_{\boldsymbol{k}n})
\end{equation}

Finally, transport coefficients can be computed from the following tensor:
\begin{equation}
K_{\alpha\beta}^{(p)}
=
\frac{g_s e^{2-p}}{(2\pi)^{3} (k_B T)^{p+1}} \int d\epsilon v^2_{\alpha\beta}(\epsilon) \rho(\epsilon) \tau(\epsilon,\mu,T) (\epsilon-\mu)^p f(\epsilon,\mu,T) [1-f(\epsilon,\mu,T)]
\end{equation}
For example, the electrical conductivity is:
\begin{equation}
\sigma_{\alpha\beta}(\mu,T) = K_{\alpha\beta}^{(0)}
\end{equation}

Note that, in order to compute the velocity term, we still need a form of interpolation of the electronic band structure.
In order to keep the spirit of avoiding the Wannierization procedure, we use the @ref thFourierElectrons.

Furthermore, we note that the density of states is integrated using the tetrahedron method.

In terms of computational parameters (besides temperature and doping concentration) the EPA requires tuning a few parameters, namely
1. the size of the coarse k/q grid used in an ab-initio code to compute the coupling; 
2. the energy bins used to approximate the electron-phonon coupling (which are set when converting data from Quantum ESPRESSO to Phoebe), 
3. the energy bins used to integrate lifetimes (typically, one should use a better energy sampling than the one for the coupling)
4. the grid of wavevectors, used to average velocities and density of states.
