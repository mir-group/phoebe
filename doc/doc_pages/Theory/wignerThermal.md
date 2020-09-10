@page Theory Theory
@section WIGNERPH Wigner correction to phonon thermal conductivity.

The theory is fully described in the Reference available at this [link](https://www.nature.com/articles/s41567-019-0520-x).

In extreme synthesis, the thermal conductivity is estimated as:
\begin{equation}
k_{\alpha\beta} = k^{BTE}_{\alpha\beta} +  \frac{k_BT^2}{\Omega N_k} \sum_{\boldsymbol{q}} \sum_{s\neq s'} \frac{\omega_{\boldsymbol{q}j}+\omega_{\boldsymbol{q}j'}}{2}   V_{jj'}^{\alpha}(\boldsymbol{q}) V_{j'j}^{\beta}(\boldsymbol{q}) \frac{ ( \frac{\partial n_{\boldsymbol{q}j}}{\partial T} + \frac{\partial n_{\boldsymbol{q}j'}}{\partial T})(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'}) }{4(\omega_{\boldsymbol{q}j}-\omega_{\boldsymbol{q}j'})^2 + (\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})^2} 
\end{equation}

where \f$k^{BTE}_{\alpha\beta}\f$ is the thermal conductivity estimated by the Boltzmann transport equation discussed above, and \f$\Gamma_{\boldsymbol{q}j} = \frac{1}{\tau_{\boldsymbol{q}j}}\f$ is the phonon linewidth, i.e. a diagonal element of the scattering matrix.