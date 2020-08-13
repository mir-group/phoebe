@page Theory Theory
@section WIGNERPH Wigner correction to phonon thermal conductivity.

The theory is fully described in the Reference available at this [link](https://www.nature.com/articles/s41567-019-0520-x).

In extreme synthesis, the thermal conductivity is estimated as:
\begin{equation}
k_{\alpha\beta} = k^{BTE}_{\alpha\beta} +  \frac{k_BT^2}{\Omega N_k} \sum_{q} \sum_{s\neq s'} \frac{\omega_{qj}+\omega_{qj'}}{2}   V_{jj'}^{\alpha}(q) V_{j'j}^{\beta}(q) \frac{ ( \frac{\partial n_{qj}}{\partial T} + \frac{\partial n_{qj'}}{\partial T})(\Gamma_{qj}+\Gamma_{qj'}) }{4(\omega_{qj}-\omega_{qj'})^2 + (\Gamma_{qj}+\Gamma_{qj'})^2} 
\end{equation}

where \f$k^{BTE}_{\alpha\beta}\f$ is the thermal conductivity estimated by the Boltzmann transport equation discussed above, and \f$\Gamma_{qj} = \frac{1}{\tau_{qj}}\f$ is the phonon linewidth, i.e. a diagonal element of the scattering matrix.