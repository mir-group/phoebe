@page thEpa Electron-phonon averaging (EPA)

\section thEppa Electron-phonon averaging

Something

\f{eqnarray}
\tau_{n\boldsymbol{k}}^{-1}(\mu,T)
&=&
\frac{V}{(2\pi)^2 \hbar} \sum_{m\nu}
\int_{BZ} d\boldsymbol{q}
|g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\bigg[ \big(n(\omega_{\boldsymbol{q}\nu},T) + f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} + \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})  \\
&& \big(n(\omega_{\boldsymbol{q}\nu},T) + 1 - f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} - \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m}) \bigg]
\f}


\begin{equation}
v^2_{\alpha\beta} (\epsilon) \rho (\epsilon)
=
\sum_n \int_{BZ} d\boldsymbol{k} v_{n\boldsymbol{k}\alpha} v_{n\boldsymbol{k}\beta} \delta(\epsilon-\epsilon_{\boldsymbol{k}n})
\end{equation}


\f{eqnarray}
\tau^{-1}(\epsilon,\mu,T)
=
\frac{2\pi V}{g_s \hbar} \sum_{\nu}
&& g^2_{\nu}(\epsilon,\epsilon+\bar{\omega}_{\nu})
\big(n(\bar{\omega}_{\nu},T) + f(\epsilon + \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon + \bar{\omega}_{\nu})  +  \\
&& g^2_{\nu}(\epsilon,\epsilon-\bar{\omega}_{\nu})
\big(n(\bar{\omega}_{\nu},T) + f(\epsilon - \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon - \bar{\omega}_{\nu})
\f}


\begin{equation}
K_{\alpha\beta}^{(p)}
=
\frac{g_s e^{2-p}}{(2\pi)^{3} (k_B T)^{p+1}} \int d\epsilon v^2_{\alpha\beta}(\epsilon) \rho(\epsilon) \tau(\epsilon,\mu,T) (\epsilon-\mu)^p f(\epsilon,\mu,T) [1-f(\epsilon,\mu,T)]
\end{equation}

\begin{equation}
\sigma_{\alpha\beta}(\mu,T) = K_{\alpha\beta}^{(0)}
\end{equation}



\begin{equation}
|g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\to
g^2_{\nu} (\epsilon_{\boldsymbol{k}n}, \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})
\end{equation}

\begin{equation}
\omega_{\boldsymbol{q}\nu}
\to
\bar{\omega}_{\nu}
\end{equation}


\begin{equation}
g^2_{\nu} (\epsilon_1,\epsilon_2)
=
\frac{1}{V_1}
\sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}} |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
\end{equation}


\begin{equation}
V_1 = \sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}}
\end{equation}


\begin{equation}
\sum_{mn} \int d\boldsymbol{q} d\boldsymbol{k} ( g^2_{\nu} (\epsilon_1,\epsilon_2) - |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2 )^2
\exp\bigg( -\frac{(\epsilon_{\boldsymbol{k}n}-\epsilon_1)^2+(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m}-\epsilon_2)^2}{2\sigma^2_{gauss}} \bigg)
\end{equation}
