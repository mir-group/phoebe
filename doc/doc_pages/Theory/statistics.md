@page Theory Theory
@section Statistics Statistics

Equations available in class Particle.
The Fermi--Dirac distribution (for electrons) for a Bloch state (k,b) (k,b are Bloch indeces) is

\begin{equation}
f_{k,b} = \frac{1}{e^{(\frac{\epsilon_{k,b}-\mu}{k_BT})}+1}
\end{equation}

The Bose--Einstein distribution (for phonons) is:
\begin{equation}
n_{k,b} = \frac{1}{e^{(\frac{\epsilon_{k,b}}{k_BT})}-1}
\end{equation}

A computationally stable way to evaluate derivatives is:

\begin{equation}
\frac{\partial n_{k,b}}{\partial T} = \frac{\epsilon_{k,b}}{4k_BT^2} \sinh( \frac{\epsilon_{k,b}}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial f_{k,b}}{\partial T} = \frac{\epsilon_{k,b}-\mu}{4k_BT^2} \cosh( \frac{\epsilon_{k,b}-\mu}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial n_{k,b}}{\partial \epsilon} = - \frac{1}{4k_BT} \sinh( \frac{\epsilon_{k,b}}{2T} ) 
\end{equation}

\begin{equation}
\frac{\partial f_{k,b}}{\partial \epsilon} = - \frac{1}{4k_BT} \cosh( \frac{\epsilon_{k,b}-\mu}{2T} ) 
\end{equation}
