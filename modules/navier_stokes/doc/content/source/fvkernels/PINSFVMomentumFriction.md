# PINSFVMomentumFriction

This kernel adds the friction term to the porous media Navier Stokes momentum
equations. This kernel must be used with the canonical PINSFV variable set,
e.g. pressure and superficial velocity. This kernel supports Darcy and
Forchheimer friction models:

Darcy drag model
\begin{equation}
\epsilon F_i = - f_i \frac{\rho v_i}{\epsilon}
\end{equation}
Forchheimer drag model
\begin{equation}
\epsilon F_i = - f_i \frac{\rho v_i}{\epsilon}
\end{equation}
where $F_i$ is the i-th component of the friction force (denoted by $\mathbf{F_f}$ in [!eqref](pinsfv.md#eq:pinsfv_mom)), $f_i$ the friction factor, which may be anisotropic,
$\epsilon$ the porosity and $\rho$ the fluid density and $v_i$ the i-th
component of the fluid
superficial velocity. We have used a negative sign to match the notation used in
[!eqref](pinsfv.md#eq:pinsfv_mom) where the friction force is on the
right-hand-side of the equation. When moved to the left-hand side, which is done
when setting up a Newton scheme, the term becomes positive which is what is
shown in the source code itself.
Though the functional forms above are identical and their
treatment in this kernel is identical, there is a fundamental conceptual
difference. Darcy is meant to represent viscous effects and should
have a linear dependence on the fluid velocity, whereas Forchheimer is meant to
represent inertial effects and should have a quadratic dependence on velocity
(as we will show below).

To outline how friction factors for Darcy and Forchheimer may be calculated,
let's consider a specific example. We'll draw from the Ergun equation, which is
outlined [here](https://en.wikipedia.org/wiki/Ergun_equation). Let's consider
the form:

\begin{equation}
\Delta p = \frac{150\mu L}{d_p^2} \frac{(1-\epsilon)^2}{\epsilon^3} v_i + \frac{1.75 L \rho}{d_p} \frac{(1-\epsilon)}{\epsilon^3} |\vec{v}| v_i
\end{equation}

where $L$ is the bed length, $\mu$ is the fluid dynamic viscosity and $d_p$ is
representative of the diameter of the pebbles in the pebble bed. We can divide
the equation through by $L$, multiply the equation through by $\epsilon$, and do
some term manipulation in order to yield:

\begin{equation}
0 = -\epsilon \nabla p + \left(-150 \frac{\mu\epsilon}{\rho} \frac{1-\epsilon}{\epsilon d_p}^2 \frac{\rho v_i}{\epsilon} - 1.75 \frac{1-\epsilon}{\epsilon d_p} |\vec{v}| \frac{\rho v_i}{\epsilon}\right)
\end{equation}

If we define the hydraulic diameter as $D_h = \frac{\epsilon d_p}{1 -
\epsilon}$, then the above equation can be rewritten as:

\begin{equation}
0 = -\epsilon \nabla p + \left(-150 \frac{\mu\epsilon}{\rho} \frac{1}{D_h^2} \frac{\rho v_i}{\epsilon} - 1.75 \frac{1}{D_h} |\vec{v}| \frac{\rho v_i}{\epsilon}\right)
\end{equation}

From this equation we can see that the Darcy coefficient is computed via
\begin{equation}
150 \frac{\mu\epsilon}{\rho} \frac{1}{D_h^2}
\end{equation}

and the Forchheimer coefficient is computed via
\begin{equation}
1.75 \frac{1}{D_h} |\vec{v}|
\end{equation}

!syntax parameters /FVKernels/PINSFVMomentumFriction

!syntax inputs /FVKernels/PINSFVMomentumFriction

!syntax children /FVKernels/PINSFVMomentumFriction
