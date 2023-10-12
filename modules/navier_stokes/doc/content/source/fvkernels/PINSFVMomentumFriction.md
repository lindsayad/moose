# PINSFVMomentumFriction

This kernel adds the friction term to the porous media Navier Stokes momentum
equations. This kernel must be used with the canonical PINSFV variable set,
e.g. pressure and superficial velocity. This kernel supports Darcy and
Forchheimer friction models:

Darcy drag model
\begin{equation}
\epsilon F_i = - f_i \rho \frac{v_{D,i}}{\epsilon}
\end{equation}
Forchheimer drag model
\begin{equation}
\epsilon F_i = - f_i \rho \frac{v_{D,i}}{\epsilon}\frac{|v_D|}{\epsilon}
\end{equation}
where $F_i$ is the i-th component of the friction force (denoted by $\mathbf{F_f}$ in [!eqref](pinsfv.md#eq:pinsfv_mom)), $f_i$ the friction factor, which may be anisotropic,
$\epsilon$ the porosity and $\rho$ the fluid density and $v_{D,i}$ the i-th
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
represent inertial effects and should have a quadratic dependence on velocity. These additional dependencies, as well as other prefactors commonly present in the Darcy and Forchheimer models, are baked into the definition of the friction factors.

## Computation of friction factors and pre-factors

To outline how friction factors for Darcy and Forchheimer may be calculated,
let's consider a specific example. We'll draw from the Ergun equation, which is
outlined [here](https://en.wikipedia.org/wiki/Ergun_equation). Let's consider
the form:

\begin{equation}
\Delta p = \frac{150\mu L}{d_p^2} \frac{(1-\epsilon)^2}{\epsilon^3} v_D + \frac{1.75 L \rho}{d_p} \frac{(1-\epsilon)}{\epsilon^3} |v_D| v_D
\end{equation}

where $L$ is the bed length, $\mu$ is the fluid dynamic viscosity and $d_p$ is
representative of the diameter of the pebbles in the pebble bed. We can divide
the equation through by $L$, recognize that $\Delta p$ denotes $p_0 - p_L$ such
that $\Delta p/L \leftarrow -\nabla p$, multiply the equation through by
$-\epsilon$, move all terms to the left-hand-side, and do
some term manipulation in order to yield:

\begin{equation}
\epsilon \nabla p + 150 \frac{\mu\epsilon}{\rho}
\frac{(1-\epsilon)^2}{\epsilon^2 d_p^2} \frac{\rho v_D}{\epsilon} + 1.75\epsilon
\frac{1-\epsilon}{\epsilon d_p} \rho \frac{|\vec{v}_D|}{\epsilon} \frac{v_D}{\epsilon} = 0
\end{equation}

If we define the hydraulic diameter as $D_h = \frac{\epsilon d_p}{1 -
\epsilon}$, then the above equation can be rewritten as:

\begin{equation}
\epsilon \nabla p + \frac{150\mu\epsilon}{\rho D_h^2}
 \frac{\rho v_D}{\epsilon} + \frac{1.75\epsilon}{D_h}
\rho \frac{|\vec{v}_D|}{\epsilon} \frac{v_D}{\epsilon} = 0
\end{equation}

From this equation we can see that the Darcy coefficient is computed via
\begin{equation}
\frac{150\mu\epsilon}{\rho D_h^2}
\end{equation}

and the Forchheimer coefficient is computed via
\begin{equation}
\frac{1.75\epsilon}{D_h}
\end{equation}

!syntax parameters /FVKernels/PINSFVMomentumFriction

!syntax inputs /FVKernels/PINSFVMomentumFriction

!syntax children /FVKernels/PINSFVMomentumFriction
