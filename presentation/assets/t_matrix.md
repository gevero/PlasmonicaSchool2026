# The Transfer Matrix Method for Optical Multilayers

The Transfer Matrix Method (TMM) is a powerful, elegant, and widely used mathematical formalism for calculating the optical properties (reflectance, transmittance, and absorptance) of multilayered structures.

Unlike the Scattering Matrix method (which maps incoming waves to outgoing waves), the Transfer Matrix Method maps the total wave field from one side of a structure directly to the other side. By multiplying a sequence of $2 \times 2$ matrices representing the interfaces and the layers, we can relate the macroscopic optical fields of the ambient medium to the fields in the substrate.

This document provides a detailed, self-contained derivation of the TMM for an optical multilayer structure.

## 1. Defining the System and the Waves

Consider a stack of $N$ homogeneous, isotropic, and parallel layers sandwiched between a semi-infinite ambient medium (Layer $0$) and a semi-infinite substrate (Layer $N+1$).

Let the $z$-axis be perpendicular to the layers, pointing downwards. The boundaries between the layers are located at $z_1, z_2, \dots, z_N$. The thickness of the $i$-th layer is $d_i = z_i - z_{i-1}$.

In any layer $i$, the electric field for a specific polarization (TE or TM) is decomposed into a forward-propagating wave $a_i(z)$ (traveling in the $+z$ direction) and a backward-propagating wave $b_i(z)$ (traveling in the $-z$ direction):

$$
E_i(z) = a_i(z) + b_i(z) = A_i^+ e^{i k_{z,i} (z - z_{i-1})} + A_i^- e^{-i k_{z,i} (z - z_{i-1})}
$$

where $k_{z,i}$ is the normal ($z$-directed) component of the wavevector in layer $i$, and we assume a time dependence of $e^{-i\omega t}$.

## 2. The Concept of the Transfer Matrix

The Transfer Matrix ($M$) maps the total field amplitudes $(a, b)$ from the incident side of a structure to the transmission side (or vice versa). By convention, we typically define the matrix to map the fields from the left side (Ambient, Layer 0) to the right side (Substrate, Layer $N+1$):

$$
\begin{pmatrix} a_0 \\ b_0 \end{pmatrix} = \begin{pmatrix} M_{11} & M_{12} \\ M_{21} & M_{22} \end{pmatrix} \begin{pmatrix} a_{N+1} \\ b_{N+1} \end{pmatrix}
$$

To find the total Transfer Matrix of the entire multilayer stack, we define individual transfer matrices for the interfaces and the propagating layers, and simply multiply them together sequentially.

## 3. The Interface Matrix

Consider an infinitely thin interface between Medium $i$ and Medium $j$. We want to relate the fields immediately to the left of the interface $(a_i, b_i)$ to the fields immediately to the right of the interface $(a_j, b_j)$.

From the S-matrix definitions or fundamental boundary conditions, the transmitted and reflected waves at the interface are governed by the Fresnel coefficients $t_{ij}$ and $r_{ij}$:

$$
a_j = t_{ij} a_i + r_{ji} b_j \quad (1)
$$
$$
b_i = r_{ij} a_i + t_{ji} b_j \quad (2)
$$

To construct a Transfer Matrix, we must express $(a_i, b_i)$ purely in terms of $(a_j, b_j)$.

From equation (1), knowing that $r_{ji} = -r_{ij}$ (Stokes' relations):

$$
t_{ij} a_i = a_j - r_{ji} b_j = a_j + r_{ij} b_j \implies a_i = \frac{1}{t_{ij}}(a_j + r_{ij} b_j)
$$

Substitute $a_i$ into equation (2), using the Stokes' relation $t_{ij} t_{ji} + r_{ij}^2 = 1$:

$$
b_i = r_{ij} \left[ \frac{1}{t_{ij}}(a_j + r_{ij} b_j) \right] + t_{ji} b_j = \frac{r_{ij}}{t_{ij}} a_j + \frac{r_{ij}^2 + t_{ij} t_{ji}}{t_{ij}} b_j = \frac{r_{ij}}{t_{ij}} a_j + \frac{1}{t_{ij}} b_j
$$

Writing these two derived equations in matrix form yields the Interface Transfer Matrix $W_{ij}$:

$$
\begin{pmatrix} a_i \\ b_i \end{pmatrix} = \underbrace{\frac{1}{t_{ij}} \begin{pmatrix} 1 & r_{ij} \\ r_{ij} & 1 \end{pmatrix}}_{W_{ij}} \begin{pmatrix} a_j \\ b_j \end{pmatrix}
$$

For non-magnetic materials, the Fresnel coefficients are:

*   **TE (s-polarization):** $r_{ij} = \frac{k_{z,i} - k_{z,j}}{k_{z,i} + k_{z,j}}$, \quad $t_{ij} = \frac{2k_{z,i}}{k_{z,i} + k_{z,j}}$
*   **TM (p-polarization):** $r_{ij} = \frac{\varepsilon_j k_{z,i} - \varepsilon_i k_{z,j}}{\varepsilon_j k_{z,i} + \varepsilon_i k_{z,j}}$, \quad $t_{ij} = \frac{2\varepsilon_j k_{z,i}}{\varepsilon_j k_{z,i} + \varepsilon_i k_{z,j}} \sqrt{\frac{\varepsilon_i}{\varepsilon_j}}$

## 4. The Propagation Matrix

Next, consider the propagation of waves strictly within layer $i$ across its thickness $d_i$. We need to relate the fields at the left boundary of the layer, $z_{i-1}$, to the fields at the right boundary, $z_i$.

The forward wave $a_i$ travels from left to right. Thus, the field at the left boundary is simply the field at the right boundary "de-propagated" backwards by a phase shift:

$$
a_i(z_{i-1}) = a_i(z_i) e^{-i k_{z,i} d_i}
$$

The backward wave $b_i$ travels from right to left. The field at the left boundary is the field at the right boundary advanced by a phase shift:

$$
b_i(z_{i-1}) = b_i(z_i) e^{i k_{z,i} d_i}
$$

Let the phase thickness of the layer be $\phi_i = k_{z,i} d_i$. Writing the equations in matrix form yields the Propagation Transfer Matrix $P_i$:

$$
\begin{pmatrix} a_i(z_{i-1}) \\ b_i(z_{i-1}) \end{pmatrix} = \underbrace{\begin{pmatrix} e^{-i \phi_i} & 0 \\ 0 & e^{i \phi_i} \end{pmatrix}}_{P_i} \begin{pmatrix} a_i(z_i) \\ b_i(z_i) \end{pmatrix}
$$

*(Note on TMM instability: If $k_{z,i}$ is complex, $e^{i \phi_i}$ grows exponentially. For very thick, lossy layers, this term causes standard floating-point arithmetic to overflow, which is why the Scattering Matrix method is sometimes preferred. For most thin films, however, TMM is perfectly stable and extremely fast).*

## 5. Deriving the Total Multilayer Transfer Matrix

To find the optical response of the entire $N$-layer stack, we sequentially multiply the matrices from the ambient medium ($0$) down to the substrate ($N+1$).

The total field in the ambient medium $\vec{v}_0$ is related to the total field in the substrate $\vec{v}_{N+1}$ via the cascade of interface and propagation matrices:

$$
\begin{pmatrix} a_0 \\ b_0 \end{pmatrix} = W_{01} \cdot P_1 \cdot W_{12} \cdot P_2 \cdot \dots \cdot P_N \cdot W_{N,N+1} \begin{pmatrix} a_{N+1} \\ b_{N+1} \end{pmatrix}
$$

Let the product of all these matrices be the Total Transfer Matrix $M$:

$$
M = \prod_{i=1}^{N} \left( W_{i-1, i} \cdot P_i \right) \cdot W_{N, N+1} = \begin{pmatrix} M_{11} & M_{12} \\ M_{21} & M_{22} \end{pmatrix}
$$

Thus, the global relationship is:

$$
\begin{pmatrix} a_0 \\ b_0 \end{pmatrix} = \begin{pmatrix} M_{11} & M_{12} \\ M_{21} & M_{22} \end{pmatrix} \begin{pmatrix} a_{N+1} \\ b_{N+1} \end{pmatrix}
$$

## 6. Extracting Reflectance and Transmittance (The Single Incident Beam)

To extract physical reflection and transmission coefficients from the matrix $M$, we must strictly define our experimental boundary conditions.

In a standard measurement, light is incident from one side only (typically the ambient Layer $0$). There is no light source hidden inside the substrate firing backwards into the multilayer.

This leads to the following crucial boundary conditions:

1.  **The Incident Beam:** Light enters from Layer $0$. We normalize the amplitude of this incoming forward wave to $1$. Therefore, $a_0 = 1$.
2.  **The Transmitted Beam:** Light exits into Layer $N+1$ and continues traveling forward to infinity. This amplitude is the total transmission coefficient $t$. Therefore, $a_{N+1} = t$.
3.  **The Reflected Beam:** Light bounces off the structure and travels backward into Layer $0$. This amplitude is the total reflection coefficient $r$. Therefore, $b_0 = r$.
4.  **No Backward Incident Beam:** There is no light entering from the substrate side traveling backward. Therefore, $b_{N+1} = 0$.

Substituting these boundary conditions into our global TMM equation gives:

$$
\begin{pmatrix} 1 \\ r \end{pmatrix} = \begin{pmatrix} M_{11} & M_{12} \\ M_{21} & M_{22} \end{pmatrix} \begin{pmatrix} t \\ 0 \end{pmatrix}
$$

This matrix equation yields two simple linear equations:

1.  $1 = M_{11} t + M_{12} (0) \implies \mathbf{1 = M_{11} t}$
2.  $r = M_{21} t + M_{22} (0) \implies \mathbf{r = M_{21} t}$

From Equation 1, we can solve for the Total Amplitude Transmission ($t$):

$$
t = \frac{1}{M_{11}}
$$

Substituting $t$ into Equation 2, we solve for the Total Amplitude Reflection ($r$):

$$
r = \frac{M_{21}}{M_{11}}
$$

Finally, the measurable macroscopic quantities—Reflectance ($R$) and Transmittance ($T$)—are calculated from these complex amplitudes:

$$
R = |r|^2 = \left| \frac{M_{21}}{M_{11}} \right|^2
$$

$$
T = \frac{\text{Re}(k_{z, N+1})}{\text{Re}(k_{z, 0})} |t|^2 = \frac{\text{Re}(k_{z, N+1})}{\text{Re}(k_{z, 0})} \left| \frac{1}{M_{11}} \right|^2
$$

By conservation of energy, the Absorptance ($A$) of the multilayer is:

$$
A = 1 - R - T
$$

This concludes the Transfer Matrix Method derivation, demonstrating how a simple sequence of $2 \times 2$ matrix multiplications naturally yields the global optical properties of any layered structure.
