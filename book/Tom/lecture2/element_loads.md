# Element loads

As the matrix method is a discrete approach, nodal loads were treated with ease. However, what to do with continuous loads or loads which are not applied at the nodes?

::::::{topic} Learning objective
You'll look into how to model element loads using differential equations and conservation of work and how these are combined in the matrix formulation.
::::::

## Force-displacement relations using differential equations
As [before](../lecture1/single_element.md), we can derive the force-displacement relations of a single extension element. However, now let's include the loads, for example a continuous load $q$

```{figure} extensionelementq.svg
:name: extensionelementq
:align: center

Single extension element with distributed load $q$
```

The same approach is used as in [](../lecture1/recap.ipynb). This results in:

- $C_1 = \cB{\cfrac{q\ell}{2EA}}+\cfrac{u_2-u_1}{\ell}$
- $C_2 = u_1$

The continuous distributions for the displacement and section force can be evaluated too:
- $u(x) = \cB{\cfrac{q}{2EA}\left(\ell x - x^2\right)} + u_1\left(1-\cfrac{x}{\ell}\right) + u_2\cfrac{x}{\ell}$
- $N(x) = \cB{\cfrac{q}{2}\left(\ell-2x \right)} -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell}u_2$

These are extended results in comparison to [before](../lecture1/single_element.md)

## Combine elements
As before, we can glue elements together by applying nodal equilibrium:

This leads to:
- $ F_1 = - N = \cfrac{EA}{\ell}u_1-\cfrac{EA}{\ell}u_2 \cB{ -\cfrac{q\ell}{2}}$
- $ F_2 = N = -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell} u_2 \cB{-\cfrac{q\ell}{2}}$

or in matrix notation:

$$\cfrac{EA}{\ell}\begin{bmatrix} 1&-1\\-1&1 \end{bmatrix}\begin{bmatrix}u_1\\u_2\end{bmatrix} \cB{- \begin{bmatrix} \cfrac{q\ell}{2}\\ \cfrac{q\ell}{2}\end{bmatrix}}= \begin{bmatrix}F_1\\F_2\end{bmatrix}$$

Effectively, we converted the continuous load to an equivalent nodal load.

This influences the nodal equilibrium (in the global coordinate system) too:

$$\begin{align} -\sum_e\mathbf{f}^e + \mathbf{f}_\text{nodal}& = \mathbf{0}\\
-\sum_e\left(\mathbf{K}^e\mathbf{u}^e \cA{-\mathbf{f}_\mathrm{eq}^e}\right) + \mathbf{f}_\mathrm{nodal}& = \mathbf{0}\\
\sum_e\mathbf{f}^e& = \mathbf{f}_\text{nodal} + \cA{\sum_e\mathbf{f}_\mathrm{eq}^e} \end{align}$$

This means that we can calculate all the equivalent nodal loads separately and add them to the nodal loads. Please note that all these forces should be in the global coordinate system: $\cA{\mathbf{f}_\mathrm{eq}} = \mathbf{T}^\mathrm{T}\cB{\bar{\mathbf{f}}_\mathrm{eq}}$

## Force-displacement relations using conservation of work

### Point load
Let's consider another example in which the application of the differential equations are not trivial: we'll introduce a discontinuous force in the form of a point load halfway the element. We'll compare this with the same element loaded by vertical forces and bending moments at the ends:

```{figure} eqwork3.svg
:name: eqwork3
:align: center

Single extension element with point load $P$ halfway in comparison to element loaded by vertical forces and bending moments at the ends 
```

The conservation of work states that the work done by the force $P$ should be equal to the force done by the forces $\cA{F_1^{eq}}$, $\cE{F_2^{eq}}$, $\cB{T_1^{eq}}$ and $\cI{T_2^{eq}}$: $W_P$ = $W_{eq}$

To ease the calculation, we'll split the displacement in four cubic shape functions (cubic is justified because with $q=0$ a cubic displacement function is found). Each shape function will have a displacement of $1$ in direction of each of the four edge forces ($\cA{w_1}$, $\cE{\varphi_1}$, $\cB{w_2}$ and $\cI{\varphi_2}$):

```{figure} shape_functions.svg
:name: shape_functions
:align: center

Four shape functions
```

The shape functions have the following function:

$$
w(x) = 
	    \cA{\underbrace{\left(\cfrac{2x^3}{\ell^3}-\cfrac{3x^2}{\ell^2}+1\right)}_{s_1}w_1} +
	    \cE{\underbrace{\left(-\cfrac{x^3}{\ell^2}+\cfrac{2x^2}{\ell}-x\right)}_{s_2}\varphi_1} +
	    \cB{\underbrace{\left(-\cfrac{2x^3}{\ell^3}+\cfrac{3x^2}{\ell^2}\right)}_{s_3}w_2} +
	    \cI{\underbrace{\left(-\cfrac{x^3}{\ell^2}+\cfrac{x^2}{\ell}\right)}_{s_4}\varphi_2}
$$

Consequently, the work $W_{eq}$ can be splitted in four independent terms easing the calculation, which can be compared to the same terms in $W_P$.

The work performed by the edge forces equals:

$$ W_{eq} = 
	    \cA{F^\mrm{eq}_1w_1} +
	    \cE{T^\mrm{eq}_1\varphi_1} +
	    \cB{F^\mrm{eq}_2w_2} + 
	    \cI{T^\mrm{eq}_2\varphi_2}$$

The work performed by $P$ (under the same displacement) is:

$$W_P = P\ w\left(\frac{\ell}{2}\right) = 
	  \cA{P\ s_1\left(\frac{\ell}{2}\right) \ w_1} +
	  \cE{P\ s_2\left(\frac{\ell}{2}\right) \ \varphi_1} +
	  \cB{P\ s_3\left(\frac{\ell}{2}\right) \ w_2} +
	  \cI{P\ s_4\left(\frac{\ell}{2}\right) \ \varphi_2}$$

Enforcing $W_F = W_q$ and isolating terms gives:

$$
\mathbf{f}_\mrm{eq}
	    =
	    \begin{bmatrix}
	      \cA{F^\mrm{eq}_1}\\[10pt]
	      \cE{T^\mrm{eq}_1}\\[10pt]
	      \cB{F^\mrm{eq}_2}\\[10pt]
	      \cI{T^\mrm{eq}_2}\\
	    \end{bmatrix}
	    =
	    \begin{bmatrix}
	      \cA{\cfrac{P}{2}}     \\
	      \cE{-\cfrac{P\ell}{8}}\\
	      \cB{\cfrac{P}{2}}     \\
	      \cI{\cfrac{P\ell}{8}} \\
	    \end{bmatrix}
$$

In the local coordinate system.

### Distributed load
The same can be done for a distributed load:

$$W_q = \int_\ell q \ w(x)\mrm{d}x = 
	  \cA{\int_\ell q \ s_1(x)\mrm{d}x \ w_1} +
	  \cE{\int_\ell q \ s_2(x)\mrm{d}x \ \varphi_1} +
	  \cB{\int_\ell q \ s_3(x)\mrm{d}x \ w_2} +
	  \cI{\int_\ell q \ s_4(x)\mrm{d}x \ \varphi_2}$$

Leading to:

$$
\mathbf{f}_\mrm{eq}
	    =
	    \begin{bmatrix}
	      \cA{F^\mrm{eq}_1}\\[10pt]
	      \cE{T^\mrm{eq}_1}\\[10pt]
	      \cB{F^\mrm{eq}_2}\\[10pt]
	      \cI{T^\mrm{eq}_2}\\
	    \end{bmatrix}
	    =
	    \begin{bmatrix}
	      \cA{\cfrac{q\ell}{2}}     \\
	      \cE{-\cfrac{q\ell^2}{12}} \\
	      \cB{\cfrac{q\ell}{2}}   \\
	      \cI{\cfrac{q\ell^2}{12}} \\
	    \end{bmatrix}
$$

## Example

Let us use what we have just learned on a simple example:

```{figure} twobarqdiri.svg
:name: twobarqdiri
:align: center
```

This example has the same two-element bar model as in [](../lecture1/directly.md), so the unconstrained stiffness matrix is unchanged:

$$
\begin{bmatrix}
	\cA{\cfrac{EA_1}{\ell_1}} & \cA{-\cfrac{EA_1}{\ell_1}} & 0\\
	\cA{-\cfrac{EA_1}{\ell_1}} & \cA{\cfrac{EA_1}{\ell_1}} + \cB{\cfrac{EA_2}{\ell_2}} & \cB{-\cfrac{EA_2}{\ell_2}}\\
	0 & \cB{-\cfrac{EA_2}{\ell_2}} & \cB{\cfrac{EA_2}{\ell_2}}\\
	  \end{bmatrix}
$$

Both supports will introduce a support reaction in the form of a Neumann boundary condition:

```{figure} assembly5_3.svg
:name: assembly5_3
:align: center
```

$$
	  \begin{bmatrix}
	\cA{\cfrac{EA_1}{\ell_1}} & \cA{-\cfrac{EA_1}{\ell_1}} & 0\\
	\cA{-\cfrac{EA_1}{\ell_1}} & \cA{\cfrac{EA_1}{\ell_1}} + \cB{\cfrac{EA_2}{\ell_2}} & \cB{-\cfrac{EA_2}{\ell_2}}\\
	0 & \cB{-\cfrac{EA_2}{\ell_2}} & \cB{\cfrac{EA_2}{\ell_2}}\\
	  \end{bmatrix}
	  \begin{bmatrix}
	u_1\\[12pt]u_2\\[12pt]u_3
	  \end{bmatrix}
	  =
	  \begin{bmatrix}
	H_1\\[12pt]0\\[12pt]H_3
	  \end{bmatrix}
$$

The distributed load can be added directly as equivalent load vector ${\mathbf{f}_\mathrm{eq}}$ to the global force vector as the local coordinate system aligns with the global coordinate system. Doing so for both elements leads to:

$$
	  \begin{bmatrix}
	\cA{\cfrac{EA_1}{\ell_1}} & \cA{-\cfrac{EA_1}{\ell_1}} & 0\\
	\cA{-\cfrac{EA_1}{\ell_1}} & \cA{\cfrac{EA_1}{\ell_1}} + \cB{\cfrac{EA_2}{\ell_2}} & \cB{-\cfrac{EA_2}{\ell_2}}\\
	0 & \cB{-\cfrac{EA_2}{\ell_2}} & \cB{\cfrac{EA_2}{\ell_2}}\\
	  \end{bmatrix}
	  \begin{bmatrix}
	u_1\\[12pt]u_2\\[12pt]u_3
	  \end{bmatrix}
	  =
	  \begin{bmatrix}
	H_1 \cA{ + \cfrac{q\ell_1}{2}}\\ \cA{\cfrac{q\ell_1}{2}} \cB{+ \cfrac{q\ell_2}{2}}\\H_3 \cB{ + \cfrac{q\ell_2}{2}}
	  \end{bmatrix}
$$