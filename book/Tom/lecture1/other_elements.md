---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: base
  language: python
  name: python3
---

+++ {"vscode": {"languageId": "plaintext"}}

# Local stiffness matrix Euler-Bernoulli element

```{custom_download_link} other_elements.md
:text: ".md:myst"
:replace_default: "False"
```

In [](./single_element.md) and [](./matrix.md) you've seen how to derive the local stiffness matrix for a simple extension element. But can you do the same for other elements?


::::::{topic} Learning objective
You'll look into deriving the local stiffness matrix for the Euler-Bernoulli element.
::::::

Previously, the local stiffness matrix for a simple extension element was found as:

$$\mathbf{K}^{(e)} = \cfrac{EA}{\ell}\begin{bmatrix} 1&-1\\-1&1 \end{bmatrix}$$

The same procedure can be followed for other element, like a combined extension and Euler-Bernoulli element:

```{figure} elemtypes.svg
:name: elemtypes
:align: center

Combined extension and Euler-Bernoulli element
```

The amount of degrees of freedom increases, as both ends of the element can translate in two directions and rotate. However, the approach is exactly the same, leading to the following element stiffness matrix:

$$
	\mathbf{K}^{(e)} = \begin{bmatrix}
	  \cfrac{EA}{\ell} & 0 & 0 & -\cfrac{EA}{\ell} & 0 & 0\\
	  0 & \cfrac{12EI}{\ell^3} & -\cfrac{6EI}{\ell^2} & 0 & -\cfrac{12EI}{\ell^3} & -\cfrac{6EI}{\ell^2}\\
	  0 & -\cfrac{6EI}{\ell^2} & \cfrac{4EI}{\ell} & 0 & \cfrac{6EI}{\ell^2} & \cfrac{2EI}{\ell}\\
	  -\cfrac{EA}{\ell} & 0 & 0 & \cfrac{EA}{\ell} & 0 & 0\\
	  0 & -\cfrac{12EI}{\ell^3} & \cfrac{6EI}{\ell^2} & 0 & \cfrac{12EI}{\ell^3} & \cfrac{6EI}{\ell^2}\\
	  0 & -\cfrac{6EI}{\ell^2} & \cfrac{2EI}{\ell} & 0 & \cfrac{6EI}{\ell^2} & \cfrac{4EI}{\ell}\\
	\end{bmatrix}$$

for

$$
  \mathbf{u}^{(e)} =\begin{bmatrix}
	  u_1\\
	  w_1 \\
	  \varphi_1 \\
	  u_2 \\
	  w_2\\
	  \varphi_2\\
	\end{bmatrix} 
	$$

## Derivation using SymPy
We can make use of software like SymPy, as we did before in [](sympy-ode) to do the calculations in this derivation:

```{code-cell} ipython3
:tags: [thebe-remove-input-init]

import sympy as sym
sym.init_printing()
```

```{code-cell} ipython3
import sympy as sym
```

```{code-cell} ipython3
EI, x, L = sym.symbols('EI, x, L')
w = sym.Function('w')

ODE_bending = sym.Eq(w(x).diff(x, 4) * EI, 0)
display(ODE_bending)
```

```{code-cell} ipython3
w = sym.dsolve(ODE_bending, w(x)).rhs
display(w)
```

```{code-cell} ipython3
phi = -w.diff(x)
kappa = phi.diff(x)
M = EI * kappa
V = M.diff(x)
```

```{code-cell} ipython3
w_1, w_2, phi_1, phi_2 = sym.symbols('w_1, w_2, phi_1, phi_2')

eq1 = sym.Eq(w.subs(x,0),w_1)
eq2 = sym.Eq(w.subs(x,L),w_2)
eq3 = sym.Eq(phi.subs(x,0),phi_1)
eq4 = sym.Eq(phi.subs(x,L),phi_2)

sol = sym.solve([eq1, eq2, eq3, eq4 ], sym.symbols('C1, C2, C3, C4'))
for key in sol:
    display(sym.Eq(key, sol[key]))
```

```{code-cell} ipython3
F_1_z, F_2_z, T_1_y, T_2_y = sym.symbols('F_1_z, F_2_z, T_1_y, T_2_y')

eq5 = sym.Eq(-V.subs(sol).subs(x,0), F_1_z)
eq6 = sym.Eq(V.subs(sol).subs(x,L), F_2_z)
eq7 = sym.Eq(-M.subs(sol).subs(x,0), T_1_y)
eq8 = sym.Eq(M.subs(sol).subs(x,L), T_2_y)
```

```{code-cell} ipython3
K_e, f_e = sym.linear_eq_to_matrix([eq5,eq7, eq6, eq8], [w_1, phi_1, w_2, phi_2])
display(K_e)
```

To use the stiffness matrix without manually copying it over, you can make use of the `lambdify` which converts a symbolic SymPy object in a python function. This allows you to evaluate it for specific numerical values and continue using it in the numerical framework of the matrix method.

```{code-cell} ipython3
K = sym.lambdify((L, EI), K_e)
print(K.__doc__)
```

```{code-cell} ipython3
print('Example of K with L=5 and EI=1000:\n',K(5,1000))
```
