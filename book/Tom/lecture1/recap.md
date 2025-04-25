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

# Recap differential equations for structures

```{custom_download_link} recap.md
:text: ".md:myst"
:replace_default: "False"
```

Let's recap how to solve a structure using differential equations.

::::::{topic} Learning objective
We'll investigate the differences between solving structures with differential equations and the matrix method. Furthermore, you'll need differential equations in the matrix method itself to define some default elements. 
::::::

Given is a cantilever beam with a distributed load:

```{figure} bending1field.svg
:name: bending1field
:align: center

Cantilever beam with distributed load
```

The differential equation for the Euler-Bernoulli model can be derived leading to:
- Kinematic relations:
  - $\varphi=-\cfrac{\mrm{d}w}{\mrm{d}x}$
  - $\kappa=\cfrac{\mrm{d}\varphi}{\mrm{d}x}$
- Constitutive relation:
  - $M=EI\kappa$
- Equilibrium relations:
  - $\cfrac{\mrm{d}V}{\mrm{d}x}=-q$
  - $\cfrac{\mrm{d}M}{\mrm{d}x}=V$

These relations can be combined into one fourth order differential equation:

$$ EI\frac{\mrm{d}^4w}{\mrm{d}x^4}=q $$

+++

This differential equation can be solved directly to solve structures.

+++

## Solving differential equation of one field by hand

This differential equations can be solved by integrating:

- $V(x) = -qx + \bar C_{1}$
- $M(x) = -\cfrac{qx^2}{2} + \bar C_1 x + \bar C_2$
- $\kappa(x) = -\cfrac{qx^2}{2EI}+ \tilde C_1 x + \tilde C_2$
- $\varphi(x) = -\cfrac{qx^3}{6EI} + \cfrac{\tilde C_1x^2}{2} + \tilde C_2 x+ \tilde C_3$
- $w(x) = \cfrac{qx^4}{24EI} + \cfrac{C_1x^3}{6} + \cfrac{C_2x^2}{2}+C_3x+C_4 $

The boundary conditions follow from the clamped side at $x=0$ and free end at $x=\ell$:

- $w(0) = 0$
- $\varphi(0) = 0$
- $	M(\ell) = 0$
- $	V(\ell) = 0$

Solving these four equations for the integration constants gives:

- $C_1 = -\cfrac{q\ell}{EI}$
- $C_2 = \cfrac{q\ell^2}{2EI}$
- $C_3 = 0$
- $C_4 = 0$

Substituting these constants, a final solution for $w$ can be found:

$w(x) = \cfrac{qx^4}{24EI} - \cfrac{q\ell x^3}{6EI} + \cfrac{q\ell^2x^2}{4EI}$

+++

(sympy-ode)=
## Solving differential equation of one field using SymPy

The differential equation can also be solved using SymPy:

```{code-cell} ipython3
:tags: [thebe-remove-input-init]

import sympy as sym
sym.init_printing()
```

```{code-cell} ipython3
import sympy as sym
```

```{code-cell} ipython3
EI, q, x, L = sym.symbols('EI, q, x, ell')
C_1, C_2, C_3, C_4 = sym.symbols('C_1, C_2, C_3, C_4')
```

```{code-cell} ipython3
V = -sym.integrate(q,x) + C_1
M = sym.integrate(V,x) + C_2
kappa = M/EI
phi = sym.integrate(kappa,x) + C_3
w = -sym.integrate(phi,x) + C_4
display(w)
```

```{code-cell} ipython3
eq1 = sym.Eq(w.subs(x,0),0)
eq2 = sym.Eq(phi.subs(x,0),0)
eq3 = sym.Eq(M.subs(x,L),0)
eq4 = sym.Eq(V.subs(x,L),0)
C_sol = sym.solve([eq1, eq2, eq3, eq4 ], [C_1, C_2, C_3, C_4])
for key in C_sol:
    display(sym.Eq(key, C_sol[key]))
```

```{code-cell} ipython3
w.subs(C_sol)
```

## Solving differential equations of two fields

A similar approach can be taken when solving two fields.

Let's investigate the following structure, consisting of two fields.

```{figure} extension2fields.svg
:name: extension2fields
:align: center

Extension bar with nodal load
```

As this structure is loaded along its axis, the differential equation for extension is used.

For the first field this gives:

- $\cA{EA_1\cfrac{\mrm{d}^2u_1}{\mrm{d}x^2}=0}$
- $\cA{N_1 = C_1}$
- $\cA{u_1(x) = \cfrac{C_1}{EA}x + C_2}$
- Boundary conditions: $\cA{u_1(0) = 0}$

For the second field it gives:

- $\cB{EA_2\cfrac{\mrm{d}^2u_2}{\mrm{d}x^2}=0}$
- $\cB{N_2 = C_3}$
- $\cB{u_2(x) = \cfrac{C_3}{EA}x + C_4}$
- Boundary conditions: $\cB{N_2(\ell_1+\ell_2) = F}$

The two remaining integration constants can be solved by specifying interface conditions:
- $\cA{u_1(\ell_1)} = \cB{u_2(\ell_1)}$
- $\cA{N_1} = \cB{N_2}$

+++

## Solving differential equations of more fields

The same approach can be taken to tackle problems with more field, like the one below:

```{figure} bigframe.svg
:name: bigframe
:align: center

Frame structure with many fields
```

How many integration constants should be solved for here? How many boundary- and interface conditions would be needed for that? It gets annoying very quickly as each of these conditions need to be defined carefully.

## Equivalence with matrix method

While both methods segment the structure in different parts, the matrix method applies a different principle in solving the structure than when directly solving differential equations: instead of solving for integration constants, nodal displacements are solved for. This shows big potential because setting up all the boundary- and interface conditions can be tedious and is problem-specific. The matrix method applies a generic algorithmic approach to combine all unknown nodal displacements

The similarities and differences are shown in the table below.

:::{table} Equivalence solving differential equations and matrix method
:widths: auto
:align: center

|Solving differential equations|Matrix method|
|:-:|:-:|
|Segment structure in separate fields|Segment structure in mostly repetitive elements|
|Define all boundary- and interface conditions|Define relations in generic algorithmic manner|
|Solve for integration constants $C_1, C_2, ...$|Solve nodal displacements $u_1, u_2, ...$|

:::
