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

# Apply

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)
::::::

::::::{versionadded} v2025.1.0 After workshop 1
Solutions workshop 1 in text and downloads 
::::::

In this notebook you will work on a homework assignment involving a Vierendeel frame.

+++

Our matrix method implementation is now completely stored in a local package, consisting of three classes.

+++

```{custom_download_link} ./Workshop_1_Apply_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./Workshop_1_Apply_stripped_sol.ipynb
:text: ".ipynb solution"
:replace_default: "False"
```

```{custom_download_link} ./Workshop_1_Apply.md
:text: ".md:myst"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments/tree/solution_workshop_1
:text: "All files practice assignments solutions workshop 1"
:replace_default: "False"
```

+++

## Vierendeel frame

+++

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/vierendeel.png
:align: center
:width: 400
```

With:

- $h = 1$
- $b = 1$
- $EI_r = 10000$
- $EI_k = 1000$
- $EA  = 1\cdot 10^{10}$
- $H = 100$

+++

In the first half of this course last quarter, you have learned that the deformation of Vierendeel frames (an example of which is shown above) can be obtained in a simplified way by assuming the global deformation can be described by a shear beam with equivalent stiffness given by:

$$
k = \frac{24}{h\left(\displaystyle\frac{h}{EI_k}+\frac{b}{EI_r}\right)}
$$

```{exercise-start} Workshop 1 - Apply
:label: exercise_ws_1
:nonumber: true

Now that you have the tools to solve the original frame problem using the Matrix Method, your task in this assignment is to investigate the validity of this equivalent shear beam model.

Note that the checks only had a single element. For this model you need to obtain $\mathbf{K}$ and $\mathbf{f}$ of all elements and add them to the correct locations in the global stiffness matrix and force vector. To do that, make use of the `global_dofs` function of the Element class and the `np.ix_` Numpy utility function. (Tip: refer back to what you did in the `constrain` function).

Once you have a solution, use SymPy / Maple / pen and paper to solve a shear beam problem with the equivalent stiffness given above (It is very similar to the simple extension problem above) and compare the horizontal displacement at the point of application of $H$ for the two models.

Investigate how the two models compare for different values of $EA$, ranging from very small (*e.g.* $1\cdot 10^{-5}$) to very large (*e.g.* $1\cdot10^{10}$). What explains the behavior you observe?
```

```{code-cell} ipython3
:tags: [thebe-remove-input-init]

import matplotlib as plt
import numpy as np
sys.path.insert(1, '/matrixmethod_solution')
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [remove-cell]

import matplotlib as plt
import numpy as np
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [disable-execution-cell]

import numpy as np
import matplotlib as plt
import matrixmethod as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
mm.Node.clear()
mm.Element.clear()
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{code-cell} ipython3
global_k #= np.zeros(YOUR CODE HERE
global_f #= np.zeros(YOUR CODE HERE

for elem in elems:
    elmat = elem#.YOUR CODE HERE
    idofs = elem#.YOUR CODE HERE
    
    #YOUR CODE HERE

for node in nodes:
    #YOUR CODE HERE
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{code-cell} ipython3
#provided in case you want to solve the shear beam problem using SymPy
import sympy as sym
x, k, L, H = sym.symbols('x, k, L, H')
w = sym.Function('w')

ODE_shear = #YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise_ws_1
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()
```

```{code-cell} ipython3
:tags: [thebe-init]

h = 1
b = 1
EIr = 10000
EIk = 1000
EA  = 1e10
H = 100

nodes = []

nodes.append(mm.Node(0,0))
nodes.append(mm.Node(b,0))
nodes.append(mm.Node(b,-h))
nodes.append(mm.Node(0,-h))

elems = []

elems.append(mm.Element(nodes[0], nodes[1]))
elems.append(mm.Element(nodes[1], nodes[2]))
elems.append(mm.Element(nodes[2], nodes[3]))
elems.append(mm.Element(nodes[0], nodes[3]))

beams = {}
columns = {}
beams['EI'] = EIr
beams['EA'] = EA
columns['EI'] = EIk
columns['EA'] = EA

elems[0].set_section (beams)
elems[1].set_section (columns)
elems[2].set_section (beams)
elems[3].set_section (columns)

for elem in elems:
    print(elem)

con = mm.Constrainer()

con.fix_dof (nodes[0], 0)
con.fix_dof (nodes[0], 1)
con.fix_dof (nodes[1], 1)

nodes[3].add_load ([H,0,0])
```

```{code-cell} ipython3
:tags: [thebe-init]

global_k = np.zeros ((3*len(nodes), 3*len(nodes)))
global_f = np.zeros (3*len(nodes))

for elem in elems:
    elmat = elem.stiffness()
    idofs = elem.global_dofs()
    
    global_k[np.ix_(idofs,idofs)] += elmat

for node in nodes:
    global_f[node.dofs] += node.p
```

```{code-cell} ipython3
:tags: [thebe-init]

Kff, Ff = con.constrain ( global_k, global_f )
u = np.matmul ( np.linalg.inv(Kff), Ff )
print(u)
```

```{code-cell} ipython3
:tags: [thebe-init]

#provided in case you want to solve the shear beam problem using SymPy
import sympy as sym

x, k, L, H = sym.symbols('x, k, L, H')
w = sym.Function('w')

ODE_shear = sym.Eq(w(x).diff(x, 2) *k, 0)
w = sym.dsolve(ODE_shear, w(x)).rhs

gamma = w.diff(x)
V = k * gamma
eq1 = sym.Eq(w.subs(x,0),0)
eq2 = sym.Eq(V.subs(x,L),H)
C_sol = sym.solve([eq1, eq2], sym.symbols('C1, C2'))
display(w.subs(C_sol))
```

As derived above, the displacement of a shear beam equals $w = \cfrac{Hx}{k}$. With $k = \cfrac{24}{h\left(\cfrac{h}{EI_k}+\cfrac{b}{EI_r}\right)} = \cfrac{24}{h\left(\cfrac{1}{1000}+\cfrac{1}{10000}\right)} \approx 21818$ this gives $w = \cfrac{100\cdot 1}{28181} \approx 0.0045833$ which is equal to `u[6]`, corresponding to the horizontal displacement of the top left node.

```{solution-end}
```
