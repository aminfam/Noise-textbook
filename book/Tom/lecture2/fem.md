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

# Finite element method vs. Matrix Method

> This page reuses <license> content from {cite:t}`MUDE`. {fa}`quote-left`{ref}`Find out more here.<external_resources>`

You've seen the finite element method before, which could be used to solve similar problems. But what are the differences?

::::::{topic} Learning objective
We'll investigate the differences and equivalence between solving structures with the finite element method and the matrix method.
::::::

+++

## Equivalence with matrix method

Although the two method can give the same results, the methods are different.

The matrix method solves the strong from of the differential equation, as derived in [](../lecture1/single_element.md). The finite element method solves the weak form by multiplying the strong form by a test function {cite:p}`MUDE2`. In doing so, the choice for the shape function of the test-functions and approximate solution matters. The two methods end up with the same solution if the "approximation" assumed by FEM (linear shape functions for extension, cubic for bending) turn out to be the exact ODE solution.

In terms of global and local coordinate systems, there's an additional difference. Where the matrix method solves the nodal displacements and support reactions globally, using locally derived force-displacement relations. On the contrary, the finite element method solves the weak form globally with shape functions defined globally.

Finally, the matrix method has limitations. It turns out to be impossible to glue element through equilibrium for twodimensional elements. Furthermore, exact solution for the differential equations, required for defining the local stiffness matrix, generally do not exist for twodimensional elements.

:::{table} Equivalence finite element- and matrix method
:widths: auto
:align: center

|Finite element method|Matrix method|
|:-:|:-:|
|Solves weak form|Solves the strong form|
|Solves everything globally|Solves discrete solution globally, with locally derived equations|
|Generally applicable to differential equations| Only applicable for 1D elements|

:::

+++

## Example with finite element method

Let's consider the examples from [](../lecture1/displacement.md):

```{figure} ../lecture1/extension2fieldsdisp.svg
:name: extension2fieldsq2
:align: center

Statically indeterminate extension bar
```

We'll apply the finite element method by using the matrix implementation from [MUDE](https://mude.citg.tudelft.nl/2024/book/fem/matrix.html) {cite:p}`MUDE`. We use linear shape functions: $ N_a(x) = \cfrac{x_b-x}{x_b-x_a}=\cfrac{x_b-x}{\Delta x}$ and $N_b(x)=\cfrac{x-x_a}{x_b-x_a}=\cfrac{x-x_a}{\Delta x}$. For this specific example, this matches the linear normal force distribution, leading to similar results:

```{code-cell} ipython3
:tags: [hide-cell, thebe-init]

import numpy as np
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
:tags: [hide-cell, thebe-init]

# Returns the evaluated N matrix
# - The local coordinate of evaluation "x_local"
# - The element size "dx"
def evaluate_N(x_local, dx):
    return np.array([[1-x_local/dx, x_local/dx]])

# Returns the evaluated B vector
# - The local coordinate of evaluation "x_local"
# - The element size "dx"
def evaluate_B(x_local, dx):
    return np.array([[-1/dx, 1/dx]])

# Returning the matrix of a single element
# - The stiffness of the rod "EA"
# - The length of an element "dx"
def get_element_matrix(EA, dx):
    
    # Defining integration locations and weights
    integration_locations = [(dx - dx/(3**0.5))/2, (dx + dx/(3**0.5))/2]
    integration_weights = [dx/2, dx/2]
    n_ip = len(integration_weights)

    # Setting up the local element matrix
    n_node = 2
    K_loc = np.zeros((n_node,n_node))

    # Evaluation of the matrix in a loop over integration points
    for x_i, w_i in zip(integration_locations, integration_weights):
        B = evaluate_B(x_i, dx)
        K_loc += EA*np.dot(np.transpose(B), B)*w_i

    return K_loc

# Returning the force vector of a single element (distrubted load only)
# - The magnitude of constant distributed load "q"
# - The length of the element "dx"
def get_element_force(q, dx):
    # Creating a matrix with the required size
    n_node = 2
    N = np.zeros((1,n_node))    # Defining in 2 dimensions is nessecary for transpose
    
    # Defining integration locations and weights
    integration_locations = [(dx - dx/(3**0.5))/2, (dx + dx/(3**0.5))/2]
    integration_weights = [dx/2, dx/2]
    
    # Setting up the local element force vector
    f_loc = np.zeros((n_node,1))
    
    for x_i, w_i in zip(integration_locations, integration_weights):
        N = evaluate_N(x_i,dx)
        f_loc += np.transpose(N)*q*w_i

    return f_loc

def get_nodes_for_element(ie):
    return np.array([ie,ie+1])

# - The length of the rod "rod_length"
# - The number of elements "n_el"
# - The stiffness of the rod "EA"
def assemble_global_K(rod_length, n_el, EA):
    n_DOF = n_el+1
    dx = rod_length/n_el
    K_global = np.zeros((n_DOF, n_DOF))
    
    for i in range(n_el):
        elnodes = get_nodes_for_element(i)
        K_global[np.ix_(elnodes,elnodes)] += get_element_matrix(EA, dx)
    
    return K_global

# Returns the global f vector with continuous forces on a rod
# - The length of the rod "rod_length"
# - The number of elements "n_el"
# - The force on the rod "EA"
def assemble_global_f(rod_length, n_el, q):
    n_DOF = n_el+1
    dx = rod_length/n_el
    f_global = np.zeros((n_DOF,1))
    
    for i in range(n_el):
        elnodes = get_nodes_for_element(i) 
        f_global[elnodes] += get_element_force(q, dx)
        
    return np.squeeze(f_global)
```

```{code-cell} ipython3
:tags: [hide-input, thebe-init]

def simulate(n_element):
    length = 3              # Length
    EA = 1e3                # Stiffness EA
    n_node = n_element + 1  # Number of nodes
    u_left = 0              # Displacement at the left boundary
    u_right = 0             # Displacement at the right boundary
    q_load = 10             # Distibuted load

    dx = length/n_element
    x = np.linspace(0,length,n_node)

    # Assmemble K for the rod problem
    K = assemble_global_K(length, n_element, EA)

    # Assemble f first with distributed load only
    f = assemble_global_f(length, n_element, q_load)

    # Initialize vector u with displacements
    u = np.zeros(n_node)

    # Solve the system for unknown displacements
    K_inv = np.linalg.inv(K[1:n_node-1, 1:n_node-1])
    u[1:n_node-1] = np.dot(K_inv, f[1:n_node-1])

    return x, u
```

```{code-cell} ipython3
:tags: [thebe-init]

x, u = simulate(2)
print(u[1])
```

## Example with matrix method

Applying the matrix method (including implementations you'll implement during [](../workshop2.md)) yields the same result:

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
:tags: [disable-execution-cell, hide-cell]

import numpy as np
import matplotlib as plt
import matrixmethod as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [hide-input, thebe-init]

mm.Node.clear()
mm.Element.clear()

nodes = []

nodes.append(mm.Node(0,0))
nodes.append(mm.Node(1.5,0))
nodes.append(mm.Node(3,0))

elems = []

elems.append(mm.Element(nodes[0], nodes[1]))
elems.append(mm.Element(nodes[1], nodes[2]))

for elem in elems:
    elem.set_section({'EA': 1e3})
    elem.add_distributed_load([10,0])

con = mm.Constrainer()

con.fix_node(nodes[0])
con.fix_node(nodes[2])

global_k = np.zeros ((3*len(nodes), 3*len(nodes)))
global_f = np.zeros (3*len(nodes))

for elem in elems:
    elmat = elem.stiffness()
    idofs = elem.global_dofs()
    
    global_k[np.ix_(idofs,idofs)] += elmat

for node in nodes:
    global_f[node.dofs] += node.p

Kff, Ff = con.constrain ( global_k, global_f )
u_free = np.matmul ( np.linalg.inv(Kff), Ff )
```

```{code-cell} ipython3
:tags: [thebe-init]

print(u_free[0])
```

Which is the same result (not true in general).
