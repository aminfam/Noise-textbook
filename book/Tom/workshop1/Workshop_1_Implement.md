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

# Implementation

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)
::::::

::::::{versionadded} v2025.1.0 After workshop 1
Solutions workshop 1 in text and downloads 
::::::

::::::{versionchanged} v2025.0.3 2025-02-10 13:33, before workshop 1
Fixed typo in [Exercise 2.6](exercise2.6)
::::::

In this notebook you will implement the matrix method and check it with some sanity checks.

+++

Our matrix method implementation is now completely stored in a local package, consisting of three classes. If you need a refresher on how to code with Classes and Objects, refer to the [section on Object Oriented Programming in the MUDE-book](https://mude.citg.tudelft.nl/2024/book/external/learn-programming/book/python/oop/classes.html), with additionally [programming assignment 1.7](https://mude.citg.tudelft.nl/2024/files/Week_1_7/PA_1_7_classy_distributions.html)

+++

```{custom_download_link} ./Workshop_1_Implement_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./Workshop_1_Implement_stripped_sol.ipynb
:text: ".ipynb solution"
:replace_default: "False"
```

```{custom_download_link} ./Workshop_1_Implement.md
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

## 1. The Node class
This class is stored in [`./matrixmethod/node.py`](../matrixmethod/node.md)

The purpose of this class is to store node information and keep track of the total number of DOFs of the problem. Note the automatic bookkeeping we introduce in `__init__`. This simple but efficient way of keeping track of which DOFs belong to which nodes will make life much easier when we need to assemble matrices from multiple elements. The Node class doesn't need any modification.

```{exercise-start} Workshop 1 - 1.1
:label: exercise1.1
:nonumber: true
:class: exercise

To test whether you understand how the class works, create two nodes on coordinates ($0$,$0$) and ($3$,$4$) and print the string representation of both nodes. The `clear` function is called to restart the node and DOF counters. Make sure this is done whenever you start solving a new problem.
```

```{code-cell} ipython3
mm.Node.clear()

node1 #= mm.Node(YOUR CODE HERE

print(node1)
#YOUR CODE HERE
```

```{exercise-end}
```

+++

::::::{hint}
:class: dropdown

Your output should look like this:

```
This node has:
 - x coordinate=0,
 - z coordinate=0,
 - degrees of freedom=[0, 1, 2],
 - load vector=[0. 0. 0.]
This node has:
 - x coordinate=3,
 - z coordinate=4,
 - degrees of freedom=[3, 4, 5],
 - load vector=[0. 0. 0.]
```
::::::

+++

```{solution-start} exercise1.1
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()

node1 = mm.Node(0,0)
node2 = mm.Node(3,4)
print(node1)
print(node2)
```

```{solution-end}
```

+++

## 2. The Element class
This class is stored in [`./matrixmethod/elements.py`](../matrixmethod/elements.md)

This class keeps track of each element in the model, including:
- Cross-section properties
- Element orientation (for coordinate system transformations)
- Which Nodes make up each element, and in turn (with help of the Node class) which DOFs belong to each element

Apart from bookkeeping element data, the other main task of this class is to provide the element stiffness matrix in the global coordinate system (for subsequent assembly) and postprocess element-level fields. For now we keep postprocessing for next week and focus only on getting the correct stiffness matrix.

Here the class describes an element combining extension and Euler-Bernoulli bending. A similar (or inherited) class could also be implemented for different element types (*e.g.* shear beam, Timoshenko beam, cable elements, etc). Here we also keep it simple by assuming elements are all arranged in a 2D plane.

However, the implementation is incomplete:
- The transformation matrix is missing in `__init__`, which is given in [](../lecture1/transformations.md). Make sure you take into account that a positive $\Delta z$ with a positive $\Delta x$ gives a negative angle $\alpha$. Make use of `numpy.arctan2` to return the angle between $-\pi$ and $\pi$, `numpy.arctan` returns an angle between $-\cfrac{\pi}{2}$ and $-\cfrac{\pi}{2}$, and therefore cannot distinguish between all four quadrants.
- The correct stiffness matrix for this extension-bending element coordinate system is missing in `stiffness`. You can derive the stiffness matrix yourself using pen and paper, SymPy or Maple, or copy the given stiffness matrix from [](../lecture1/other_elements.ipynb).
- We keep the functions which add a distributed load and compute the moments / displacements untouched for this week. Next week we'll implement those as well.


```{exercise} Workshop 1 - 2.1
:label: exercise2.1
:nonumber: true

Add the missing pieces to the code in `./matrixmethod/elements.py`, before you perform the checks below. Do you specify your stiffness matrix in the global or local coordinate system?
```

(exercise2_1_text_name)=
````{solution} exercise2.1
:class: dropdown

The stiffness matrix is specified in the global coordinate system.

For the code implementations see `./matrixmethod/elements.py`:
- [`__init__`](exercise2_1_py)
- [`stiffness`](exercise2_1_2_py)
````

+++

Whenever you make changes to your code in the `./matrixmethod/` folder, you need to reimport those. Instead of restarting the kernel, we use some magic ipython commands. Run the cell below once. Consequently, whenever you save your changes in one of the `.py`-files, it's automatically reloaded.

```{code-cell} ipython3
%load_ext autoreload
%autoreload 2
```

::::::{error}
Note that in this online book, you cannot make changes to multiple files simultaneously. These instructions are only applicable when you're working on this assignment locally; please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments).
::::::

```{exercise-start} Workshop 1 - 2.2
:label: exercise2.2
:nonumber: true

First, let's check the stiffness matrix for a beam which doesn't require rotation. Create a horizontal element with length $2$ and $EI=4$ and print both the transformation matrix and the stiffness matrix.

Do the matrices match with what you'd expect?
```

```{code-cell} ipython3
mm.Node.clear()
mm.Element.clear()
```

```{code-cell} ipython3
#YOUR CODE HERE

elem #= mm.Element(#YOUR CODE HERE

section = {}
section['EI'] #= YOUR CODE HERE
elem.set_section(section)

print(elem.T)
print(elem.stiffness())
```

```{exercise-end}
```

```{solution-start} exercise2.2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()
```

```{code-cell} ipython3
:tags: [thebe-init]

node1 = mm.Node(0,0)
node2 = mm.Node(2,0)

elem = mm.Element ( node1, node2 )

section = {}
section['EI'] = 4

elem.set_section (section)

print(elem)
print(elem.T)
print(elem.stiffness())
```

- The transformation matrix is identity which should be because the local coordinate system is aligned with the global coordinate system
- The values of the stiffness matrix are manually checked with the stiffness matrix from the slides and are correct.

```{solution-end}
```

```{exercise-start} Workshop 1 - 2.3
:label: exercise2.3
:nonumber: true

Now, create a vertical element with length $2$ and $EI=4$ and print the transformation and stiffness matrix.

Do the matrices match with what you'd expect?
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise2.3
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()

node1 = mm.Node(0,0)
node2 = mm.Node(0,1)
elem = mm.Element ( node1, node2 )

section = {}
section['EI'] = 1

elem.set_section (section)

print(elem)
print(elem.T)
print(elem.stiffness())
```

- The transformation matrix moves the extension terms to the column corresponding to the vertical displacements, which is correct.
- The transformation matrix moves the deflection terms to the column corresponding to the horizontal displacements, which is correct.
- The rotation terms are untouched, which is correct

```{solution-end}
```

```{exercise-start} Workshop 1 - 2.4
:label: exercise2.4
:nonumber: true

Now, create an element rotated in $120^{\circ}$ with length $2$ and print the transformation matrix.

Do the matrices match with what you'd expect?
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise2.4
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()

node1 = mm.Node(0,0)
node2 = mm.Node(-1,-np.sqrt(3))
elem = mm.Element ( node1, node2 )

print(elem)
print(elem.T)
```

- Again, the rotation term is not transformed, which is expected.
- A rotation of $120^{\circ}$ corresponds with a cosine term of $-0.5$ which is correctly calculated
- A rotation of $120^{\circ}$ corresponds with a sin term of $0.5\sqrt{3}$ which is correctly calculated

```{solution-end}
```

```{exercise-start} Workshop 1 - 2.5
:label: exercise2.5
:nonumber: true

Now, create an element rotated in $60^{\circ}$ with length $2$ and print the transformation matrix.

Do the matrices match with what you'd expect?
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise2.5
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()

node1 = mm.Node(0,0)
node2 = mm.Node(1,-np.sqrt(3))
elem = mm.Element ( node1, node2 )

print(elem)
print(elem.T)
```

- Again, the rotation term is not transformed, which is expected.
- A rotation of $60^{\circ}$ corresponds with a cosine term of $0.5$ which is correctly calculated
- A rotation of $60^{\circ}$ corresponds with a sin term of $0.5\sqrt{3}$ which is correctly calculated

```{solution-end}
```

```{exercise-start} Workshop 1 - 2.6
:label: exercise2.6
:nonumber: true

For the previous element, a global displacement vector $\mathbf{u}^{(e)} = \begin{bmatrix} 0 \\0 \\ 0 \\ \sqrt{3} \\ 1 \\ 0 \end{bmatrix}$ is given. What would be the local displacement vector $\bar{\mathbf{u}}^{(e)}$?

Check your answer using pen and paper. Tip: make a drawing instead of doing all the algebra.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise2.6
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

print(np.matmul(elem.T,np.array([0,0,0,np.sqrt(3),1,0])))
```

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/sketch.png
:align: center
```

```{solution-end}
```

## 3. The Constrainer class
This class is stored in [`./matrixmethod/constrainer.py`](../matrixmethod/constrainer.md)

This small class keeps track of which DOFs have prescribed displacements and takes care of applying these constraints to the global $\mathbf{K}$ and $\mathbf{f}$. For now we keep it simple and assume all constraints fix the DOF values to zero. Next week we will deal with non-zero prescribed values. 

However, the implementation is incomplete:
- The `constrain` function is incomplete, which should mimic the process of striking rows/columns of constrained DOFs and reduce the size of the system to be solved. Remember that `Constrainer` stores which DOFs are constrained in `self.dofs`, so **all the others** should be free. After gathering the free DOFs in an array, you will need to select the correct blocks of $\mathbf{K}$ and $\mathbf{f}$. For the stiffness matrix you will need the `np.ix_()` helper function (check its documentation [here](https://numpy.org/doc/stable/reference/generated/numpy.ix_.html))
- We keep the function which calculates supports reaction untouched for this week. Next week we'll implement that one as well.

+++

```{exercise} Workshop 1 - 3.1
:label: exercise3.1
:nonumber: true

Add the missing pieces to the code, before you perform the check below
```

````{solution} exercise3.1
:class: dropdown

For the code implementations see `./matrixmethod/constrainer.py`: [`constrain`](exercise3_1_py)
````

```{exercise-start} Workshop 1 - 3.2
:label: exercise3.2
:nonumber: true

Take the inclined element of exercise 2.5 and a bending stiffness of $1$. What happens if you invert $\mathbf{K}$? Now fix all degrees of freedom of the first node. What happens when you invert your 'constrained' $\mathbf{K}$? Are the dimensions of the 'constrained' $\mathbf{K}$ correct?
```

```{code-cell} ipython3
#YOUR CODE HERE

con = mm.Constrainer()

#YOUR CODE HERE

f = np.zeros (6) #empty load vector
Kff, Fff = con.constrain( K, F )
print(np.shape(np.linalg.inv(Kff)))
```

```{exercise-end}
```

```{solution-start} exercise3.2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

section = {}
section['EI'] = 1
elem.set_section (section)

k = elem.stiffness()
#print(np.linalg.inv(k))

con = mm.Constrainer()

con.fix_node (node1)

print(con)

f = np.zeros (6) #empty load vector
Kff, Fff = con.constrain( k, f )
print(np.shape(np.linalg.inv(Kff)))
```

- The inverted matrix of the original $\mathbf{K}$ gives a singular matrix, which is correct as the structure can perform any rigid body translation / rotation.
- The constrained matrix reduces $\mathbf{K}$ to only the free terms, which is not singular any more.
- The dimensions should match the number of free dofs, which is correct.

```{solution-end}
```

## 4. Full implementation extension bar

Having made our implementations, we now check them with two simple examples that serve as sanity checks. The first is a simple bar undergoing extension:

+++

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/extpointload.png
:align: center
:width: 200
```

With $EA = 1000$, $F = 100$ and $L = 1$.

+++

Use the code blocks below to set up and solve this problem using the classes above. The steps to follow are outlined below and short explanations/hints are given. Once you have a solution for the horizontal displacement of the node at the right end of the bar, compare it to the analytical solution you obtained in the first half of the course.

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()
```

```{exercise-start} Workshop 1 - 4.1
:nonumber: true
:label: exercise4.1

Create two nodes here. You can store them on a `list` or simply create them as two separate objects (*e.g.* `node1` and `node2`).
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.1
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

L  = 1

node1 = mm.Node (0,0)
node2 = mm.Node (L,0)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 4.2
:nonumber: true
:label: exercise4.2

Here we only have a single element, so there is no need to store it in a `list` yet. You are also going to need a `dict` defining the cross-section of the element.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

elem = mm.Element ( node1, node2 )

EA = 1000
section = {}
section['EA'] = EA

elem.set_section (section)
print(elem)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 4.3
:nonumber: true
:label: exercise4.3

Let's define the boundary conditions. We create an instance of the `Constrainer` class to deal with prescribed displacements. Take a look at its functions and inform if Node 1 is fully fixed.

You also need to pass the load $F$ on to Node 2. Check the member functions of `Node` to infer how that should be done.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.3
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

con = mm.Constrainer()

con.fix_node (node1)

F  = 100
node2.add_load ([F,0,0])

print(node2)
print(con)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 4.4
:nonumber: true
:label: exercise4.4

Now assemble the global stiffness matrix and force vector. Since we only have one element, there is no real assembly to be performed other than getting the stiffness matrix of the single element and storing the load at Node 2 in the correct positions of $\mathbf{f}$.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.4
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

global_k = elem.stiffness()
global_f = np.zeros (6)

global_f[3:6] = node2.p
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 4.5
:nonumber: true
:label: exercise4.5

Constrain the problem and solve for nodal displacements.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.5
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

Kff, Ff = con.constrain ( global_k, global_f )
u = np.matmul ( np.linalg.inv(Kff), Ff )
print(u)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 4.6
:nonumber: true
:label: exercise4.6

Finally, compare the displacement at the end of the bar with the one coming from the ODE solution. Note that since our element is already suitable for frames combining extension and bending, $\mathbf{u}$ has three entries. Which one is the entry that matters to us here? Did your solutions match? If so, that is a sign your implementation is correct. Can you use the function `full_disp` to obtain a vector of all displacements?
```

```{code-cell} ipython3
#EVENTUALLY YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise4.6
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

con.full_disp(u)
```

- The ODE solution is $ \cfrac{F L}{EA} = \cfrac{100 \cdot 1}{1000} = 0.1$ which equals the solution from the matrix method.
- Only the first term (displacement in horizontal direction) is relevant

```{solution-end}
```

## 5. Full implementation bending beam

In the first example above we tested our model under extension. But that does not really guarantee it will behave correctly in bending. That is the goal of this second sanity check. Let's solve the following problem:

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/cantilever.png
:align: center
:width: 200
```

Choose appropriate values yourself

When setting up and solving your model, note that we are now interested in $w$ displacements, our load is now vertical and the cross-section property driving our deformation is now $EI$. Good luck!

```{code-cell} ipython3
:tags: [thebe-init]

mm.Node.clear()
mm.Element.clear()
```

```{exercise-start} Workshop 1 - 5.1
:label: exercise5.1
:nonumber: true

Create nodes
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.1
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

L  = 1
node1 = mm.Node (0,0)
node2 = mm.Node (L,0)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 5.2
:label: exercise5.2
:nonumber: true

Create element
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

elem = mm.Element ( node1, node2 )

EI = 1000
section = {}
section['EI'] = EI

elem.set_section (section)
print(elem)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 5.3
:label: exercise5.3
:nonumber: true

Set boundary conditions
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.3
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

con = mm.Constrainer()

con.fix_node (node1)
F  = 100
node2.add_load ([0,F,0])
print(con)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 5.4
:label: exercise5.4
:nonumber: true

Assemble the system of equations.
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.4
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

global_k = elem.stiffness()
global_f = np.zeros (6)

global_f[0:3] = node1.p
global_f[3:6] = node2.p
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 5.5
:label: exercise5.5
:nonumber: true

Constrain the problem and solve for nodal displacements
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.5
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

Kff, Ff = con.constrain ( global_k, global_f )
u = np.matmul ( np.linalg.inv(Kff), Ff )
print(u)
```

```{solution-end}
```

```{exercise-start} Workshop 1 - 5.6
:label: exercise5.6
:nonumber: true

Check with the analytical solution

Did your solutions match? If so, your implementation is correct!
```

```{code-cell} ipython3
#EVENTUALLY YOUR CODE HERE
```

```{exercise-end}
```

```{solution-start} exercise5.5
:class: dropdown
```

The solution from a forget-me-not is $\cfrac{FL^3}{3EI} = \cfrac{100 \cdot 1^3}{3\cdot 1000} \approx 0.0333$ for the deflection and $\left|\cfrac{FL^2}{2EI}\right|= \left|\cfrac{100 \cdot 1^2}{2\cdot 1000}\right| = 0.05$ for the rotation which equals the solution from the matrix method.

```{solution-end}
```
