# Combine elements

In the previous [chapter](./single_element.md) you've seen how to set up the force-displacement relations for a single element. However, to solve for complete structures you'll need to combine multiple elements.

::::::{topic} Learning objective
You'll look into how to combine elements to represent a full structure.
::::::

Let's reconsider the extension problem with two fields:

```{figure} assembly5.svg
:name: assembly5
:align: center

Extension bar with nodal load
```

Both elements have the same solution as derived in [](./single_element.md). To keep track of the different forces, we'll use the subscript $F_1$ and $F_2$ for defining the left- and right-end force of an element, the superscript $F^{(1)}$ and $F^{(2)}$ to define the element itself. The nodal displacement are numbered with a subscript $u_1$, $u_2$ and $u_3$.

Now, let's draw free body diagrams of the nodes and the elements itself. The forces acting on the ends of the elements, act in the opposite direction on the nodes. On node $1$ and $3$ additional forces are present: $H$ for the reaction force in $1$ which is assumed in the positive direction and $F$ for the external force:

```{figure} nodaleq2.svg
:name: nodaleq2
:align: center

Free body diagrams of nodes and elements
```

Horizontal equilibrium of each of the nodes gives:

- $\sum F_1 = 0 \Rightarrow \cA{-\cfrac{EA_1}{\ell_1}}u_1 \cA{ + \cfrac{EA_1}{\ell_1}}u_2 + H=0$
- $\sum F_2 = 0 \Rightarrow \cA{\cfrac{EA_1}{\ell_1}}u_1 \cA{ - \cfrac{EA_1}{\ell_1}}u_2 \cB{ -\cfrac{EA_2}{\ell_2}}u_2 \cB{ + \cfrac{EA_2}{\ell_2}}u_3 =0$
- $\sum F_3 = 0 \Rightarrow \cB{\cfrac{EA_2}{\ell_2}}u_2 \cB{ - \cfrac{EA_2}{\ell_2}}u_3 + F=0$ 

The algoritmic approach should become visible now!

These equations could also be regarding in a vector formulation. If $\mathbf{f}^e$ represents a vector with all the forces acting on each of the nodes coming from a single element and $\mathbf{f}_\text{nodal}$ represents a vector of all the nodal forces not coming from the elements (the reaction force $H$ and external force $F$), these equations can be simplified as:

$$\begin{align} -\sum_e\mathbf{f}^e + \mathbf{f}_\text{nodal}& = \mathbf{0}\\
\sum_e\mathbf{f}^e& = \mathbf{f}_\text{nodal} \end{align}$$

All that's needed now is to solve our linear set of equations for our unknown nodal displacements. Luckily the amount of equations equals the amount of unknowns, so you should have no problem solving this! Take into account that one displacement is already known: $u_1 = 0$. Without this, the matrix is singular.

However, solving our vector formulation is not trivial. $\mathbf{f}^e$ still contains our unknown nodal displacements hidden inside the vector. Let's split the vector on the next page to reach our final form of the matrix method formulation.