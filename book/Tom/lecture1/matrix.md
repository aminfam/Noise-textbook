# Combine elements using matrix formulation

In the previous [chapter](./combine.md) you've seen how to combine multiple elements to arrive to a linear system of equation to solve the nodal displacements. However, the nodal displacements were part of $\mathbf{f}^e$, which doesn't allow us to solve the vector formulation directly.

::::::{topic} Learning objective
You'll look how to define local and global matrices and vectors of separately: unknown displacements vector, external force vector and elements stiffness matrices.
::::::

## Local stiffness matrix, displacement vector and force vector
In [](./single_element.md) you've derived the following two equations relating the nodal forces and displacements:

- $ F_1 = \cfrac{EA}{\ell}u_1-\cfrac{EA}{\ell}u_2$
- $ F_2 = -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell}u_2$

These two equations can be rewritten in a matrix formulation:

$$\cfrac{EA}{\ell}\begin{bmatrix} 1&-1\\-1&1 \end{bmatrix}\begin{bmatrix}u_1\\u_2\end{bmatrix} = \begin{bmatrix}F_1\\F_2\end{bmatrix}$$

This can represented as follows too:

$$
\mathbf{K}^{(e)}\mathbf{u}^{(e)}=\mathbf{f}^{(e)}
$$

$\mathbf{K}^{(e)}$ is known as the local stiffness matrix, $\mathbf{u}^{(e)}$ the local displacement vector and $\mathbf{f}^{(e)}$ the local force vector.

## Global stiffness matrix, displacement vector and force vector
Now let's see how to combine multiple elements using our new matrix formulation. In [](./combine.md) you've derived the following three equation relating the nodal displacements with external forces:

- $\sum F_1 = 0 \Rightarrow \cA{-\cfrac{EA_1}{\ell_1}}u_1 \cA{ + \cfrac{EA_1}{\ell_1}}u_2 + H=0$
- $\sum F_2 = 0 \Rightarrow \cA{\cfrac{EA_1}{\ell_1}}u_1 \cA{ - \cfrac{EA_1}{\ell_1}}u_2 \cB{ -\cfrac{EA_2}{\ell_2}}u_2 \cB{ + \cfrac{EA_2}{\ell_2}}u_3 =0$
- $\sum F_3 = 0 \Rightarrow \cB{\cfrac{EA_2}{\ell_2}}u_2 \cB{ - \cfrac{EA_2}{\ell_2}}u_3 + F=0$ 

These equation can be represented in a matrix formulation:

$$\begin{bmatrix}\cA{\cfrac{EA_1}{\ell_1}}&\cA{-\cfrac{EA_1}{\ell_1}}&0\\\cA{-\cfrac{EA_1}{\ell_1}}&\cA{\cfrac{EA_1}{\ell_1}}+\cB{\cfrac{EA_2}{\ell_2}}&\cB{-\cfrac{EA_2}{\ell_2}}\\0&\cB{-\cfrac{EA_2}{\ell_2}}&\cB{\cfrac{EA_2}{\ell_2}}\end{bmatrix}\begin{bmatrix}u_1\\[12pt]u_2\\[12pt]u_3\end{bmatrix}=\begin{bmatrix}H\\[12pt]0\\[12pt]F\end{bmatrix} $$

This can represented as follows too:

$$\mathbf{K}\mathbf{u}=\mathbf{f}$$


$\mathbf{K}$ is known as the global stiffness matrix, $\mathbf{u}$ the global displacement vector and $\mathbf{f}$ the glocal force vector.

As you might see, the local matrices are visible as 'blocks' in the global matrix, while the displacement and force vector are very clean. It seems to be possible to define the global stiffness matrix directly without manually evaluating the force equilibrium at every node! You'll look at that in the next chapter!
