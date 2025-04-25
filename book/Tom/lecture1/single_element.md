# Force-displacement relations single extension element

In the previous [chapter](./displacement.md) you've seen how you can solve structure using nodal displacements. However, the approach was still very problem-dependent. As proposed, the matrix method solves this by defining a default element which can be solved for a priori.

::::::{topic} Learning objective
Let's define the force-displacement relations for a single element loading in extension in terms of nodal displacements.
::::::

We'll do that using differential equations. Later you'll see how to get the same result using shape functions.

Let's consider the most simple extension element:

```{figure} extensionelement.svg
:name: extensionelement
:align: center

Single extension element
```

The same approach is used as in [](./recap.ipynb). However, the boundary conditions now are defined in terms of unknown nodal displacements:

- $u(0) = u_1$
- $u(\ell)=u_2$

This results in:

- $C_1 = \cfrac{u_2-u_1}{\ell}$
- $C_2 = u_1$

Effectively, we replaced the unknown integration constants by unknown nodal displacements. However, this will prove to be useful because these nodal displacements have a clear physical meaning and it will allow us to 'glue' elements together as other elements will be connected with the same nodal displacement.

Using our new formulation of unknown, the continuous distributions for the displacement and section force can be evaluated too:
- $u(x) = u_1\left(1-\cfrac{x}{\ell}\right) + u_2\cfrac{x}{\ell}$
- $N = -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell}u_2$

We'll combine elements using force equilibrium, therefore, the force at the ends of the elements are of main interst. The section force $N$ is derived above, which is not the same as the force $F_1$ and $F_2$ defined as positive to the right. For the node on the right, the direction of $N$ coincides with the positive direction of $F_1$. However, on the left-hand-side, the sign flips. Leading to our force-displacement relation:
- $ F_1 = - N = \cfrac{EA}{\ell}u_1-\cfrac{EA}{\ell}u_2$
- $ F_2 = N = -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell}u_2$