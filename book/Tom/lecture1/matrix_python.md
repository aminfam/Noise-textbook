# Implementation in Python

The matrix method is implemented in an (incomplete) python package which you'll extend in the workshops and graded assignments. 

::::::{topic} Learning objective
You'll look into the structure of the matrix method python package.
::::::

The method is broken down as follows:
- A list of nodes floating in space with loads and DOFs associated to them

```{figure} nodes.svg
:name: nodes
:align: center
```

- A list of elements defined by linking two nodes together

```{figure} elements.svg
:name: elements
:align: center
```

- A constrainer to apply Dirichlet boundary conditions

```{figure} constrain.svg
:name: constrain
:align: center
```

With this in mind, the python package is set up as an object-oriented code which each of the items above as a class, containing multiple attributes and functions. The docstring provide an introduction to the classes, attributes and functions. These will be treated later in the workshops

The three classes work together: nodes defined with the `Node` class are an input for the `Element` class. The `Constrainer` takes the nodes defined with the `Node` class and constrains specific degree of freedom to comply with Dirichlet boundary conditions.

```{figure} code.svg
:name: code
:align: center

Structure of python package
```

The global stiffness matrix, nodal displacement vector and nodal force vector are not defined anywhere in each of the classes (only local stiffness matrices are defined in the `Element` class). So, these need to defined separately.

+++

## Setting up global equations with python package
The steps provide in the [previous page](./directly.md) can be implemented using the matrix method package as follows:

### 1. Identify degrees of freedom
The `Node` class automatically keeps track of the amount of degrees of freedom based on the amount of instances created with that class. For every node it will count 3 degrees of freedom. So far you've only seen one degree of freedom per node, but soon more will follow!

### 2. Initialize the system with zeros
As said before, the global stiffness matrix, nodal displacement vector and nodal force vector are not defined anywhere in each of the classes. Therefore, you'll to create empty matrices and vectors yourself in the form of numpy arrays.

### 3. Assemble stiffness, element by element
The `Element` class can create a local stiffness matrix based on the `Node`-instances which are provided as an input. You'll need some smart indexing to map the local stiffness matrix to the correct location of the global stiffness matrix. Again, this is not part of any of the provided classes.

### 4. Apply external loads
Loads can be added to the instances of the `Node` class. Together with the indexing of degrees of freedom inside these instances, the loads can be mapped to the correct location in the global nodal force vector.

### 5. Apply prescribed displacements
The prescribed displacements can be applied using the `Constrain` class. It allows you to fix the degrees of freedom, which you'll need in the next step. It doesn't take care of the reaction force belonging to the prescribed displacement.

### 6. Solve for the unknown nodal displacements
Again, the `Constrain` class can help you here. As it knows which degrees of freedom are fixed, it can give you a reduced global stiffness matrix and global nodal force vector which you can use yourself to solve for the free global nodal displacements
