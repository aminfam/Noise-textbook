---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# `constrainer.py`

::::::{versionadded} v2025.2.0 After workshop 2
Solutions workshop 2 and additional assignments in text and downloads 
::::::

::::::{versionadded} v2025.1.0 After workshop 1
Solutions workshop 1 in text and downloads 
::::::

::::::{attention}
This page shows a preview of the `matrixmethod` package. Please fork and clone the practice assignments to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)

After each workshop, the solution will be added to this preview and to the [GitHub-repository](https://github.com/CIEM5000-2025/practice-assignments)
::::::

```{custom_download_link} constrainer.py
:text: ".py"
:replace_default: "False"
```

```{custom_download_link} ./matrixmethod_solution/constrainer.py
:text: ".py solution workshop 1, 2 and additonal assignments"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments/tree/solution_workshop_1
:text: "All files practice assignments with solutions workshop 1"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments/tree/solution_workshop_2
:text: "All files practice assignments with solutions workshop 1"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments/tree/solution_additional_exercises
:text: "All files practice assignments with solutions additional exercises"
```

```{code-cell} ipython3
import numpy as np
```

```{code-cell} ipython3
class Constrainer:
    """
    A class that represents a constrainer for fixing degrees of freedom in a structural analysis.

    Attributes:
        cons_dofs (list): A list of constrained degrees of freedom.
        cons_vals (list): A list of corresponding constraint values.

    Methods:
        fix_dof: Fixes a degree of freedom at a specific value.
        fix_node: Fixes all degrees of freedom of a node.
        full_disp: Combines the displacements of free and constrained degrees of freedom.
        constrain: Applies the constraints to the stiffness matrix and load vector.
        support_reactions: Calculates the support reactions based on the constrained displacements.
    """

    def __init__(self):
        """
        Initializes a new instance of the Constrainer class.

        Attributes:
            cons_dofs (list): A list of constrained degrees of freedom.
            cons_vals (list): A list of corresponding constraint values.  
        """
        self.cons_dofs = []
        self.cons_vals = []

    def fix_dof (self, node, dof, value = 0):
        """
        Fixes a degree of freedom  at a specific value.

        Args:
            node (Node): The node object.
            dof (int): The index of the degree of freedom to fix.
            value (float, optional): The value to fix the degree of freedom at. Defaults to 0.
        """
        self.cons_dofs.append(node.dofs[dof])
        assert value == 0, "Only zero values are supported for now."
        self.cons_vals.append(value)
 
    def fix_node (self, node):
        """
        Fixes all degrees of freedom of a node.

        Args:
            node (Node): The node object.
        """
        for dof in [0,1,2]:
            self.fix_dof (node, dof)    

    def full_disp (self,u_free):
        """
        Combines the displacements of free and constrained degrees of freedom.

        Args:
            u_free (numpy.ndarray): The displacements of the free degrees of freedom.

        Returns:
            numpy.ndarray: The combined displacements of all degrees of freedom.
        """
        u_full = np.zeros(len(self.free_dofs) + len(self.cons_dofs))
        
        u_full[self.free_dofs] = u_free
        u_full[self.cons_dofs] = self.cons_vals
        
        return u_full
    
    def constrain (self, k, f):
        """
        Applies the constraints to the stiffness matrix and load vector.

        Args:
            k (numpy.ndarray): The stiffness matrix.
            f (numpy.ndarray): The load vector.

        Returns:
            tuple: A tuple containing the stiffness matrix corresponding to free dofs and the corresponding load vector.
        """
        self.free_dofs = [i for i in range(len(f)) if i not in self.cons_dofs]
        
        Kff #= k[np.ix_(YOUR CODE HERE)]
        Ff # YOUR CODE HERE

        return Kff, Ff
```

+++

(exercise3_1_py)=
```{solution-start} exercise3.1
:class: dropdown
```

```{code-cell} ipython3

        self.free_dofs = [i for i in range(len(f)) if i not in self.cons_dofs]
        
        Kff = k[np.ix_(self.free_dofs,self.free_dofs)]
        Ff = f[self.free_dofs]

        return Kff, Ff
```

```{solution-end}
```

+++

(2_exercise3_1_py_1)=
```{solution-start} 2_exercise3.1
:class: dropdown
```

```{code-cell} ipython3
        self.free_dofs = [i for i in range(len(f)) if i not in self.cons_dofs]
        
        Kff = k[np.ix_(self.free_dofs,self.free_dofs)]
        Kfc = k[np.ix_(self.free_dofs,self.cons_dofs)]
        Ff = f[self.free_dofs]

        return Kff, Ff - np.matmul(Kfc,self.cons_vals)
```

```{solution-end}
```

+++
        
```{code-cell} ipython3
    def support_reactions (self,k,u_free,f):       
        """
        Calculates the support reactions based on the constrained displacements.

        Args:
            k (numpy.ndarray): The stiffness matrix.
            u_free (numpy.ndarray): The displacements of the free degrees of freedom.
            f (numpy.ndarray): The load vector.

        Returns:
            numpy.ndarray: The support reactions.
        """
        #YOUR CODE HERE
        
        return #YOUR CODE HERE
```

+++

(2_exercise3_1_py_2)=
```{solution-start} 2_exercise3.1
:class: dropdown
```

```{code-cell} ipython3
        Kcf = k[np.ix_(self.cons_dofs,self.free_dofs)]
        Kcc = k[np.ix_(self.cons_dofs,self.cons_dofs)]

        return np.matmul(Kcf,u_free) + np.matmul(Kcc,self.cons_vals) - f[self.cons_dofs]
```

```{solution-end}
```

+++
        
```{code-cell} ipython3
    def __str__(self):
        """
        Returns a string representation of the Constrainer object.

        Returns:
            str: A string representation of the Constrainer object.
        """
        return f"This constrainer has constrained the degrees of freedom: {self.cons_dofs} with corresponding constrained values: {self.cons_vals})"
```
