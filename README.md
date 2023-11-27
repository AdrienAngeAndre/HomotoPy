# HomotoPy
Python package for manipulating the aromatic forms on the aromatic bicomplex

This package contains some code in Julia and Matlab for obtaining the results and figures presented in the article:

[ ] A. Laurent, Robert I. McLachlan, Hans Z. Munthe-Kaas, and Olivier Verdier,
      The aromatic bicomplex for the description of divergence-free aromatic forms and volume-preserving integrators.
      Forum of Mathematics Sigma.
      
The parts of the code related to the Lie derivative are discussed in:

[ ] A. Laurent,
      The Lie derivative and Noether's theorem on the aromatic bicomplex.
      Journal of Computational Dynamics.

Please cite the above papers when using this package for research ! :-)


The Python package HomotoPy allows to manipulate aromatic forests and forms with the operators defined in the paper.

The file classes.py defines the classes corresponding to the aromatic forests and forms, as well as the operators of the paper.
The file functions.py gives a handful of practical functions to manipulate the sets Omega_{n,p}^N.
The file homotopy.py gives some examples to manipulate the forms, compute their derivatives, their homotopy, in both the standard context and the divergence-free context.
The file aromatic_bicomplex.py prints the aromatic bicomplex in both contexts, the augmented column, as well as the dimensions of the spaces involved.

The forests are implemented in the following way. We number the vertices from 1 to N-p, the covertices from -1 to -p.
Then, an aromatic forest is associated to a dictionary of the form gdic = {1: ['r', 1],2: ['r', 2],3: ['v', 1],-1: ['v', 3],-2: ['r', 3]}.
In this example, the node 1 is the first root (signaled by the character 'r'), -2 is the third root, and 3 is a vertex linked to 1.
The associated aromatic forest is g = aromatic_forest(gdic).
The class aromatic form is a list of aromatic forest with coefficients (a vector space). All the operations allowed on the forests are extended to the forms.
We refer to homotopy.py for examples and to classes.py for the list of all available functions.

Note: due to the complexity in N! of the implementation of forms in the software, there is a limit around N=8 for manipulating the aromatic forms at the present time.
