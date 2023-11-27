# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 09:20:04 2021

@author: Adrien Ange Andre Laurent
"""

# Packages
from classes import aromatic_forest, aromatic_form
import functions as fun




print("1) Computing the derivatives")
# Definition of a form
gdic = {1: ['r', 1],2: ['r', 2],3: ['v', 1],-1: ['v', 3]}
g = aromatic_forest(gdic).wedge()
print('Form g')
g.aro_print()


# Derivatives of g
print('dH g')
g.d_H().aro_print()

print('dV g')
g.d_V().aro_print()

print('dH_div_free g')
g.d_H_div_free().aro_print()

print('dV_div_free g')
g.d_V_div_free().aro_print()




print("\n")
print("2) Testing the horizontal homotopy identity")
# Definition of a form
gdic = {1: ['r', 1],2: ['r', 2],3: ['v', 1],-1: ['v', 3]}
g = aromatic_forest(gdic).wedge()
print('Form g')
g.aro_print()


# Test horizontal homotopy -- divergence-free
g.simplify_1_loop()
gaux=aromatic_forest(gdic).wedge()
gaux.ext_mult(-1)
gaux.simplify_1_loop()

print('hH o dH + dH o hH - Id + correction')
# It returns the empty aromatic form as expected.
res1=g.d_H_div_free().h_H_div_free()
res2=g.h_H_div_free().d_H_div_free()
corr=g.correction_horizontal_homotopy_div_free()
gaux.add(res1)
gaux.add(res2)
gaux.add(corr)
gaux.simplify()
gaux.aro_print()




print("\n")
print("3) Comparison of the two horizontal homotopy operators on Omega_0")
# Alternative homotopy operator for a form in Omega_0
gdic = {1: ['v', 1],2: ['v', 1],3: ['v', 1],4: ['v', 1]}
g = aromatic_forest(gdic).wedge()
print('Form g')
g.aro_print()

# We recall that the two homotopy operators do not coincide in general.
#They coincide up to a divergence-free term.
print('Homotopy operator with the Euler operators')
g.h_H().aro_print()

print('Homotopy operator with the IBP process')
g.alternative_h_H().aro_print()




print("\n")
print("4) Generators of the solenoidal forms")
N=4
L_forests=fun.List_forests(N,2,0)
for x in L_forests.alist:
    print("\n",x[0].dic,"\n")
    x[0].d_H().aro_print()




print("\n")
print("5) Example of a non-trivial element of Ker(dH)")
alist=[
        [aromatic_forest({1: ['v', 1],2: ['v', 1],3: ['v', 2],4: ['r', 1],5: ['r', 2],6: ['v', 5]}),1],
        [aromatic_forest({1: ['v', 1],2: ['v', 1],3: ['r', 1],4: ['r', 2],5: ['v', 4],6: ['v', 5]}),-1],
        [aromatic_forest({1: ['v', 1],2: ['r', 1],3: ['v', 2],4: ['r', 2],5: ['v', 4],6: ['v', 5]}),1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 2],4: ['v', 2],5: ['v', 4],6: ['v', 5]}),1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 1],4: ['v', 3],5: ['v', 2],6: ['v', 2]}),-1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 2],4: ['v', 2],5: ['v', 3],6: ['v', 4]}),-1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 2],4: ['v', 3],5: ['v', 3],6: ['v', 5]}),-1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 1],4: ['v', 2],5: ['v', 2],6: ['v', 5]}),1],
        [aromatic_forest({1: ['r', 1],2: ['r', 2],3: ['v', 1],4: ['v', 2],5: ['v', 4],6: ['v', 4]}),1],
        [aromatic_forest({1: ['r', 1],2: ['v', 1],3: ['v', 2],4: ['r', 2],5: ['v', 6],6: ['v', 5]}),1],
        [aromatic_forest({1: ['r', 2],2: ['v', 1],3: ['r', 1],4: ['v', 5],5: ['v', 4],6: ['v', 5]}),1],
        [aromatic_forest({1: ['r', 2],2: ['v', 1],3: ['r', 1],4: ['v', 6],5: ['v', 4],6: ['v', 5]}),1]
        ]
g=aromatic_form(alist).wedge()
print("Aromatic form g")
g.aro_print()
print("\n","d_H(g)=0")
g.d_H().aro_print()




print("\n")
print("6) Lie derivatives of the first aromatic vector fields and Lagrangians")
fun.total_list_Lie_derivative(2,0,0,2)








