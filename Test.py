# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:37:10 2023

@author: ala111
"""

# Packages
from classes import aromatic_forest, aromatic_form
import functions as fun


import sys
original_stdout = sys.stdout # Save a reference to the original standard output




#  Make sure that tau has its root in position 1!!!






def list_Lie_derivative(N,n,p,tau):
    print("List of the Lie derivatives")
    L_forests=fun.List_forests(N,n,p)
    for x in L_forests.alist:
        print(x[0].dic)
        g = aromatic_forest(x[0].dic).wedge()
        res1=g.Lie(tau)
        res1.aro_print()
        print()


def total_list_Lie_derivative(N,n,p,M):
    print("List of the Lie derivatives - (N,M)=",N,M)
    L_forests=fun.List_forests(N,n,p)
    L_tau=fun.List_forests(M,1,0)
    for y in L_tau.alist:
        print("tau=",end='')
        y[0].draw()
        for x in L_forests.alist:
            print("gamma=",end='')
            x[0].draw()
            print()
            g = x[0].wedge()
            res=g.Lie(y[0].wedge()).delta_V()
            # res.draw()
            res.aro_print()
            print()


def search_conservation_law(N,n,p,M):
    print("List of the Lie derivatives - (N,M)=",N,M)
    L_forests=fun.List_forests(N,n,p)
    L_tau=fun.List_forests(M,1,0)
    for x in L_forests.alist:
        print("gamma=",end='')
        x[0].draw()
        for y in L_tau.alist:
            print("tau=",end='')
            y[0].draw()
            print()
            g = x[0].wedge()
            res=g.Lie(y[0].wedge())
            # res.simplify_1_loop()
            res.draw()
            print()


with open('results.txt', 'w') as f:
    sys.stdout = f
# L_tau=fun.List_forests_div_free(4,1,0)
# for y in L_tau.alist:
#     print("tau=",end='')
#     y[0].draw()
#     print()
#     g = aromatic_forest({1: ['r', 1], 2: ['r',2], 3: ['v',2]}).wedge().d_H()
#     g.ext_mult(2)
#     res=g.Lie(y[0].wedge())
#     res.ext_mult(-1)
#     res.add(y[0].Lie(g))
#     # res.d_H().draw()
#     # print()
#     res.d_H_div_free().draw()
#     print()
    for M in [1,2,3]:
        total_list_Lie_derivative(2,0,0,M)
        print()
        total_list_Lie_derivative(3,0,0,M)
        print()
    sys.stdout = original_stdout











# total_list_Lie_derivative(2,0,0,3)







def list_H_Lie(g,M):
    print("List of the Lie homotopy")
    L_tau=fun.List_forests(M,1,0)
    for y in L_tau.alist:
        print("tau=",end='')
        y[0].draw()
        print()
        res=g.h_H_Lie(y[0].wedge())
        res.draw()
        print()


# g = aromatic_forest({1: ['r', 1], 2: ['r',2], 3: ['v',2]}).wedge().d_H()
# g.ext_mult(2)
# list_H_Lie(g,4)



def list_H_Lie_test(N,M):
    L_forests=fun.List_forests(N,2,0)
    L_tau=fun.List_forests(M,1,0)
    for x in L_forests.alist:
        print("gamma=",end='')
        x[0].draw()
        g = x[0].wedge()
        dg=g.d_H()
        for y in L_tau.alist:
            print("tau=",end='')
            y[0].draw()
            print()
            res=dg.h_H_Lie(y[0].wedge())
            res_aux=g.Lie(y[0].wedge())
            res_aux.ext_mult(-1)
            res.add(res_aux)
            res.draw()
            print()


# list_H_Lie_test(5,3)




# L_forests=fun.List_forests(6,2,0)
# for x in L_forests.alist:
#     print("gamma=",end='')
#     x[0].draw()
#     print()
#     g = x[0].wedge()
#     g.h_H().d_H().draw()
#     print()

































