# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:59:06 2022

@author: Adrien Ange Andre Laurent
"""

# Packages
import itertools as it
import numpy as np
import math
from classes import aromatic_forest, aromatic_form




# Auxiliary functions
def sub_lists(l): 
    base = []   
    lists = [base] 
    for i in range(len(l)): 
        orig = lists[:] 
        new = l[i] 
        for j in range(len(lists)): 
            lists[j] = lists[j] + [new] 
        lists = orig + lists
    lists.remove([])
    return lists 




# Main functions
def List_forests(N,n,p):
    '''
    List of all forests in F_{n,p}^N (up to permutation of covertices or roots).
    Output as an aromatic form containing all such forests.
    '''
    res=aromatic_form()
    resf=aromatic_form()
    for k in range (max(0,n+p-N),min(n,p)+1):
        list_roots=[j for j in range(-k,0)]+[j for j in range(1,n-k+1)]
        list_vertices=[j for j in range(n-k+1,N-p+1)]
        list_covertices=[j for j in range(-p,0)]
        list_nodes=list_covertices+[j for j in range(1,N-p+1)]
        list_non_roots=[j for j in range(-p,-k)]+list_vertices
        graftings=[item for item in it.product(list_nodes,repeat=N-n)]
        for f in graftings:
            gdic={}
            for j in range(1,n+1):
                gdic[list_roots[j-1]]=['r',j]
            for j in range(0,N-n):
                gdic[list_non_roots[j]]=['v',f[j]]
            if aromatic_forest(gdic).wedge().alist!=[]:
                res.add_one_elmt(aromatic_forest(gdic),1)
        res.simplify()
    leng=len(res.alist)
    for i in range(0,leng-1):
        boo=True
        for j in range(i+1,leng):
            aux=aromatic_form()
            aux.add_one_elmt(res.alist[i][0], 1)
            aux.add_one_elmt(res.alist[j][0], 1)
            if aux.wedge().alist==[]:
                boo=False
                break
            aux.add_one_elmt(res.alist[j][0], -2)
            if aux.wedge().alist==[]:
                boo=False
                break
        if boo:
            resf.add_one_elmt(res.alist[i][0], 1)
    if leng>0:
        resf.add_one_elmt(res.alist[leng-1][0], 1)
    return resf


def List_forests_div_free(N,n,p):
    '''
    List of all forests in F_{n,p}^N in the divergence-free case (up to permutation of covertices or roots).
    Output as an aromatic form containing all such forests.
    '''
    res=List_forests(N,n,p)
    res.simplify_1_loop()
    return res


def aromatic_bicomplex(N):
    '''
    Prints the generators of the aromatic bicomplex.
    '''
    print("Aromatic bicomplex")
    for p in range(0,N+1):
        for n in range(0,N+1):
            L_forests=List_forests(N,n,p)
            print("\n","F(",n,",",p,")")
            for x in L_forests.alist:
                print(x[0].dic)
                
                
def aromatic_bicomplex_div_free(N):
    '''
    Prints the generators of the aromatic bicomplex in the divergence-free context.
    '''
    print("Divergence-free aromatic bicomplex")
    for p in range(0,N+1):
        for n in range(0,N+1):
            L_forests=List_forests_div_free(N,n,p)
            print("\n","F(",n,",",p,")")
            for x in L_forests.alist:
                print(x[0].dic)
    
    
def dimension_aromatic_bicomplex(N):
    '''
    Prints the dimensions of the space involved in the aromatic bicomplex in a matrix.
    '''
    A=np.zeros((N+1,N+1))
    print("Dimensions of the aromatic bicomplex")
    for p in range(0,N+1):
        for n in range(0,N+1):
            L_forests=List_forests(N,n,p)
            A[N-p][N-n]=len(L_forests.alist)
    print(A)


def dimension_aromatic_bicomplex_div_free(N):
    '''
    Prints the dimensions of the space involved in the divergence-free aromatic bicomplex in a matrix.
    '''
    A=np.zeros((N+1,N+1))
    print("Dimensions of the divergence-free aromatic bicomplex")
    for p in range(0,N+1):
        for n in range(0,N+1):
            L_forests=List_forests_div_free(N,n,p)
            A[N-p][N-n]=len(L_forests.alist)
    print(A)
    
    
def augmented_column(N):
    '''
    Prints the generators of the extra column of the augmented bicomplex. They do not form a basis in general.
    '''
    print("Augmented column")
    for p in range(1,N+1):
        L_forests=List_forests(N,0,p)
        print("\n","I(",p,")")
        for x in L_forests.alist:
            print("Forest",x[0].dic)
            print("Image by I")
            x[0].wedge().I().aro_print()


def list_Lie_derivative(N,n,p,tau):
    print("List of the Lie derivatives")
    L_forests=List_forests(N,n,p)
    for x in L_forests.alist:
        print(x[0].dic)
        g = aromatic_forest(x[0].dic).wedge()
        res1=g.Lie(tau)
        res1.aro_print()
        print()


def total_list_Lie_derivative(N,n,p,M):
    print("List of the Lie derivatives - (N_gamma,N_tau)=",N,M)
    L_forests=List_forests(N,n,p)
    L_tau=List_forests(M,1,0)
    for y in L_tau.alist:
        print("tau=",end='')
        y[0].draw()
        for x in L_forests.alist:
            print("gamma=",end='')
            x[0].draw()
            print()
            g = x[0].wedge()
            res=g.Lie(y[0].wedge()).delta_V()
            res.draw()
            # res.aro_print()
            print()








