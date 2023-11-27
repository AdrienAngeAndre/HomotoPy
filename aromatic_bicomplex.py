# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:52:05 2022

@author: Adrien Ange Andre Laurent
"""

# Packages
from classes import aromatic_forest, aromatic_form
import functions as fun




# Parameters
N=3


print("1) Standard case")
fun.aromatic_bicomplex(N)

print("\n")
fun.augmented_column(N)

print("\n")
fun.dimension_aromatic_bicomplex(N)




print("\n")
print("2) Divergence-free case")
fun.aromatic_bicomplex_div_free(N)

print("\n")
fun.dimension_aromatic_bicomplex_div_free(N)



