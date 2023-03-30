#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:55:04 2022

@author: Casper Vindahl Jensen
Department of Chemistry
University of Copenhagen

Note about script:
    
This script will displace the H in the chosen OH-bond and create files for each displacement.


All electronic structure calculations are done using Gaussian 16. To use this script you need
an optimized equilibrium structure of the molecule and a inputfile for a single point energy calculation
on that structure.

"""
#%%
import numpy as np
import pandas as pd
import os
#%%

filename = "filename.extension" #Give filename of input file with extension 
                                #NEEDS TO HAVE "eq" IN FILENAME WITHOUT QUOTATIONS

with open(filename, "r") as file:
    filedata = file.read()
#%%
steps = np.arange(-0.5, 1.525, 0.025) # give the range and increments in Ångstrøm to displace the H
                                        #from the Eq-structure
#%%
eq_distance = a #replace a with the eq-distance in OH bond

for step in steps:
    data_i = filedata.replace(str(eq_distance), str(round(eq_distance+step,8)))

    with open(filename.replace("eq", str(round(step,3))), 'w') as file:
        file.write(data_i)
        file.close()
        
        #Will write files with eq replaces with displacement distance.
