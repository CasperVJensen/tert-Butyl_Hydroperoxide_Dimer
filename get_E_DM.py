#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:16:14 2022

@author: Casper Vindahl Jensen
Department of Chemistry
University of Copenhagen

Note about script:
    
This script will take out the HF energy and dipole moment of all the single point energy calculation
log files. The log files are from a gaussian 16 output and the word searches in this script might
need to be updated if output files change format.


All electronic structure calculations are done using Gaussian 16. To use this script you need
log files of the single point electronic energy calculations
"""
#%%
import pandas as pd
import os
import glob
import math
import matplotlib.pyplot as plt
#%%
os.chdir('directory') #change to the directory of the folder with your log files
mol = str(a)     #replace "a" with name of your molecule. will be used in written textfile
#%%

for file in glob.glob("*.log"):
    print(file) #prints all files found with molecule name
    
files=glob.glob("*.log") #make list of all files with extension .log. can be replaced to match pattern.

q = [item.split("_",-1)[-1] for item in files]
q = [float(item.split(".log",-1)[0]) for item in q] #takes out the displacement from the file names and save it in a list
#%%
E = []      #Electronic energy
mux = []    #dipole moment X
muy = []    #dipole moment Y
muz = []    #dipole moment Z

def get(data, start, stop):                     #function to get number between two patterns
    splitstring = data.split(start, -1)[1]
    splitstring = splitstring.split(stop, -1)[0]
    return splitstring

d_const = 2.541747761 #conversion of dipolemoment from atomic units to Debye

for f in files:
    with open(f, 'r') as fff:
        data = fff.read().replace('\n', '')
        data = data.replace(' ', '')        #delete space and linebreak to make the string easier to recognize
       
    counter = data.count("Normaltermination")
    
    if counter != 0:
        
        Ei = float(get(data,"HF=", "\RMSD"))*219474.6 #find HF energy and convert to cm-1
        
        mu = get(data, "\Dipole=", "\Quadrupole=") #find dipole moment. patterns might need changing depending on outputfile
        muxi = float(mu.split(",",-1)[0])*d_const
        muyi = float(mu.split(",",-1)[1])*d_const
        muzi = float(mu.split(",",-1)[2])*d_const
    
    else:               #If gaussian calc did not terminate correctly return nan values
        Ei = math.nan
        muxi = math.nan
        muyi = math.nan
        muzi = math.nan
    
    fff.close()
    
    E.append(Ei)
    mux.append(muxi)
    muy.append(muyi)
    muz.append(muzi)
#%%
#Collect data in dataframe and write it to csv file. We collect displacement q, E and the dipole moment in each orientation
mydata = pd.DataFrame({"q" : q,
                       "E" : E,
                       "mux" : mux,
                       "muy" : muy,
                       "muz" : muz})
mydata = mydata[mydata['E'].notna()] #remove values that are nan
#%%
mydata = mydata[mydata['q']>=-0.45] #these two lines are hardcoded to chose maximum displacements manually.
mydata = mydata[mydata['q']<=1.475] #can be outcommented

mydata = mydata.sort_values(by=["q"])
mydata = mydata.reset_index(drop=True)

minE = mydata["E"].min()        #change energy into DELTA E
mydata["E"] = mydata["E"]-minE
#%%
mydata.to_csv(mol+".txt", index=None, sep='\t', mode='w')
#%%
plt.plot(mydata["q"],mydata["E"])

