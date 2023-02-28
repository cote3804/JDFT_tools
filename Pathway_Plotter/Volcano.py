# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 11:23:45 2022

@author: coopy
"""

import sys
sys.path.append('C:\\Users\\coopy\\OneDrive - UCB-O365\\github\\JDFT_tools\\Pathway_Plotter')

import os
import matplotlib.pyplot as plt
import numpy as np
from Pathway import Pathway
from scipy.stats import linregress

class Volcano:
    def __init__(self ,materials: list, bias='0.00', reaction='N2R', path=None, filename='all_data.json'):
        self.materials = materials
        self.bias = bias
        self.reaction = reaction
        self.filename = filename
        self.path = path
        if path is None:
            self.path = os.getcwd()
            
        
        
    def call_pathway(self):
        #instantiate pathway object
        pathway = Pathway(path=self.path, filename=self.filename, reaction='N2R')
        
        #loop through materials and add them to the pathway 
        for material in self.materials:
            pathway.add_pathway(material, 'all', self.bias, 'min', free_energy=True)
        energies = pathway.energies()
        adatoms = pathway.get_adatom_energies()
        limiting_energies = pathway.get_limiting_energies()
        self.energies = energies
        self.adatoms = adatoms
        self.limiting_energies = limiting_energies
    
    def calculate_steps(self):
        energies = self.energies
        
        #create an empty array that is size (materials, reaction steps)
        #x direction specifies the material index in self.materials
        #y direction specifies the step in the reaction pathway in energies
        step_energy_data = np.zeros((len(self.materials),len(energies[self.materials[0]])-1))
        for imat, material in enumerate(self.materials):
            for istep, step in enumerate(energies[material].keys()):
                if istep == 0:
                    last_step = step
                    continue
                else:
                    step_energy = energies[material][step] - energies[material][last_step]
                    last_step = step
                    step_energy_data[imat,istep-1] = step_energy
            
        self.step_energy_data = step_energy_data
        return step_energy_data
    
    def get_regressions(self):
        step_energy_data = self.step_energy_data
        
        adatoms = []
        for material in self.materials:
            adatoms.append(self.adatoms[material])
        #regression data is an array that has the rows denominated by the reaction step
        #and the columns are slope and intercept. The regression is fit across the materials
        #so each regression is for all materials on one particular step.
        regression_data = np.zeros((len(step_energy_data[1,:]),2))
        for istep in range(len(step_energy_data[1,:])):
            regression = linregress(adatoms, step_energy_data[:,istep])
            regression_data[istep,:] = regression.slope, regression.intercept
        self.regression_data = regression_data
    
    def plot_volcano(self):
        # colors = ["#448aff","#1565c0","#009688","#8bc34a","#ffc107","#ff9800","#f44336","#ad1457","#448aff","#1565c0","#009688","#8bc34a","#ffc107","#ff9800","#f44336","#ad1457"]
        # colors = ["#152528","#08485e","#3cb0cd","#7cd5f3","#b4e7f8","#adffbc","#1bee9a","#21c063","#229631","#0c5a23"]
        colors = ["#d00000","#ffba08","#cbff8c","#8fe388","#1b998b","#3185fc","#5d2e8c","#46237a","#ff7b9c","#ff9b85"]

        self.call_pathway()
        self.calculate_steps()
        self.get_regressions()
        regression_data = self.regression_data
        
        states = self.energies[self.materials[0]].keys()
        steps = []
        for istate,state in enumerate(states):
            if istate == 0:
                last_state = state
                continue
            else:
                steps.append(f'{last_state}' + r'$ \rightarrow $' + f'{state}')
                last_state = state
        x = np.arange(-4,4)
        plt.figure(dpi=300)
        for istep in range(regression_data.shape[0]):
            plt.plot(x, 
                     -(x*regression_data[istep,0] + regression_data[istep,1]), 
                     label=steps[istep], color=colors[istep])
        location = [0.5, 0.5, 1, 0.5, 0.5] #TODO place labels so that they don't overlap
        for imat,material in enumerate(self.materials):
            # scatter_data.append(np.max(self.step_energy_data[imat]))
            plt.scatter(self.adatoms[material], -np.max(self.step_energy_data[imat]), color='k',)
            plt.annotate(f'{material.split("_")[0]} {material.split("_")[1]}', 
                         (self.adatoms[material] + 0.1*location[imat],-np.max(self.step_energy_data[imat]) - 0.5*location[imat]), 
                         color='k')
        plt.ylabel('$U_L (V)$')
        plt.xlabel('$ \Delta \Phi _{ads} (eV)$')
        plt.legend(bbox_to_anchor=(1.1, 1.05))
        
        
        
if __name__ == '__main__':
    materials = ['Au_211', 'Ru_211', 'Co_211_mag', 'Re_211', 'Ir_211']
    file = 'combined_N2R_all_data.json'
    path = os.path.normpath('C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\N2R_Scaling')
    
    volcano = Volcano(materials, path=path, filename=file)
    # volcano.call_pathway()
    # test = volcano.calculate_steps()
    # volcano.get_regressions()
    volcano.plot_volcano()