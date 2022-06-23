# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:12:04 2022

@author: coopy
"""
import os
import json
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.electronic_structure import dos
from pymatgen.electronic_structure.core import Spin
import numpy as np


from dos_helper import set_rc_params, get_plot, fill_zeros
set_rc_params()

class DOS_visualize:
    def __init__(self, files, Atoms=[], orbitals=['all'], spins=None, path=None):
        '''
        Options for Atoms are:
        'Total', Plots total DOS
        ['Element_1','Element_2'], List of elements in target structure that need
        to be plotted
        ['Atom_1','Atom_2'], List of atoms as they are numbered in the CONTCAR
        to be plotted
                                                            
        Options for orbitals are 'all' or ['s','p','d','f','s','px','py','pz','dxy',
                                           'dxz','dyz','dz2','dx2-y2','fy3','fz3','fx3',
                                           'fyz2-x2','fzx2-y2','fxz2-y2','fxyz']
        Spins should be a list of spins of the length of the number of files.
        format: [['up','down'],['up','down']]. Leave blank to specify up and down
        for all files.
        '''
        
        self.Atoms = Atoms            

        if path is None:
            path = os.getcwd()
        else:
            path = os.path.normpath(path)
          
        if spins is None:
            spins = []
            for i in range(len(files)):
                spins.append(['up','down'])
                
        
        self.path = path    
        self.orbitals = orbitals
        self.files = files
        self.spins = spins
    
    def open_jdos(self,files):
        all_jdos = {}
        for file in files:
            with open(os.path.join(self.path,file),'r') as f:
                jdos = json.load(f)
            all_jdos[file] = jdos
        
        return all_jdos
    
    def dict_build(self):
        inputs_dict = {}
        for ifile, file in enumerate(self.files):
            inputs_dict[file] = {}
            try:
                for isite, site in enumerate(self.Atoms[ifile]):
                    inputs_dict[file][site] = self.orbitals[ifile][isite]
            except:
                inputs_dict[file] = {'Total':'all'}
        return inputs_dict
            
    def gen_plot_strings(self):
        #This function generates the strings used in visualize() to determine what
        #gets plotted. This loops through the dictionary of inputs built by 
        #dict_build()
        
        spins = ['up','down']
        #TODO
        #add spin customization for individual sites
        
        input_dict = self.dict_build()
        plot_strings = {}
        for file_key in input_dict:
            dos_to_plot = []
            if 'Total' in input_dict[file_key].keys():
                for spin in spins:
                    dos_to_plot.append(tuple(['Total',spin]))
                plot_strings[file_key] = dos_to_plot
            
            else:
                for site_key in input_dict[file_key]:
                    for orbit in input_dict[file_key][site_key]:
                        for spin in spins:
                            dos_to_plot.append(tuple([site_key,orbit,spin]))
                plot_strings[file_key] = dos_to_plot
        return plot_strings
    
    def visualize(self,plot_strings,all_jdos,smear=0.1,zero_width=0.02,elim=[-9, 1],
                 dens_range=[-75,75], dos_lines=[], colors=None, save=False, savename=None):
        
        nplots = len(plot_strings.keys())
        fig, axes = plt.subplots(1, nplots, figsize = (15, 5), sharey = 'row',
                                 dpi = 500 if save else 100)
        
        if colors is None:
            colors = ['#0D7B98','#0D7B98','#15CBAD','#15CBAD','#FDE869','#FDE869',
                      '#FC7968','#FC7968','#FA4A59','#FA4A59']
        plot_colors = []
        if savename is None:
            savename = 'Dos_Plot_{self.file[0])'

        for ifile,file in enumerate(plot_strings.keys()):
            ##Plotting initialization##
            plotter = DosPlotter(sigma = 0.3, zero_at_efermi=True)  # energies are already zeroed!
            
            if nplots == 1:
                ax = axes
            else: 
                ax = axes[ifile]
            jdos = all_jdos[file]
            efermi = all_jdos[file]['Efermi']
            names = []
            for ipdos, pdos in enumerate(plot_strings[file]): #for entries in list of plot string 
            #tuples stored in plot_strings dictionary under key [file].
                dens = []
                spin = pdos[-1]
                if pdos[0] == 'Total':
                    dens = jdos[spin]['Total']
                    energy, density = fill_zeros(jdos[spin]['Energy'], dens, zero_width)
                    energy = [e + efermi for e in energy]
                    max_dens = max(dens)
                    dens_range = [-(0.2*max_dens),(0.2*max_dens)]
                    name = f'Total {spin}'
                else:
                    if len(pdos[0].split('_')) == 1:
                        dos_type = 'element'
                    else:
                        dos_type = 'atom'
                    if dos_type == 'atom': ##################Atom
                        if pdos[1] == 'all': ######################all
                            accepted_orbitals = ['s','p','d','f']
                            initial = False
                            for orbit in jdos[spin][pdos[0]]:
                                if initial is False:
                                    dens_add = np.zeros(len(jdos[spin][pdos[0]][orbit]))
                                    initial = True
                                if orbit in accepted_orbitals:
                                    dens_add += np.array(jdos[spin][pdos[0]][orbit])
                            dens = dens_add
                            energy, density = fill_zeros(jdos[spin]['Energy'], dens, zero_width)
                            energy = [e + efermi for e in energy]
                            name = f'{pdos[0]} {pdos[1]} {spin}'
                        else:                           ######## individual orbital
                            dens = jdos[spin][pdos[0]][pdos[1]] # Atom_# + orbital
                            energy, density = fill_zeros(jdos[spin]['Energy'], dens, zero_width)
                            energy = [e + efermi for e in energy]
                            max_dens = max(dens)
                            dens_range = [-(1.3*max_dens),(1.3*max_dens)]
                            name = f'{pdos[0]} {pdos[1]} {spin}'
                    if dos_type == 'element': ###################Element
                        if pdos[1] == 'all':  ######################all
                            element = pdos[0]
                            initial = False
                            accepted_orbitals = ['s','p','d','f']
                            for atom in jdos[spin]:
                                if atom.split('_')[0] == element:
                                    for orbit in jdos[spin][atom]:
                                        if initial is False:
                                            dens_add = np.zeros(len(jdos[spin][atom][orbit]))
                                            initial = True
                                        if orbit in accepted_orbitals:
                                            dens_add += np.array(jdos[spin][atom][orbit])
                            dens = dens_add
                            energy, density = fill_zeros(jdos[spin]['Energy'], dens, zero_width)
                            energy = [e + efermi for e in energy]
                            name = f'{pdos[0]} {pdos[1]} {spin}'
                        else:     ################################ individial orbital
                            element = pdos[0]
                            initial = False
                            for atom in jdos[spin]:
                                if atom.split('_')[0] == element:
                                    if initial is False:
                                        dens_add = np.zeros(len(jdos[spin][atom][pdos[1]]))
                                        initial = True
                                    dens_add += np.array(jdos[spin][atom][pdos[1]])
                                else: continue
                            dens = dens_add
                            energy, density = fill_zeros(jdos[spin]['Energy'], dens, zero_width)
                            energy = [e + efermi for e in energy]
                            name = f'{pdos[0]} {pdos[1]} {spin}'
                            
                pmg_jdos = dos.Dos(efermi, #if center_fermi else 0.0, 
                                       energy,
                                       {Spin.up if spin == 'up' else Spin.down: density})

                plotter.add_dos(name, pmg_jdos)
                names.append(name)
                
                plot_colors.append(colors[ipdos])

            plot = get_plot(plotter, energy_lim=elim, colors=plot_colors, alpha = 0.3,
                            density_lim = dens_range,# if axid in [0,1,2,3,4] else None,
                            normalize_density=False, ax = ax, 
                            lloc = (-0.1, 1.05), #if axid == 3 else (1.0, 1.0), 
                            lcol = 1, lframe = True,
                            ylabel = True if ifile == 0 else False, density_ticks = False, 
                           # mark_fermi = efermi if axid in [0,1,2,3] else None,
                            mark_fermi = True, fill_to_efermi = False, pdos_label = True,
                            dos_lines = dos_lines, show_legend = True
                            )
        plt.show()
        
        if save:
            plt.savefig(savename)
            plt.close()
        
    def plot_dos(self,smear=0.1,save=False,savename=None,zero_width=0.02,elim=[-8,2],
                 dos_lines=[],colors=None):
        all_jdos = self.open_jdos(self.files)
        inputs_dict = self.dict_build()
        plot_strings = self.gen_plot_strings()
        self.visualize(plot_strings,all_jdos,smear=smear,save=save,
                       savename=savename,zero_width=zero_width,elim=elim,
                       dos_lines=dos_lines,colors=colors)
        
                    
dos_object = DOS_visualize(files=['jpdos.json'], 
                           Atoms=[['Total'],['Ir','Ir_2'],['Ir_20','Ir_21']],
                            orbitals=[['all'],[['s','d'],['all']],[['all'],['all']]],
                            path='C:/Users/coopy/OneDrive - UCB-O365/Desktop/temp')
                           
# dos_object = DOS_visualize(files=['jpdos.json'], Atoms=[['Total']],
#                             orbitals=[['all']],
#                             path='C:/Users/coopy/OneDrive - UCB-O365/Desktop/temp')

# test_dict = dos_object.dict_build()
# plot_dict = dos_object.gen_plot_strings()
# all_jdos = dos_object.open_jdos(['jpdos.json','jpdos _2.json'])
# dos_object.visualize(plot_dict,all_jdos)

dos_object.plot_dos()

