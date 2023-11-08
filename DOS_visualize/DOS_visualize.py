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

import sys #TODO remove from github
sys.path.append("C:\\Users\\coopy\\OneDrive - UCB-O365\\github\\JDFT_tools\\DOS_visualize")

from dos_helper import set_rc_params, get_plot, fill_zeros
set_rc_params()

class DOS_visualize:
    def __init__(self, files, Atoms=[], orbitals=[], spins=None, path=None):
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
                spins.append(['spinup','spindn'])
                
        
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
            for isite, site in enumerate(self.Atoms[ifile]):
                inputs_dict[file][site] = self.orbitals[ifile][isite]
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
    
    def visualize(self,plot_strings,all_jdos,smear=0.1,zero_width=0.02,elim=[-9, 1], dens_lim=None,
                 dens_range=[-75,75], dos_lines=[], input_colors=None, save=False, savename=None,
                 legend = True):
        nplots = len(plot_strings.keys())
        fig, axes = plt.subplots(1, nplots, figsize = (10, 5), sharey = 'row',
                                 dpi = 500 if save else 100)
        plot_colors = []
        if savename is None:
            savename = 'Dos_Plot_{self.file[0])'

        for ifile,file in enumerate(plot_strings.keys()):
            
            colors = []
            if input_colors is None:
                colors = ['#0D7B98','#0D7B98','#15CBAD','#15CBAD','#FDE869','#FDE869',
                          '#FC7968','#FC7968','#FA4A59','#FA4A59']
            elif type(input_colors) is list:
                colors = input_colors[ifile]
                        
            ##Plotting initialization##
            plotter = DosPlotter(sigma = smear, zero_at_efermi=True)  # energies are already zeroed!
            
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
                    # dens_range = [-1000,1000]
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
                            dens_range = [-(0.001*max_dens),(0.001*max_dens)]
                            # dens_range = [-.05,0.05]
                            name = f'{pdos[0]} {pdos[1]} {spin}'
                    if dos_type == 'element': ###################Element #TODO fix element
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
                                elif atom in ['Energy','Total']:
                                    continue
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
                            density_lim = dens_range if dens_lim is None else dens_lim[ifile],
                            # density_lim = [-10,10],
                            normalize_density=False, ax = ax, 
                            lloc = (-0.1, 1.05),  
                            lcol = 1, lframe = True,
                            ylabel = True if ifile == 0 else False, density_ticks = False, 
                            mark_fermi = True, fill_to_efermi = True, pdos_label = True,
                            dos_lines = dos_lines, show_legend = True if legend else False
                            )
        
        if save:
            plt.savefig(os.path.join(self.path,savename))
            # plt.show()
            # plt.close()
        elif save == False:
            plt.show()
            plt.close()
        
    def plot_dos(self,smear=0.1,save=False,savename=None,zero_width=0.02,elim=[-8,2],
                 dos_lines=[],colors=None,dens_lim=None,legend=True):
        all_jdos = self.open_jdos(self.files)
        inputs_dict = self.dict_build()
        plot_strings = self.gen_plot_strings()
        self.visualize(plot_strings,all_jdos,smear=smear,save=save,
                       savename=savename,zero_width=zero_width,elim=elim,
                       dos_lines=dos_lines,input_colors=colors,dens_lim=dens_lim,legend=legend,
                       )
        


if __name__ == '__main__':   

    path = 'C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\Binaries\\DOS'
    ads_files = ["jpdos_0V.json", "jpdos_5V.json"]
    surf_files = ["jpdos_surf_0V.json", "jpdos_surf_5V.json"]
    endon_files = ["jpdos_endon_0V.json", "jpdos_endon_5V.json"]
    NH_file = ["jpdos_NH_0V.json"]
    TaB_file = ["jpdos_TaB_surf_0V.json", "jpdos_TaB_N2_0V.json","jpdos_TaB_N2H_0V.json","jpdos_TaB_NH_0V.json", "jpdos_TaB_NH3_0V.json"]
    ############ Place test code here ###############        
    # dos_object = DOS_visualize(files=['jpdos.json'], 
    #                            Atoms=[['Total'],['Ir','Ir_2'],['Ir_20','Ir_21']],
    #                             orbitals=[['all'],[['s','d'],['all']],[['all'],['all']]],
    #                             path='C:/Users/coopy/OneDrive - UCB-O365/Desktop/temp')
    
    # Cr_dos = DOS_visualize(files=['jpdos_N2.json','jpdos_Cr_surf.json','jpdos_Cr_N2_3.json'],
    #                             Atoms=[['Total'],['Total'],['Cr_5','N_1']],
    #                             orbitals=[['p'],['all'],[['d'],[ 'p']]],
    #                             path='C:/Users/coopy/OneDrive - UCB-O365/Desktop/temp')

    # V_dos = DOS_visualize(files=['jpdos_V_N2.json'], 
    #                         Atoms=[['V_23', 'V_24', 'V_25', 'V_26', 'N_1', 'N_2']],
    #                         orbitals = [[['d'],['d'],['d'],['d'],['p'],['p']]],
    #                         path=path
    # )

    # V_N2_dos = DOS_visualize(files=['jpdos_V_N2.json'], 
    #                     Atoms=[[ 'N_1', 'N_2']],
    #                     orbitals = [[['p'],['p']]],
    #                     path=path
    # )

    # Ti_sideon_dos = DOS_visualize(files=ads_files,
    #                             Atoms=[['N_2', 'N_1','Ti_7','Ti_8'], ['N_2', 'N_1','Ti_7','Ti_8']],
    #                             orbitals=[[['p'], ['p'],['d'],['d']], [['p'], ['p'],['d'],['d']]],
    #                             path=path
    # )

    # Ti_NH_dos = DOS_visualize(files=['jpdos_surf_0V.json','jpdos_NH_0V.json','jpdos_NH_mol_0V.json'],
    #                            Atoms=[['Ti_7','Ti_8'], ['Ti_7','Ti_8','N_1','H_1'], ['N']],
    #                            orbitals=[[['d'],['d']],[['d'],['d'],['p'],['s']], [['all']]],
    #                             path=path
    # )

    # TaB_dos = DOS_visualize(files=TaB_file,
    #                         Atoms=[['Ta_4','Ta_12', 'Ta_8', 'Ta_16'],['Ta_4','Ta_12', 'Ta_8', 'Ta_16','N'],['Ta_4','Ta_12', 'Ta_8', 'Ta_16','N','H'],
    #                                ['Ta_4','Ta_12', 'Ta_8', 'Ta_16','N','H'],['Ta_4','Ta_12', 'Ta_8', 'Ta_16','N','H']],
    #                         orbitals=[[['d'],['d'],['d'],['d']], [['d'],['d'],['d'],['d'],['p']], [['d'],['d'],['d'],['d'],['p'],['s']],
    #                                   [['d'],['d'],['d'],['d'],['p'],['s']],[['d'],['d'],['d'],['d'],['p'],['s']]],
    #                         path=path
    # )
    # Ti_8_dos = DOS_visualize(files=ads_files,
    #                             Atoms=[['Ti_8'], ['Ti_8']],
    #                             orbitals=[[['dxy','dxz','dyz','dz2','dx2-y2']], [['dxy','dxz','dyz','dz2','dx2-y2']]],
    #                             path=path)

    # Ti_7_dos = DOS_visualize(files=ads_files,
    #                             Atoms=[['Ti_7'], ['Ti_7']],
    #                             orbitals=[[['dxy','dxz','dyz','dz2','dx2-y2']], [['dxy','dxz','dyz','dz2','dx2-y2']]],
    #                             path=path)

    # ads_dos = DOS_visualize(files=ads_files,
    #                             Atoms=[['Ti_8'], ['Ti_8']],
    #                             orbitals=[[['dxy','dxz','dyz','dz2','dx2-y2'],['dxy','dxz','dyz','dz2','dx2-y2']], [['dxy','dxz','dyz','dz2','dx2-y2'],['dxy','dxz','dyz','dz2','dx2-y2']]],
    #                             path=path)
                                
    # surf_dos = DOS_visualize(files=surf_files,
    #                          Atoms=[['Ti','C'], ['Ti','Ti_7','Ti_8']],
    #                             orbitals=[[['all'],['all'], ['all']], [['all'], ['all'], ['all']]],
    #                             path=path)
                             
    # endon_dos = DOS_visualize(files=endon_files,
    #                             Atoms=[['N_2', 'N_1','Ti_7','Ti_8'], ['N_2', 'N_1','Ti_7','Ti_8']],
    #                             orbitals=[[['p'], ['p'],['d'],['d']], [['p'], ['p'],['d'],['d']]],
    #                             path=path)


    # HfP_dos  = DOS_visualize(files=["jpdos_HfP_NH3_0V.json"],
    #                             Atoms=[['N_2', 'N_1','Ti_7','Ti_8'] ],
    #                             orbitals=[[['p'], ['p'],['d'],['d']], [['p'], ['p'],['d'],['d']]],
    #                             path=path)

    TiC_dos = DOS_visualize(files=["jpdos_TiC_N_0V.json","jpdos_TiC_NH_0V.json","jpdos_TiC_NH2_0V.json","jpdos_TiC_NH3_0V.json"],
                            Atoms=[['Ti_15','Ti_16','Ti_7','Ti_8','N_1'], ['Ti_15','Ti_16','Ti_7','Ti_8','N_1','H_1'],['Ti_15','Ti_16','Ti_7','Ti_8','N_1','H_1','H_2'],['Ti_15','Ti_16','Ti_7','Ti_8','N_1','H_1','H_2','H_3']],
                            orbitals=[[['d'],['d'],['d'],['d'],['p','s']], [['d'],['d'],['d'],['d'],['p','s'],['s']], [['d'],['d'],['d'],['d'],['p','s'],['s'],['s']], [['d'],['d'],['d'],['d'],['p','s'],['s'],['s'],['s']]],
                            path=path)
    # TiC_surf_dos = DOS_visualize(files=["jpdos_surf_0V.json"],
    #                             Atoms=[[]]
    
    
    # )

# Plotting
    
    # Ti_sideon_dos.plot_dos(save=True,
    #                         savename=os.path.join(path,'sideon_pdos'),
    #                         dens_lim=[[-100,100],[-100,100]],
    #                         colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                         ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
    #                         elim=[-8,1.5],smear=0.1
    #                         # legend=False
    # )


                               
                               




    '''Ti_7_dos.plot_dos(save=True, 
                      savename=os.path.join(path,'metal_only_adsorbed_pdos_Ti7'),
                        dens_lim=[[-100,100],[-100,100]],
                        colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
                            ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
                            elim=[-8,1.5],smear=0.1
                        )
    Ti_8_dos.plot_dos(save=True, 
                      savename=os.path.join(path,'metal_only_adsorbed_pdos_Ti8'),
                        dens_lim=[[-100,100],[-100,100]],
                        colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
                            ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
                            elim=[-8,1.5],smear=0.1
                        )
    surf_dos.plot_dos(save=True,
                        savename=os.path.join(path,'metal_only_surface_pdos'),
                            dens_lim=[[-600,600],[-600,600]],
                            colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
                            ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
                            elim=[-8,1.5],smear=0.1
                            # legend=False
                            )'''
    
    # endon_dos.plot_dos(save=True, 
    #                   savename=os.path.join(path,'metal_only_endon_pdos'),
    #                     dens_lim=[[-100,100],[-100,100]],
    #                     colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                         ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
    #                     elim=[-8,1.5],smear=0.1
    #                     # legend=False
    #                     )

    # V_dos.plot_dos(save=True,
    #                 savename=os.path.join(path,'V_pdos'),
    #                 dens_lim=[[-100,100],[-100,100]],
    #                     colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                         ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
    #                     elim=[-8,1.5],smear=0.1
    #                     # legend=False
    # )

    # V_N2_dos.plot_dos(save=True,
    #                   savename=os.path.join(path,'V_N2_pdos'),
    #                 dens_lim=[[-50,50],[-50,50]],
    #                     colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                         ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
    #                     elim=[-8,1.5],smear=0.1
    #                     # legend=False
    # )

    # Ti_NH_dos.plot_dos(save=True,
    #                     savename=os.path.join(path,'Ti_NH_pdos'),
    #                     dens_lim=[[-50,50],[-50,50]],
    #                         colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                             ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                             ['#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']],
    #                         elim=[-8,1.5],smear=0.1
    #                         # legend=False
    #     )

    # TaB_dos.plot_dos(save=True,
    #                  savename=os.path.join(path,'TaB_pdos'),
    #                  dens_lim=[[-50,50],[-50,50]],
    #                  colors=[['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                          ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                          ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                          ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522'],
    #                          ['#FF1F5B','#FF1F5B','#00CD6C','#00CD6C','#009ADE','#009ADE','#AF58BA','#AF58BA','#FFC61E','#FFC61E','#F28522','#F28522']]   ,
    #                          elim=[-8,1.5],smear=0.1
    # )    
    TiC_dos.plot_dos(save=True,
                    savename=os.path.join(path,'TiC_NHx_pdos_renumber'),
                    dens_lim=[[-50,50],[-50,50],[-50,50],[-50,50]],
                    colors=[['#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#FFC61E','#FFC61E','#FFC61E','#FFC61E','#F28522','#F28522','#F28522','#F28522'],
                                ['#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#FFC61E','#FFC61E','#F28522','#FFC61E','#FFC61E','#F28522','#F28522','#F28522'],
                                ['#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#FFC61E','#FFC61E','#F28522','#FFC61E','#FFC61E','#F28522','#F28522','#F28522','#F28522','#F28522'],
                                ['#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#00CD6C','#FFC61E','#FFC61E','#F28522','#FFC61E','#FFC61E','#F28522','#F28522','#F28522','#F28522','#F28522','#F28522','#F28522']],
                    elim=[-10,5],smear=0.1
)

                 