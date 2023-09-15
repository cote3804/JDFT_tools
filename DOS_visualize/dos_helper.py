# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:34:17 2021

@author: NSing
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcdefaults()
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.plotting import pretty_plot

def set_rc_params():
    """
    Args:
        
    Returns:
        dictionary of settings for mpl.rcParams
    """
    params = {'axes.linewidth' : 1.5,'axes.unicode_minus' : False,
              'figure.dpi' : 100,
              'font.size' : 16,'font.family': 'sans-serif','font.sans-serif': 'Verdana',
              'legend.frameon' : False,'legend.handletextpad' : 0.2,
              'legend.handlelength' : 0.6,'legend.fontsize' : 12,
              'legend.columnspacing': 0.8,
              'mathtext.default' : 'regular','savefig.bbox' : 'tight',
              'xtick.labelsize' : 16,'ytick.labelsize' : 16,
              'xtick.major.size' : 6,'ytick.major.size' : 6,
              'xtick.major.width' : 1.5,'ytick.major.width' : 1.5,
              'xtick.top' : False,'xtick.bottom' : True,'ytick.right' : True,'ytick.left' : True,
              'xtick.direction': 'out','ytick.direction': 'out','axes.edgecolor' : 'black'}
    for p in params:
        mpl.rcParams[p] = params[p]
    return params

def get_plot(plotter, energy_lim=[-5, 5], density_lim=None, flip_axes = True, colors = None,
             normalize_density = True, fill = True, alpha = 1, ax = None, 
             lloc = None, lcol = 1, lframe = True, ylabel = True, density_ticks = True,
             mark_fermi = True, fill_to_efermi = True, pdos_label = True,
             dos_lines = [], show_legend = True):
    """
    Taken from pymatgen.electronic_structure.plotter, DosPlotter.get_plot()
    Needed to rewrite to allow for axes flipping
    """
    if colors is None:
        import palettable
        ncolors = max(3, len(plotter._doses))
        ncolors = min(9, ncolors)
        # pylint: disable=E1101
        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
    else:
        ncolors = len(colors)
        
    y = None
    max_density = 0
    min_density = 0
    alldensities = []
    allenergies = []
    added_keys = []
    
    if ax is None:
        plot = pretty_plot(8, 12)
        plfill = plot.fill_between
        plfillx = plot.fill_betweenx
        plplot = plot.plot
        plxm = plot.xlim
        plym = plot.ylim
    else:
        plfill = ax.fill_between
        plfillx = ax.fill_betweenx
        plplot = ax.plot
        plxm = ax.set_xlim
        plym = ax.set_ylim

    for key, idos in plotter._doses.items():
        energies = idos["energies"]
        densities = idos["densities"]
        allenergies.append(energies)
        alldensities.append(densities)

    keys = list(plotter._doses.keys())
#    print(keys)
    keys.reverse()
    colors.reverse()
    if type(alpha) == list:
        alpha.reverse()
#    print(keys)
    alldensities.reverse()
    allenergies.reverse()
    allpts = []
    for i, key in enumerate(keys):
        x = []
        y = []
        for spin in [Spin.up, Spin.down]:
            if spin in alldensities[i]:
                densities = list(int(spin) * alldensities[i][spin])
                energies = list(allenergies[i])
                if spin == Spin.down:
                    energies.reverse()
                    densities.reverse()
                x.extend(energies)
                y.extend(densities)
        allpts.extend(list(zip(x, y)))
        
        maxy = max([yi for ii,yi in enumerate(y) if x[ii] <= max(energy_lim) and x[ii] >= min(energy_lim)])
        miny = min([yi for ii,yi in enumerate(y) if x[ii] <= max(energy_lim) and x[ii] >= min(energy_lim)])
        if normalize_density:
            y = [yi / maxy for yi in y]
            density_lim = [0, 1.1]
        if maxy > max_density:
            max_density = maxy
        if miny < min_density:
            min_density = miny
            
        if flip_axes:
            tmp = x
            x = y
            y = tmp
        
        red_key = str(key).replace(' up','').replace(' down','')
        if red_key != 'Total':
            red_key = red_key.split()[0] + '('+red_key.split()[1]+')'
        if red_key in added_keys:
            label = '__No_Label__' 
        else:
            label = red_key#.split()[0] + ' ('+red_key.split()[1]+')'
            added_keys.append(label)
        
        
        if fill and flip_axes:
            # x = densities, y = energies
            if not fill_to_efermi:
                plfillx(y, np.zeros_like(y), x, color=colors[i % ncolors], label=label,
                              alpha = alpha[i % ncolors] if type(alpha) == list else alpha)
            else:
                eflim = 0.0 if plotter.zero_at_efermi else plotter._doses[key]["efermi"]
                xx = [xi for ix,xi in enumerate(x) if y[ix] < eflim]  # densities
                yy = [yi for yi in y if yi < eflim]  # energies
                plfillx(yy, np.zeros_like(yy), xx, 
                        color=colors[i % ncolors], label=label,
                        alpha = alpha[i % ncolors] if type(alpha) == list else alpha)
                plplot(x, y, color=colors[i % ncolors], #label=label,
                       alpha = alpha[i % ncolors] if type(alpha) == list else alpha)
                
        elif fill:
            plfill(x, np.zeros_like(x), y, color=colors[i % ncolors], label=label,
                              alpha = alpha[i % ncolors] if type(alpha) == list else alpha)
        else:
            plplot(x, y, color=colors[i % ncolors], label=label, linewidth=1,
                      alpha = alpha[i % ncolors] if type(alpha) == list else alpha)
        
#        if not plotter.zero_at_efermi:
#            ylim = plym()
#            plplot([plotter._doses[key]["efermi"], plotter._doses[key]["efermi"]],ylim,
#                      color=colors[i % ncolors],linestyle="--",linewidth=1,)

    if not flip_axes:
        if energy_lim:
            plxm(energy_lim)
        if density_lim:
            plym(density_lim)
        else:
            xlim = plxm()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plym((min(relevanty), max(relevanty)))
    else:
        if density_lim:
            plxm(density_lim)
        else:
            plxm((1.1 * min_density, max_density * 1.1))
            print('Auto density range')
        if energy_lim:
            plym(energy_lim)

#    if plotter.zero_at_efermi:
#        if not flip_axes:
#            ylim = plym()
#            plplot([0, 0], ylim, "k--", linewidth=2)
#        else:
#            plplot(plxm(), [0, 0], "k--", linewidth=2)
    
    elabel = 'Energy (eV)' if not plotter.zero_at_efermi else 'Energy - $E_{Fermi}$ (eV)'
    pdoslabel = 'pDOS (a.u.)' if pdos_label else "Density of states (a.u.)"
    
    if ax is not None:
        if flip_axes:
            if ylabel: ax.set_ylabel(elabel)
            ax.set_xlabel(pdoslabel)            
        else:
            ax.set_xlabel(elabel)
            if ylabel: ax.set_ylabel(pdoslabel)
#        if mark_fermi is not None:
#            ax.axhline(y=mark_fermi, color="k", linestyle="--", linewidth=1)
#        else:
#            ax.axhline(y=0, color="k", linestyle="--", linewidth=1)
        
        if mark_fermi:
            if plotter.zero_at_efermi:
                ax.axhline(y=0, color="k", linestyle="--", linewidth=1)
            else:
                ax.axhline(y=plotter._doses[key]["efermi"], color="k", linestyle="--", linewidth=1)
        
        ax.axvline(x=0, color="k", linestyle="--", linewidth=1)
        
        if not density_ticks:
            ax.set_xticklabels([])
        
        if len(dos_lines) > 0:
            for line in dos_lines:
                if flip_axes:
                    ax.axhline(y=line[0], color=line[1], linestyle="-", linewidth=1)
                else:
                    ax.axvline(x=line[0], color=line[1], linestyle="-", linewidth=1)
        
        if show_legend:
            if lloc is None:
                ax.legend()
            else:
                mpl.rcParams['legend.frameon'] = lframe
                ax.legend(loc='upper center', bbox_to_anchor=lloc,
                          ncol=lcol, fancybox=lframe, shadow=lframe)
        return ax
    else:
        if flip_axes:
            if ylabel: plot.ylabel(elabel)
            plot.xlabel(pdoslabel)           
        else:
            plot.xlabel(elabel)
            if ylabel: plot.ylabel(pdoslabel)
#        if mark_fermi is not None:
#            plot.axhline(y=mark_fermi, color="k", linestyle="--", linewidth=1)
#        else:
#            plot.axhline(y=0, color="k", linestyle="--", linewidth=1)
        if mark_fermi:
            if plotter.zero_at_efermi:
                plot.axhline(y=0, color="k", linestyle="--", linewidth=1)
            else:
                plot.axhline(y=plotter._doses[key]["efermi"], color="k", linestyle="--", linewidth=1)
        plot.axvline(x=0, color="k", linestyle="--", linewidth=1)
        plot.legend()
        return plot
#    leg = plot.gca().get_legend()
#    ltext = leg.get_texts()  
#    plot.setp(ltext, fontsize=30)
#    plot.tight_layout()

def fill_zeros(energy, density, interval = 0.1, max_itt = 10000):
    new_energy, new_density = [], []
    for ii, en in enumerate(energy):
        den = density[ii]
        if ii == 0:
            new_energy.append(en)
            new_density.append(den)
            continue
        if en - new_energy[-1] < interval:
            new_energy.append(en)
            new_density.append(den)
        else:
            last_en = new_energy[-1]
            for i in range(max_itt):
                # continuously add new energies at 'interval' until reaching next energy
                assert en - new_energy[-1] > 0.0, ('METAERROR: Energy difference is negative: '+
                                                   str(en)+' - '+str(new_energy[-1])+', '+str(last_en))
                if en - new_energy[-1] <= interval+0.001:
                    new_energy.append(en)
                    new_density.append(den)
                    break
                new_energy.append(last_en + (i+1)*interval)
                new_density.append(0.0)
    return new_energy, new_density