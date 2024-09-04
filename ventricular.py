#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human ventricular
# models.
#
import os

import matplotlib
import matplotlib.pyplot as plt
import myokit
import myokit.lib.plots as mp
import numpy as np

import shared


# Update matplotlib styles
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'

# Current colors
cmap = matplotlib.colormaps['tab20']
current_colours = dict(shared.current_colours)
del(current_colours['I_f'])
del(current_colours['I_Kur'])
del(current_colours['I_CaT'])
del(current_colours['I_SK'])

# Human atrial models
model_names = {
    'bartolucci': 'bartolucci-2020.mmt',
    'carro': 'carro-2011.mmt',
    'cipa': 'ohara-cipa-v1-2017.mmt',
    'fink': 'fink-2008.mmt',
    'grandi': 'grandi-2010.mmt',
    'iyer': 'iyer-2004.mmt',
    'ohara': 'ohara-2011.mmt',
    'priebe': 'priebe-1998.mmt',
    'tnnp': 'tentusscher-2004.mmt',
    'tomek': 'tomek-2020.mmt',
    'tp': 'tentusscher-2006.mmt',
}


model_modes = {
    'bartolucci': {'cell.mode': 1},
    'carro': {'mode.epi': 1},
    'grandi': {'mode.epi': 1},
    'tnnp': {'cell.type': 1},
    'tp': {'cell.type': 1},
    'ohara': {'cell.mode': 1},
    'cipa': {'cell.mode': 1},
    'tomek': {'cell.mode': 1},
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'priebe' in name:
        currents = {
            'I_NaCa': 'inaca.i_NaCa',
            'I_to': 'ito.i_to',
            'I_Ks': 'iks.i_Ks',
            'I_Kr': 'ikr.i_Kr',
            'I_K1': 'ik1.i_K1',
            'I_NaK': 'inak.i_NaK',
            'I_CaL': 'ica.i_Ca',
            'I_Ca,B': 'icab.i_b_Ca',
            'I_Na,B': 'inab.i_b_Na',
            'I_Na': 'ina.i_Na',
        }
    elif 'iyer' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito1',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_Na': 'ina.INa',
        }
    elif 'grandi' in name:
        currents = {
            'I_Cl,B': 'iclb.IClB',
            'I_ClCa': 'iclca.IClCa',
            'I_to': 'ito.Ito',
            'I_Kp': 'ikp.IKp',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_Na': 'ina.INa',
        }
    elif 'tusscher-2004' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito',
            'I_Kp': 'ipk.IpK',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_CaL': 'ical.ICaL',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'tusscher-2006' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kp': 'ipk.IpK',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_K1': 'ik1.IK1',
            'I_NaCa': 'inaca.INaCa',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_CaL': 'ical.ICaL',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'ohara-2011' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kp': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_tot',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inacass.INaCa_tot',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'dutta-2017' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kp': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_tot',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inacass.INaCa_tot',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'tomek-2020' in name:
        currents = {
            'I_Cl,B': 'iclb.IClb',
            'I_ClCa': 'iclca.IClCa',
            'I_to': 'ito.Ito',
            'I_Kp': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_K,ATP': 'ikatp.IKatp',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_total',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'fink-2008' in name:
        currents = {
            'I_to': 'ito.i_to',
            'I_Kp': 'ipk.i_p_K',
            'I_Ks': 'iks.i_Ks',
            'I_Kr': 'ikr.i_Kr',
            'I_K1': 'ik1.i_K1',
            'I_NaCa': 'inaca.i_NaCa',
            'I_NaK': 'inak.i_NaK',
            'I_Ca,P': 'ipca.i_p_Ca',
            'I_CaL': 'ical.i_CaL',
            'I_Ca,B': 'icab.i_b_Ca',
            'I_Na,B': 'inab.i_b_Na',
            'I_Na': 'ina.i_Na',
        }
    elif 'carro' in name:
        currents = {
            'I_Cl,B': 'iclb.IClB',
            'I_ClCa': 'iclca.IClCa',
            'I_to': 'ito.Ito',
            'I_Kp': 'ikp.IKp',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_Na': 'ina.INa',
        }
    elif 'bartolucci-2020' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kp': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inacass.INaCa_total',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    else:
        currents = shared.guess_currents(model)
        print('\n'.join(currents))
        print(len(currents))
        raise NotImplementedError('Unknown model: ' + model.name())

    if colours:
        colours = [cmap(current_colours[x]) for x in currents.keys()]
        currents = list(currents.values())
        return currents, colours
    return list(currents.values())


# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = {}
for name, fname in model_names.items():
    print(f'Preparing {name}...')
    pre_pace = True
    if 'priebe' in name:
        # 2024-09-03 Priebe runs into numerical issues when prepacing
        pre_pace = False
    model = myokit.load_model(os.path.join('models', 'c', fname))
    for k, v in model_modes.get(name, {}).items():
        print(f'  Setting {k} to {v}')
        model.get(k).set_rhs(v)
    shared.prepare_model(model, protocol, current_variables(model), pre_pace)
    models[name] = model
    model.labelx('g_Kr')
print('Finished preparation.\nPreparing plots')


# Maximum time to show in plots
tmax = 800


def text(ax, x, y, t, c='w'):
    ax.text(x, y, t, color=c, transform=ax.transAxes, fontweight='bold',
            horizontalalignment='right', verticalalignment='center')


def plot(grid, code, ylabel='Relative contribution', legend=False):

    model = models[code]
    print(f'+ {model.meta["display_name"]}')
    currents, colours = current_variables(model, True)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    d = s.run(tmax)
    s.reset()
    g = model.labelx('g_Kr')
    s.set_constant(g.qname(), 0.3 * g.eval())
    e = s.run(tmax)

    # V
    v = model.labelx('membrane_potential')
    gr = grid.subgridspec(4, 1, hspace=0)
    ax = fig.add_subplot(gr[0, 0])
    ax.set_ylabel('V (mV)')
    ax.set_title(model.meta['display_name'])
    ax.set_xticklabels([])
    ax.plot(d.time(), d[v], 'k', label='baseline')
    ax.plot(e.time(), e[v], 'k--', label='30% IKr')
    ax.set_xlim(0, tmax)
    ax.set_ylim(-95, 45)
    ax.set_yticks([-80, -40, 0, 40])
    if legend:
        ax.legend(loc='upper right', frameon=False)

    # Total current
    #k = model.labelx('cellular_current').qname()
    #ax2 = ax.twinx()
    #ax2.set_ylim(-0.2, 1.6)
    #ax2.plot(d.time(), d[k], 'r')

    # Contributions
    ax = fig.add_subplot(gr[1:, 0])
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, tmax)
    ax.set_ylim(-1.02, 1.02)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels(['-1', '-0.5', '0', '0.5', '1'])
    ax.yaxis.get_majorticklabels()[-1].set_verticalalignment('top')
    ax.yaxis.get_majorticklabels()[0].set_verticalalignment('bottom')

    mp.cumulative_current(d, currents, ax, colors=colours, normalize=True)


# Create figure
fig = plt.figure(figsize=(9, 12.5))
fig.subplots_adjust(0.067, 0.035, 0.98, 0.98, hspace=0.35, wspace=0.25)
grid = fig.add_gridspec(4, 3)

plot(grid[0, 0], 'priebe', legend=True)
plot(grid[0, 1], 'iyer', ylabel=None)
plot(grid[1, 0], 'tnnp')
plot(grid[1, 1], 'tp', ylabel=None)
plot(grid[1, 2], 'fink', ylabel=None)
plot(grid[2, 0], 'grandi')
plot(grid[2, 1], 'carro', ylabel=None)
plot(grid[3, 0], 'ohara')
#plot(grid[3, 1], 'cipa', ylabel=None)
plot(grid[3, 1], 'tomek', ylabel=None)
plot(grid[3, 2], 'bartolucci', ylabel=None)

# Legend
ax = fig.add_subplot(grid[0, 2])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
labels = [shared.current_names[x] for x in current_colours]
ax.legend(lines, labels, loc=(-0.13, 0.05), ncol=2)
#ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)

# Show / store
plt.savefig('ventricular.png')
plt.savefig('ventricular.pdf')
print('Done')
