#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human atrial models.
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
del(current_colours['I_CaT'])
del(current_colours['I_K,ATP'])
del(current_colours['I_SK'])

# Human atrial models
model_names = {
    'aguilar': 'aguilar-2017.mmt',
    'bai': 'bai-2018.mmt',
    'courtemanche': 'courtemanche-1998.mmt',
    'ellinwood': 'ellinwood-2017.mmt',
    'grandi': 'grandi-2011.mmt',
    'koivumaki': 'koivumaki-2011.mmt',
    'maleckar': 'maleckar-2009.mmt',
    'ni': 'ni-2017.mmt',
    'nygren': 'nygren-1998.mmt',
    'voigt': 'voigt-2013.mmt',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'nygren' in name:
        currents = {
            'I_Kur': 'isus.Isus',
            'I_to': 'it.It',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'icap.ICaP',
            'I_Ca,B': 'ibca.IBCa',
            'I_Na,B': 'ibna.IBNa',
            'I_Na': 'ina.INa',
        }
    elif 'maleckar-' in name:
        currents = {
            'I_Kur': 'ikur.IKur',
            'I_to': 'it.It',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'icap.ICaP',
            'I_Ca,B': 'ibca.IBCa',
            'I_Na,B': 'ibna.IBNa',
            'I_K,ACh': 'ikach.IKACh',
            'I_Na': 'ina.INa',
        }
    elif 'koivumaki' in name:
        currents = {
            'I_Kur': 'ikur.IKur',
            'I_to': 'it.It',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'icap.ICaP',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_f': 'if.If',
            'I_Na': 'ina.INa',
        }
    elif 'courtemanche-1998' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'ib.IbCa',
            'I_Na,B': 'ib.IbNa',
            'I_Na': 'ina.INa',
        }
    elif 'ni-' in name:
        currents = {
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'icap.ICap',
            'I_Ca,B': 'ibca.IbCa',
            'I_Na,B': 'ibna.IbNa',
            'I_Na': 'ina.INa',
        }
    elif 'grandi-2011' in name:
        currents = {
            'I_Cl,B': 'iclb.IClB',
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_ClCa': 'iclca.IClCa',
            'I_Kp': 'ikp.IKp',
            'I_Na': 'ina.INa',
            'I_NaL': 'inal.INaL',
        }
    elif 'voigt' in name:
        currents = {
            'I_Cl,B': 'iclb.IClB',
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_K,ACh': 'ikach.IKACh',
            'I_ClCa': 'iclca.IClCa',
            'I_Kp': 'ikp.IKp',
            'I_Na': 'ina.INa',
            'I_NaL': 'inal.INaL',
        }
    elif 'bai-2018' in name:
        currents = {
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Kp': 'ipk.IpK',
            'I_Na': 'ina.INa',
        }
    elif 'ellinwood' in name:
        currents = {
            'I_Cl,B': 'iclb.IClB',
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'icab.ICaB',
            'I_Na,B': 'inab.INaB',
            'I_K,ACh': 'ikach.IKACh',
            'I_ClCa': 'iclca.IClCa',
            'I_Kp': 'ikp.IKp',
            'I_Na': 'ina.INa',
        }
    elif 'aguilar' in name:
        currents = {
            'I_Kur': 'ikur.IKur',
            'I_to': 'ito.Ito',
            'I_CaL': 'ical.ICaL',
            'I_NaCa': 'inaca.INaCa',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_Ca,B': 'ib.IbCa',
            'I_Na,B': 'ib.IbNa',
            'I_K,ACh': 'ikach.IKACh',
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
    if 'koiv' in name:
        # 2024-09-03 Koivumaki doesn't stabilise, with difference increasing
        # even after 60000 beats.
        pre_pace = False
    model = myokit.load_model(os.path.join('models', 'c', fname))
    shared.prepare_model(model, protocol, current_variables(model), pre_pace)
    models[name] = model
    model.labelx('g_Kur')
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
    g = model.labelx('g_Kur')
    try:
        s.set_constant(g.qname(), 0.5 * g.eval())
    except ValueError:
        print('  Creating second simulation')
        m2 = model.clone()
        m2.labelx('g_Kur').set_rhs(0.5 * g.eval())
        s = myokit.Simulation(m2, protocol)
        s.set_tolerance(1e-8, 1e-8)
    e = s.run(tmax)

    # V
    v = model.labelx('membrane_potential')
    gr = grid.subgridspec(4, 1, hspace=0)
    ax = fig.add_subplot(gr[0, 0])
    ax.set_title(model.meta['display_name'])
    ax.set_xticklabels([])
    ax.plot(d.time(), d[v], 'k', label='Baseline')
    ax.plot(e.time(), e[v], 'k--', label='50% IKur')
    ax.set_xlim(0, tmax)
    ax.set_ylim(-95, 45)
    ax.set_yticks([-80, -40, 0, 40])
    if legend:
        ax.legend(loc='upper right', frameon=False)

    # Total current
    #k = model.labelx('cellular_current').qname()
    #ax2 = ax.twinx()
    #ax2.set_ylim(-0.1, 0.85)
    #ax2.plot(d.time(), d[k], 'r')

    # Contributions
    ax = fig.add_subplot(gr[1:, 0])
    ax.set_xlabel('Time (s)')
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

plot(grid[0, 0], 'nygren', legend=True)
plot(grid[0, 1], 'maleckar', ylabel=None)
plot(grid[0, 2], 'koivumaki', ylabel=None)
plot(grid[1, 0], 'courtemanche')
plot(grid[1, 1], 'ni', ylabel=None)
plot(grid[1, 2], 'aguilar', ylabel=None)
plot(grid[2, 0], 'grandi')
plot(grid[2, 1], 'voigt', ylabel=None)
plot(grid[2, 2], 'ellinwood', ylabel=None)
plot(grid[3, 0], 'bai')

# Legend
ax = fig.add_subplot(grid[3, 1])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
labels = [shared.current_names[x] for x in current_colours]
ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
#ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)

# Show / store
plt.savefig('atrial.png')
plt.savefig('atrial.pdf')
print('Done')
