#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human atrial models.
#
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import myokit
import myokit.lib.plots as mp

import shared


# Update matplotlib styles
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'

# Current colors
cmap = matplotlib.cm.get_cmap('tab20')
current_colours = dict(shared.current_colours)
del(current_colours['I_Kb'])
del(current_colours['I_Kur'])
del(current_colours['I_ClCa'])
del(current_colours['I_Cl,B'])
del(current_colours['I_K,ACh'])
del(current_colours['I_K,ATP'])

# Human atrial models
model_names = {
    'paci-2013': 'paci-2013-ventricular.mmt',
    'paci-2018': 'paci-2018.mmt',
    'paci-2020': 'paci-2020.mmt',
    'kernik': 'kernik-2019.mmt',
}

fancy_names = {
    'paci-2013': 'Paci et al. 2013 (ventricular)',
    'paci-2018': 'Paci et al. 2018',
    'paci-2020': 'Paci et al. 2020',
    'kernik': 'Kernik et al. 2019',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'paci-2013' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_f': 'if.If',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_Na,B': 'ibna.IbNa',
            'I_Ca,B': 'ibca.IbCa',
            'I_Ca,P': 'ipca.IpCa',
            'I_Na': 'ina.INa',
        }
    elif 'paci-2018' in name or 'paci-2020' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_f': 'if.If',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL',
            'I_NaL': 'inal.INaL',
            'I_Na,B': 'ibna.IbNa',
            'I_Ca,B': 'ibca.IbCa',
            'I_Ca,P': 'ipca.IpCa',
            'I_Na': 'ina.INa',
        }
    elif 'kernik-' in name:
        currents = {
            'I_NaCa': 'inaca.i_NaCa',
            'I_to': 'ito.i_to',
            'I_Kr': 'ikr.i_Kr',
            'I_Ks': 'iks.i_Ks',
            'I_f': 'ifunny.i_f',
            'I_K1': 'ik1.i_K1',
            'I_NaK': 'inak.i_NaK',
            'I_CaL': 'ical.i_CaL',
            'I_CaT': 'icat.i_CaT',
            'I_Na,B': 'ibna.i_b_Na',
            'I_Ca,B': 'ibca.i_b_Ca',
            'I_Ca,P': 'ipca.i_PCa',
            'I_Na': 'ina.i_Na',
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
cl = 800
protocol = myokit.pacing.blocktrain(cl, duration=5, offset=50)

# Load and prepare models
models = {}
for name, fname in model_names.items():
    pre_pace = True
    if 'kernik' in name:
        pre_pace = False
    model = myokit.load_model(os.path.join('models', 'hipsc', fname))
    shared.prepare_model(model, protocol, current_variables(model), pre_pace)
    models[name] = model

# Maximum time to show in plots
tmax = 800

def text(ax, x, y, t, c='w'):
    ax.text(x, y, t, color=c, transform=ax.transAxes, fontweight='bold',
            horizontalalignment='right', verticalalignment='center')

# Create figure
fig = plt.figure(figsize=(9, 9))
fig.subplots_adjust(0.075, 0.05, 0.98, 0.97, hspace=0.35, wspace=0.2)
grid = GridSpec(3, 3)

#
# Top row: Paci models
#
# Paci 2013
code = 'paci-2013'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Paci 2018
code = 'paci-2018'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 1])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Paci 2020
code = 'paci-2020'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 2])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

#
# Middle row
#
# Kernik 2019
code = 'kernik'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

#
# Legend
#
ax = fig.add_subplot(grid[1, 2])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
labels = [shared.current_names[x] for x in current_colours]
#ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)


# Show / store
plt.savefig('hipsc.png')
plt.savefig('hipsc.pdf')
print('Done')
