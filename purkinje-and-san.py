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
cmap = matplotlib.colormaps['tab20']
current_colours = dict(shared.current_colours)
del(current_colours['I_Kur'])
del(current_colours['I_ClCa'])
del(current_colours['I_Cl,B'])
del(current_colours['I_K,ACh'])
del(current_colours['I_K,ATP'])

# Human atrial models
model_names = {
    'sampson': 'sampson-2010.mmt',
    'stewart': 'stewart-2009.mmt',
    'trovato': 'trovato-2020.mmt',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'sampson' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito_total',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_f': 'ihcn.IHCN',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_total',
            'I_CaT': 'icat.I',
            'I_Ca,P': 'ipca.IpCa',
            'I_NaL': 'nav11.INa1',
            'I_Na': 'nav15.INa',
        }
    elif 'stewart' in name:
        currents = {
            'I_to': 'ito.i_to_total',
            'I_Kr': 'ikr.i_Kr',
            'I_Ks': 'iks.i_Ks',
            'I_Kb': 'ipk.i_p_K',
            'I_f': 'if.i_f_total',
            'I_K1': 'ik1.i_K1',
            'I_NaK': 'inak.i_NaK',
            'I_CaL': 'ical.i_CaL',
            'I_NaCa': 'inaca.i_NaCa',
            'I_Na,B': 'ibna.i_b_Na',
            'I_Ca,B': 'ibca.i_b_Ca',
            'I_Ca,P': 'ipca.i_p_Ca',
            'I_Na': 'ina.i_Na',
        }
    elif 'trovato' in name:
        currents = {
            'I_to': 'ito.Ito_total',
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_f': 'if.If',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_tot',
            'I_CaT': 'icat.ICaT',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inacass.INaCa_tot',
            'I_Na,B': 'inab.INab',
            'I_Ca,B': 'icab.ICab',
            'I_Ca,P': 'ipca.IpCa',
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
    model = myokit.load_model(os.path.join('models', 'c', fname))
    if 'stewart' in name:
        c = model.get('ito')
        v = c.add_variable('i_to_total')
        v.set_unit(c.get('i_to').unit())
        v.set_rhs('ito.i_to + isus.i_sus')
        c = model.get('if')
        v = c.add_variable('i_f_total')
        v.set_unit(c.get('i_f_Na').unit())
        v.set_rhs('i_f_Na + i_f_K')

    elif 'sampson' in name:
        c = model.get('ito')
        v = c.add_variable('Ito_total')
        v.set_unit(c.get('Ito1').unit())
        v.set_rhs('ito.Ito1 + isus.Isus')
        c = model.get('ical')
        v = c.add_variable('ICaL_total')
        v.set_unit(c.get('ICa').unit())
        v.set_rhs('ICa + ICaK')

    elif 'trovato' in name:
        c = model.get('ito')
        v = c.add_variable('Ito_total')
        v.set_unit(c.get('Ito').unit())
        v.set_rhs('ito.Ito + isus.Isus')

    pre_pace = True
    if 'stewart' in name:
        # 2024-09-03 Stewart model destabilises when pre-paced
        pre_pace = False

    shared.prepare_model(model, protocol, current_variables(model), pre_pace)
    models[name] = model
print('Finished preparation.\nPreparing plots')


# Maximum time to show in plots
tmax = 800


def text(ax, x, y, t, c='w'):
    ax.text(x, y, t, color=c, transform=ax.transAxes, fontweight='bold',
            horizontalalignment='right', verticalalignment='center')


def plot(code, grid, i, j, d, ylabel='Relative contribution'):
    gr = grid[i, j].subgridspec(4, 1, hspace=0)

    # V and CaT
    ax = fig.add_subplot(gr[0, 0])
    ax.set_title(model.meta['display_name'])
    ax.set_xticklabels([])
    ax.plot(d.time(), d['membrane.V'], 'k')
    ax.set_xlim(0, tmax)
    ax.set_ylim(-95, 45)
    ax.set_yticks([-80, -40, 0, 40])
    #ax.set_yticklabels([None, -40, 0, 40])

    #ax = ax.secondary_yaxis()

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
grid = GridSpec(4, 3)

#
# Top row: Purkinje
#
# Stewart 2009
code = 'stewart'
if code in models:
    model = models[code]
    currents, colours = current_variables(model, True)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    d = s.run(tmax)
    plot(code, grid, 0, 0, d)

# Sampson 2010
code = 'sampson'
if code in models:
    model = models[code]
    currents, colours = current_variables(model, True)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    d = s.run(tmax)
    plot(code, grid, 0, 1, d, ylabel=None)

# Trovato 2020
code = 'trovato'
if code in models:
    model = models[code]
    currents, colours = current_variables(model, True)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    d = s.run(tmax)
    plot(code, grid, 0, 2, d, ylabel=None)

#
# Third row: SAN
#

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
plt.savefig('purkinje-and-san.png')
plt.savefig('purkinje-and-san.pdf')
print('Done')
